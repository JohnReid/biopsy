/* Copyright John Reid 2007
*/

#include "bio-pch.h"


#include "bio/defs.h"

#include "bio/bifa_algorithm.h"
#include "bio/bifa_test_case.h"
#include "bio/run_match.h"
#include "bio/adjust_hits.h"
#include "bio/biobase_binding_model.h"
#include "bio/match_binding_model.h"
#include "bio/biobase_score.h"

#include <boost/iterator/filter_iterator.hpp>
#include <boost/iterator/transform_iterator.hpp>

#include <gsl/gsl_sf_pow_int.h>

#include <iostream>



BIO_NS_START



BiFaInput::BiFaInput( const seq_t & centre_sequence )
: centre_sequence( centre_sequence )
{
}


BiFaAlgorithm::~BiFaAlgorithm() { }


BiFaAlgorithm::ptr_t
BiFaAlgorithm::get_default_algorithm()
{
	static BiFaAlgorithm::ptr_t default_algorithm;

	if( ! default_algorithm )
	{
		default_algorithm.reset( new TransfacBiFaAlgorithm( PssmMatchArgs(), true ) );
	}

	return default_algorithm;
}

TransfacBiFaAlgorithm::TransfacBiFaAlgorithm(
	const PssmMatchArgs & args,
	bool adjust_phylo,
	bool verbose)
: verbose( verbose )
, adjust_phylo( adjust_phylo )
, args( args )
{
}


TransfacBiFaAlgorithm::~TransfacBiFaAlgorithm() { }


std::string
TransfacBiFaAlgorithm::get_name() const
{
	return
		BIO_MAKE_STRING(
			"BiFa-"
			<< ( adjust_phylo ? "phylo" : "no phylo" )
			<< "-"
			<< ( args.use_bayesian ? "bayesian" : "not bayesian" )
			<< "-"
			<< ( args.use_or_better ? "or better" : "equal" )
			<< "-"
			<< args.threshold
			);
}


BindingModel *
TransfacBiFaAlgorithm::get_model_for( const boost::any & key )
{
	return Link2BiobaseBindingModel( BioEnvironment::singleton().get_tf_binding_prior(), args.use_or_better )( boost::any_cast< TableLink >( key ) );
}

void
TransfacBiFaAlgorithm::fill_model_universe( BindingModel::set_t & universe )
{
	transform_biobase_sites_and_matrices(
		args.get_filter(),
		Link2BiobaseBindingModel( BioEnvironment::singleton().get_tf_binding_prior(), args.use_or_better ),
		std::inserter( universe, universe.begin() ) );
}


BiFaOutput::ptr_t
TransfacBiFaAlgorithm::operator()(const BiFaInput & input)
{
	BiFaOutput::ptr_t result(new BiFaOutput);

	//score the sites and matrices
	score_all_biobase_pssms(
		make_sequence_scorer(
			input.centre_sequence.begin(),
			input.centre_sequence.end(),
			args.threshold,
			std::inserter( result->hits, result->hits.begin() )
		),
		args.get_filter(),
		Link2BiobaseBindingModel( BioEnvironment::singleton().get_tf_binding_prior(), args.use_or_better )
	);

	if ( verbose )
	{
		std::cout << result->hits.size() << " matches over threshold\n";
	}

	//are we making adjustments for the phylogenetic sequences
	if ( adjust_phylo )
	{
		//raise the threshold to the power of the number of sequences
		const float_t phylo_threshold = float_t( gsl_sf_pow_int( args.threshold, input.conserved_sequences.size() + 1 ) );

		adjust_hits(
			result->hits,
			input.conserved_sequences,
			phylo_threshold);

		if ( verbose )
		{
			std::cout
				<< "Adjusted for " << input.conserved_sequences.size()
				<< " sequences, still have "
				<< std::count_if(
					result->hits.begin(),
					result->hits.end(),
					boost::bind(
						std::greater< double >(),
						boost::bind(
							&BindingModel::hit_t::get_p_binding,
							_1 ),
						phylo_threshold ) )
				<< " hits above the threshold\n";
		}
	}

	return result;
}



std::string
MatchBiFaAlgorithm::get_name() const
{
	return "Match";
}


struct MatchMapFilter
{
	bool operator()( const TableLink & link ) const
	{
		return get_min_fp_match_map().find( link ) != get_min_fp_match_map().end();
	}

	bool operator()( Matrix::map_t::value_type matrix ) const
	{
		return ( *this )( matrix.first );
	}

	bool operator()( Site::map_t::value_type site ) const
	{
		return false;
	}
};

BindingModel *
MatchBiFaAlgorithm::get_model_for( const boost::any & key )
{
	const TableLink & link = boost::any_cast< const TableLink & >( key );

	return
		MatchMapFilter()( link )
			? Link2MatchBindingModel()( link )
			: 0
			;
}

void
MatchBiFaAlgorithm::fill_model_universe( BindingModel::set_t & universe )
{
	transform_biobase_sites_and_matrices(
		MatchMapFilter(),
		Link2MatchBindingModel(),
		std::inserter( universe, universe.begin() ) );
}


BiFaOutput::ptr_t
MatchBiFaAlgorithm::operator()(const BiFaInput & input)
{
	BiFaOutput::ptr_t result(new BiFaOutput);

	//score the sites and matrices
	score_all_biobase_pssms(
		make_sequence_scorer(
			input.centre_sequence.begin(),
			input.centre_sequence.end(),
			0.5,					// our hits for match algorithm are 0 or 1 so 0.5 is a good enough threshold
			std::inserter( result->hits, result->hits.begin() )
		),
		MatchMapFilter(),
		Link2MatchBindingModel()
	);

	return result;
}





ROCPoint::ROCPoint( double specificity, double sensitivity )
	: specificity( specificity )
	, sensitivity( sensitivity )
{
}


std::ostream &
operator<<( std::ostream & os, const ROCPoint & roc_point )
{
	os << "(" << roc_point.specificity << "," << roc_point.sensitivity << ")";
	return os;
}

BinaryTestResults::BinaryTestResults(
	unsigned true_positives,
	unsigned false_positives,
	unsigned true_negatives,
	unsigned false_negatives)
	: true_positives( true_positives )
	, false_positives( false_positives )
	, true_negatives( true_negatives )
	, false_negatives( false_negatives )
{
}



BinaryTestResults &
BinaryTestResults::operator+=( const BinaryTestResults & rhs )
{
	true_positives += rhs.true_positives;
	false_positives += rhs.false_positives;
	true_negatives += rhs.true_negatives;
	false_negatives += rhs.false_negatives;

	return *this;
}



BinaryTestResults &
BinaryTestResults::operator()(bool tested_positive, bool should_be_positive)
{
	if (tested_positive)
	{
		if (should_be_positive)
		{
			++true_positives;
		}
		else
		{
			++false_positives;
		}
	}
	else
	{
		if (should_be_positive)
		{
			++false_negatives;
		}
		else
		{
			++true_negatives;
		}
	}

	return *this;
}


ROCPoint
BinaryTestResults::get_roc_point() const
{
	return
		ROCPoint(
			double( true_negatives ) / double( true_negatives + false_positives ),
			double( true_positives ) / double( true_positives + false_negatives ) );
}



double
BiFaTestCase::calculate_next_threshold(
	BiFaTestCase::threshold_result_map_t::const_iterator begin,
	BiFaTestCase::threshold_result_map_t::const_iterator end,
	double min_threshold,
	double max_threshold)
{
	//either we use the specificity or sensitivity to choose the next threshold
	static bool use_specificity = ! use_specificity;

	//check args
	if (min_threshold > max_threshold)
	{
		throw std::invalid_argument( "Min must be <= max threshold" );
	}

	//if we have no results try the average
	if ( end == begin )
	{
		return (min_threshold + max_threshold) / 2.0;
	}

	//look for the largest gap between successive iterators
	double largest_gap = 0.0;
	threshold_result_map_t::const_iterator max_gap = end;
	for (threshold_result_map_t::const_iterator left = begin;
		end != left;
		++left)
	{
		threshold_result_map_t::const_iterator right = left;
		++right;
		if (end == right)
		{
			break;
		}

		//is there a bigger gap between these two than any two before?
		const double left_value = use_specificity ? left->second.specificity : left->second.sensitivity;
		const double right_value = use_specificity ? right->second.specificity : right->second.sensitivity;
		const double gap = fabs( left_value - right_value );
		if ( gap > largest_gap )
		{
			largest_gap = gap;
			max_gap = left;
		}
	}

	//check the gaps between the ends of the sequence and the min/max thresholds

	return (min_threshold + max_threshold) / 2.0;
}


BIO_NS_END
