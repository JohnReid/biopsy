/* Copyright John Reid 2007
*/

#include "bio-pch.h"




#include "bio/application.h"
#include "bio/remo.h"
#include "bio/remo_analysis.h"
#include "bio/pair_analysis.h"
#include "bio/biobase_binding_model.h"
#include "bio/biobase_score.h"
USING_BIO_NS

#include <boost/foreach.hpp>
#include <boost/regex.hpp>
namespace po = boost::program_options;

#include <iostream>
#include <fstream>
#include <string>
using namespace std;




struct BinderRegexFilter
{
	const boost::regex r;
	BinderRegexFilter( const std::string & pattern )
		: r( pattern )
	{
	}

	bool
	operator()( const BindingModel::hit_t & hit )
	{
		return regex_search( hit.get_binder()->get_name(), r );
	}

	bool
	operator()( const Matrix::map_t::value_type & v )
	{
		return regex_search( BiobaseDb::singleton().get_entry< MATRIX_DATA >( v.first )->get_name(), r );
	}

	bool
	operator()( const Site::map_t::value_type & v )
	{
		Site * site = BiobaseDb::singleton().get_entry< SITE_DATA >( v.first );

		//make sure it is a consensus sequence
		return "CONS" == site->id.factor && regex_search( site->get_name(), r );
	}
};



struct PairAnalysisApp : Application, AnalysisVisitor
{
	unsigned num_output;
	PairStatistics< BindingModel > pair_statistics;
	std::string binder_regex_pattern;
	std::string input;
	std::string output;
	bool analyse_distances;
	bool analyse_singles;
	bool analyse_clusters;
	boost::scoped_ptr< BinderRegexFilter > binder_filter;

	PairAnalysisApp()
	{
		get_options().add_options()
			( "input,i", po::value( &input )->default_value( "" ), "serialised statistics" )
			( "regex,r", po::value( &binder_regex_pattern )->default_value( "." ), "regex to match binder names" )
			( "distance,d", po::value( &pair_statistics.distance) ->default_value( 15 ), "distance apart pairs can be" )
			( "num_output,n", po::value( &num_output )->default_value( 16 ), "how many pairs to output" )
			( "output,o", po::value( &output )->default_value( "" ), "where to serialise statistics to" )
			( "distances", po::value( &analyse_distances )->default_value( true ), "analyse distances between pairs" )
			( "clusters", po::value( &analyse_clusters )->default_value( true ), "analyse clusters" )
			( "singles", po::value( &analyse_singles )->default_value( true ), "analyse singles" )
			;

		add_analysis_options(get_options());
	}

	bool visit_sequence_group(const std::string & seq_group_name)
	{
		return true;
	}

	void visit_remo(
		const std::string & seq_group_name,
		ReMoLocation location,
		const ReMoRange & range,
		const std::string & remo_name,
		bifa_hits_t & results,
		const seq_t & sequence)
	{
		typedef bifa_hits_t::index< BindingHitSet< BindingModel >::position >::type by_position_t;
		const by_position_t & by_position = results.get< BindingHitSet< BindingModel >::position >();
		pair_statistics.add_hits( 
			seq_group_name, 
			make_filter_iterator( *binder_filter, by_position.begin(), by_position.end() ), 
			make_filter_iterator( *binder_filter, by_position.end(), by_position.end() ), 
			sequence.length() );
	}

	int task()
	{
		cout << "Finding doubles not further than " << pair_statistics.distance << " apart\n";
		cout << "Outputting " << num_output << " most significant doubles\n";

		if( "." != binder_regex_pattern )
		{
			cout << "Using filter on binders: \"" << binder_regex_pattern << "\"\n";
		}
		binder_filter.reset( new BinderRegexFilter( binder_regex_pattern ) );
		//check we match some pssms
		{
			BindingModel::set_t model_universe;
			transform_biobase_sites_and_matrices(
				*binder_filter,
				Link2BiobaseBindingModel( BioEnvironment::singleton().get_tf_binding_prior(), false ),
				std::inserter( model_universe, model_universe.begin() ) );
			if( model_universe.empty() )
			{
				throw std::logic_error( "Regex pattern does not match any PSSMs" );
			}
			cout << "PSSMS that match filter:\n";
			BOOST_FOREACH( BindingModel * model, model_universe )
			{
				cout << model->get_name() << "\n";
			}
		}


		//deserialise analysis and examine unless have already done so
		if ("" == input)
		{
			deserialise_analysis();

			visit_remo_analysis();
			pair_statistics.num_remos = num_remos;
			pair_statistics.num_bases = num_bases;
		}
		else
		{
			cout << "Deserialising statistics from \"" << input << "\"\n";
			std::ifstream stream(input.c_str(), std::ios::binary);

			boost::archive::binary_iarchive(stream) >> pair_statistics;
		}

		if( analyse_distances )
		{
			cout << "Analysing distances\n";
			pair_statistics.analyse_distances( num_output );
		}

		if( analyse_singles )
		{
			cout << "Analysing singles\n";
			pair_statistics.analyse_singles( num_output );
		}

		if( analyse_clusters )
		{
			cout << "Analysing clusters\n";
			pair_statistics.analyse_clusters( num_output );
		}

		if ("" != output)
		{
			cout << "Serialising statistics to \"" << output << "\"\n";
			std::ofstream stream(output.c_str(), std::ios::binary);

			boost::archive::binary_oarchive(stream) << const_cast< const PairStatistics< BindingModel > & >(pair_statistics);
		}

		return 0;
	}
};

int
main(int argc, char * argv[])
{
	return PairAnalysisApp().main(argc, argv);
}

