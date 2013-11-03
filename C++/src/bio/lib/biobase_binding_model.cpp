/* Copyright John Reid 2007
*/

#include "bio-pch.h"


#include "bio/defs.h"

#include "bio/pssm_bayesian_binding_model.h"
#include "bio/bayesian_binding_model.h"
#include "bio/biobase_likelihoods.h"
#include "bio/biobase_binding_model.h"
#include "bio/match_binding_model.h"
#include "bio/biobase_db.h"
#include "bio/biobase_match.h"
#include "bio/pssm_cache.h"
#include "bio/cache.h"
#include "bio/singleton.h"
#include "bio/matrix_match.h"

BIO_NS_START




BiobaseBindingModel::BiobaseBindingModel( const BiobaseBindingModel::parameter_t & parameters )
	: BayesianBindingModel< PssmScorer, QuantisedScores, QuantisedScores >(
		BiobaseDb::singleton().get_pssm_entry( parameters.link )->get_name(),
		PssmScorer( PssmCache::singleton()( parameters.link ) ),
		get_biobase_quantised_scores(
			parameters.link,
			true,
			parameters.or_better ),
		get_biobase_quantised_scores(
			parameters.link,
			false,
			parameters.or_better ),
		parameters.p_Hb_prior )
	, parameters( parameters )
{
}



const BindingModel::parameter_t *
BiobaseBindingModel::get_parameters() const
{
	return boost::addressof( parameters );
}




/**
Creates binding models from biobase tablelink references.
*/
struct BiobaseBindingModelCreator
	: std::unary_function< BiobaseBindingModel::parameter_t, BindingModel::ptr_t >
{
	BindingModel::ptr_t operator()( const BiobaseBindingModel::parameter_t & parameters ) const
	{
		return BindingModel::ptr_t( new BiobaseBindingModel( parameters ) );
	}
};




/**
Caches biobase binding models.
*/
struct BiobaseBindingModelCache
	: Cache< BiobaseBindingModelCreator >
	, Singleton< BiobaseBindingModelCache >
{
};





BindingModel *
BiobaseBindingModel::parameter_t::get_model() const
{
	return BiobaseBindingModelCache::singleton()( *this ).get();
}

BiobaseBindingModel::parameter_t::parameter_t(
	const TableLink & link,
	double p_Hb_prior,
	bool or_better )
	: link( link )
	, p_Hb_prior( p_Hb_prior )
	, or_better( or_better )
{
}

bool
BiobaseBindingModel::parameter_t::operator<( const BiobaseBindingModel::parameter_t & rhs ) const
{
	if( link < rhs.link ) return true;
	else if( ! ( rhs.link < link ) ) {
		if( p_Hb_prior < rhs.p_Hb_prior ) return true;
		else if( ! ( rhs.p_Hb_prior < p_Hb_prior ) ) {
			return or_better < rhs.or_better;
		}
	}
	return false;
}



MatchBindingModel::MatchBindingModel( const MatchBindingModel::parameter_t & parameters )
: parameters( parameters )
, pssm( make_pssm( parameters.link ) )
{

	MatrixMatch::map_t::const_iterator mm = get_min_fp_match_map().find( parameters.link );
	if( get_min_fp_match_map().end() == mm )
	{
		throw std::logic_error( BIO_MAKE_STRING( "Could not find threshold for: " << parameters.link ) );
	}
	threshold = mm->second.threshold;
}


MatchBindingModel::~MatchBindingModel()
{
}


std::string
MatchBindingModel::get_name() const
{
	return
		BIO_MAKE_STRING(
			"Biobase TRANSFAC Match algorithm: "
			<< parameters.link );
}

unsigned
MatchBindingModel::get_num_bases() const
{
	return pssm.size();
}

const MatchBindingModel::parameter_t *
MatchBindingModel::get_parameters() const
{
	return &parameters;
}

double
MatchBindingModel::operator()(
	seq_t::const_iterator begin,
	bool match_complement,
	BindingModelContext * context) const
{
	const double score = pssm.score( begin, match_complement );
	const bool above = score > threshold;
	return
		above
			? 1.0
			: 0.0 ;
}




/**
Creates match binding models from biobase tablelink references.
*/
struct MatchBindingModelCreator
	: std::unary_function< MatchBindingModel::parameter_t, BindingModel::ptr_t >
{
	BindingModel::ptr_t operator()( const MatchBindingModel::parameter_t & parameters ) const
	{
		return BindingModel::ptr_t( new MatchBindingModel( parameters ) );
	}
};




/**
Caches match binding models.
*/
struct MatchBindingModelCache
	: Cache< MatchBindingModelCreator >
	, Singleton< MatchBindingModelCache >
{
};




MatchBindingModel::parameter_t::parameter_t(
	TableLink link )
	: link( link )
{
}

MatchBindingModel::parameter_t::~parameter_t()
{
}

BindingModel *
MatchBindingModel::parameter_t::get_model() const
{
	return MatchBindingModelCache::singleton()( *this ).get();
}

bool
MatchBindingModel::parameter_t::operator<( const MatchBindingModel::parameter_t & rhs ) const
{
	return link < rhs.link;
}

BIO_NS_END

BOOST_CLASS_EXPORT( BIO_NS::MatchBindingModel::parameter_t )
BOOST_CLASS_EXPORT( BIO_NS::BiobaseBindingModel::parameter_t )

