/**
@file

Copyright John Reid 2006

*/

#include "biopsy/python.h"
//#include <boost/python.hpp>
#include "bio/pssm_match_args.h"
#include "bio/svg_match.h"

using namespace boost;
using namespace boost::python;
using namespace std;


namespace biopsy {

BIO_NS::PssmMatchArgs & nullPssmMatchArgs()
{
	static BIO_NS::PssmMatchArgs pma;
	return pma;
}

BIO_NS::BuildSvgArgs & nullBuildSvgArgs()
{
	static BIO_NS::BuildSvgArgs bsa;
	return bsa;
}

void
pssm_match(
	const sequence & sequence,
	const sequence_vec & phylo_seqs,
	BIO_NS::PssmMatchArgs & pssm_match_args = nullPssmMatchArgs(),
	BIO_NS::BuildSvgArgs & build_svg_args = nullBuildSvgArgs(),
	BIO_NS::float_t phylo_threshold = 0.03f
)
{
	USING_BIO_NS;
	
	match_result_vec_t results;
	pssm_match_args.run_pssm_match(
		sequence.begin(),
		sequence.end(),
		std::inserter( results, results.begin() ) );

	adjust_hits_for_phylo_sequences(
		results.begin(),
		results.end(),
		phylo_seqs,
		phylo_threshold,
		pssm_match_args.use_bayesian ? BAYESIAN_SCORE_ALGORITHM : OTT_SCORE_ALGORITHM );

	build_svg_args.min_threshold = pssm_match_args.threshold;
	build_svg_args.run_build_svg(
		sequence,
		results );
}


void export_pssm_match()
{
	using boost::python::arg;

	/** These are here to mimic pssm_match functionality. */
	USING_BIO_NS;
	class_< 
		PssmMatchArgs
	>( 
		"PssmMatchArgs",
		"Parameters for old style BiFa analysis (pssm_match)",
		init<
			BIO_NS::float_t,
			bool,
			std::string,
			std::string,
			bool,
			bool
		>(
			( 
				arg( "threshold" ) = 0.05f, 
				arg( "use_consensus_sequences" ) = true,
				arg( "species" ) = "V",
				arg( "matrix_regex" ) = ".",
				arg( "bayesian" ) = true,
				arg( "or_better" ) = false
			),
			"Initialise parameters for running a pssm_match"
		)
	)
		.def_readwrite(
			"threshold",
			&PssmMatchArgs::threshold,
			"Threshold for centre sequence" )
		.def_readwrite(
			"use_consensus_sequences",
			&PssmMatchArgs::use_consensus_sequences,
			"Use consensus sequences" )
		.def_readwrite(
			"species_filter",
			&PssmMatchArgs::species_filter,
			"Filter for species" )
		.def_readwrite(
			"species_filter",
			&PssmMatchArgs::species_filter,
			"Filter for species" )
		.def_readwrite(
			"matrix_regex",
			&PssmMatchArgs::matrix_regex,
			"Filter (regex) for matrix names" )
		.def_readwrite(
			"species_filter",
			&PssmMatchArgs::threshold,
			"Filter for species" )
		.def_readwrite(
			"use_bayesian",
			&PssmMatchArgs::use_bayesian,
			"Use Bayesian analysis" )
		.def_readwrite(
			"use_or_better",
			&PssmMatchArgs::use_or_better,
			"Use \"or better\" when calculating binding/not binding log odds ratio" )
	;	

	def(
		"pssm_match",
		pssm_match,
		( 
			arg( "centre_sequence" ), 
			arg( "phylo_sequences" ),
			arg( "pssm_match_args" ) = PssmMatchArgs(),
			arg( "build_svg_args" ) = BuildSvgArgs(),
			arg( "phylo_threshold" ) = 0.03f
		),
		"Old style BiFa analysis" 
	);

}



} //namespace biopsy
