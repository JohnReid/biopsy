/* Copyright John Reid 2007
*/

#include "bio-pch.h"


#include "bio/defs.h"
#include "bio/pssm_match_args.h"


#include <boost/program_options.hpp>
namespace po = boost::program_options;



BIO_NS_START




PssmMatchArgs::PssmMatchArgs(
	float_t threshold,
	bool use_consensus_sequences,
	std::string species_filter,
	std::string matrix_regex,
	bool use_bayesian,
	bool use_or_better)
: threshold( threshold )
, species_filter( species_filter )
, use_consensus_sequences( use_consensus_sequences )
, matrix_regex( matrix_regex )
, use_bayesian( use_bayesian )
, use_or_better( use_or_better )
{
}

BiobasePssmFilter
PssmMatchArgs::get_filter() const
{
	return
		BiobasePssmFilter(
			use_consensus_sequences,
			species_filter,
			matrix_regex );
}


void
PssmMatchArgs::add_options(boost::program_options::options_description & options)
{
	options.add_options()
		("threshold,t", po::value(&threshold), "threshold for hits")
		("species_filter", po::value(&species_filter), "species filter")
		("consensus", po::value(&use_consensus_sequences), "use consensus_sequences")
		("matrix_regex,r", po::value(&matrix_regex), "RegEx to match matrix names")
		("use_bayesian", po::value(&use_bayesian), "use Bayesian scoring algorithm")
		;
}


std::ostream &
operator<<(std::ostream & os, const PssmMatchArgs & args)
{
	os
		<< "Threshold: " << args.threshold << "\n"
		<< "Consensus: " << (args.use_consensus_sequences ? "yes" : "no") << "\n"
		<< "Bayesian: " << (args.use_bayesian ? "yes" : "no") << "\n"
		<< "Or better: " << (args.use_or_better ? "yes" : "no") << "\n"
		<< "Species filter: " << args.species_filter << "\n"
		<< "Matrix regex: " << args.matrix_regex << "\n";

	return os;
}

BIO_NS_END

