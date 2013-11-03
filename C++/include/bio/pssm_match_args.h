
#ifndef BIO_PSSM_MATCH_ARGS_H_
#define BIO_PSSM_MATCH_ARGS_H_

#include "bio/defs.h"
#include "bio/sequence.h"
#include "bio/biobase_filter.h"
#include "bio/run_match.h"




BIO_NS_START




/** Arguments for pssm matching algorithm. */
struct PssmMatchArgs
{
	float_t threshold;
	std::string species_filter;
	bool use_consensus_sequences;
	std::string matrix_regex;
	bool use_bayesian;
	bool use_or_better;

	PssmMatchArgs(
		float_t threshold = 0.05f,
		bool use_consensus_sequences = true,
		std::string species_filter = "V",
		std::string matrix_regex = ".",
		bool use_bayesian = true,
		bool use_or_better = false);

	void add_options(boost::program_options::options_description & options);

	BiobasePssmFilter get_filter() const;

	/** Run the algorithm with the arguments. */
	template< typename It, typename OutIt >
	void
	run_pssm_match(
		It seq_begin,
		It seq_end,
		OutIt output)
	{
		const BiobasePssmFilter filter = get_filter();

		pssm_match(
			seq_begin,
			seq_end,
			threshold,
			use_bayesian ? BAYESIAN_SCORE_ALGORITHM : OTT_SCORE_ALGORITHM,
			use_or_better,
			filter,
			filter,
			output);
	}
};

/** Display the arguments on the stream. */
std::ostream &
operator<<(std::ostream & os, const PssmMatchArgs & args);



BIO_NS_END

#endif //BIO_PSSM_MATCH_ARGS_H_

