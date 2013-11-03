
#ifndef BIO_ENVIRONMENT_H_
#define BIO_ENVIRONMENT_H_

#include "bio/defs.h"
#include "bio/singleton.h"

#include <string>

#ifdef WIN32
# define DIR_SEP "\\"
#else //WIN32
# define DIR_SEP "/"
#endif //WIN32

BIO_NS_START

struct BioEnvironment
	: Singleton< BioEnvironment >
{
	BioEnvironment();

	/** The name of the file we write logging output to. Defaults to "" and then we use std::cout. */
	std::string get_log_file_name() const;

	/** The stream we write logging output to. Defaults to std::cout. */
	std::ostream * get_log_stream() const;

	/** The URI for the KEGG web service. */
	std::string get_kegg_wsdl_uri() const;

	/** The prior for TF binding probability. */
	float_t get_tf_binding_prior() const;

	/** The directory where the serialised files are stored. */
	std::string get_serialised_dir() const;

	/** The directory where the ensembl files are stored. */
	std::string get_ensembl_dir() const;

	/** The directory where the chromosome files are stored. */
	std::string get_chromosome_dir() const;

	/** The directory where the biobase files are stored. */
	std::string get_biobase_dir() const;

	/** The directory where the transfac files are stored. */
	std::string get_transfac_dir() const;

	/** The directory where the transcompel files are stored. */
	std::string get_transcompel_dir() const;

	/** The directory where the transpath files are stored. */
	std::string get_transpath_dir() const;
	
	/** The file where the site test cases are stored. */
	std::string get_site_test_cases_file() const;
	
	/** The file where the default remo archive is stored. */
	std::string get_default_remo_archive_file() const;
	
	/** The file where the matrix match parameters are stored. */
	std::string get_matrix_match_file() const;

	/** The file where the min false positive thresholds are stored. */
	std::string get_matrix_min_fp_file() const;

	/** The file where the min false negative thresholds are stored. */
	std::string get_matrix_min_fn_file() const;

	/** The file where the min sum thresholds are stored. */
	std::string get_matrix_min_sum_file() const;

	/** The file where the likelihoods are stored. */
	std::string get_likelihoods_cache_file() const;

	/** The file where the pssm likelihoods are stored. */
	std::string get_pssm_likelihoods_cache_file() const;

	/** The file where the HMMs for different species are stored. */
	std::string get_species_hmm_file() const;

	/** The file where the script for the output svg is stored. */
	std::string get_svg_script_file() const;

	/** The file where the script for version 2 of the output svg is stored. */
	std::string get_svg_script_file_ver_2() const;

	/** The file where the tss estimates are stored. */
	std::string get_serialised_tss_estimates_file() const;

	/** The file where the new tss format is stored. */
	std::string get_tss_file_new_format() const;

	/** The file where the tss file is stored. */
	std::string get_tss_file() const;

	/** The file where the tss clones are stored. */
	std::string get_tss_clones_file() const;

	/** The file where the factor synonyms is stored. */
	std::string get_factor_synonyms_file() const;

	/** The directory where the custom PSSM files are stored. */
	std::string get_custom_pssm_dir() const;

	/** The prior for TF binding probability. */
	float_t tf_binding_prior;

	/** Where the data is stored. */
	std::string data_dir;
	
	/** Which version of TRANSPATH are we using. */
	size_t transpath_major_version;
	size_t transpath_minor_version;
	size_t transcompel_major_version;
	size_t transcompel_minor_version;
	size_t transfac_major_version;
	size_t transfac_minor_version;
	std::string custom_PSSM_version;
	
	/** How to link to biobase. */
	std::string biobase_url_prefix;

	/** The names of the species we use. */
	std::vector<std::string> species_prefixes;
	
	/** The port we run the embedded web server on. */
	size_t http_port;

	/** The number of quanta we use for normalisation. */
	size_t num_normalisation_quanta;

	/** The number of distinct scores we allow before approximating when calculating pssm likelihoods. */
	size_t max_pssm_likelihood_map_size;

	/** The number of Hidden Markov Model states we use by default. */
	size_t num_hmm_states;

	/** The log output stream. */
	mutable boost::scoped_ptr< std::ofstream > log_os;

	/** The name of the file we write the logs to. "" is the default, then we use std::cout. */
	std::string log_file_name;
};



BIO_NS_END


#endif //BIO_ENVIRONMENT_H_
