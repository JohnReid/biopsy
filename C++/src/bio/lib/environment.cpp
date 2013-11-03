/* Copyright John Reid 2007
*/

#include "bio-pch.h"


#include "bio/defs.h"


#include "bio/environment.h"

#include <sstream>

BIO_NS_START


BioEnvironment::BioEnvironment()
: tf_binding_prior(float_t(4.0 / (100.0 * 2 * 1000.0))) //100 bases, 1000 TFs, 4 hits expected in conserved region
, data_dir(".")
, transpath_major_version(6)
, transpath_minor_version(2)
, transcompel_major_version(9)
, transcompel_minor_version(3)
, transfac_major_version(9)
, transfac_minor_version(3)
, custom_PSSM_version("X")
, biobase_url_prefix("https://portal.biobase-international.com/cgi-bin/build_t/idb/1.0/")
//, biobase_url_prefix("http://www.biobase-international.com/cgi-bin/biobase/")
, http_port(8080)
, num_normalisation_quanta(100)
, max_pssm_likelihood_map_size(50000)
, num_hmm_states(16)
, log_file_name("")
{
}

std::ostream *
BioEnvironment::get_log_stream() const
{
	//do we need to open the stream
	if( "" != get_log_file_name() && 0 == log_os )
	{
		//yes
		using namespace boost::filesystem;
		path log_file( get_log_file_name() );
		std::cout << "Opening log file: " << log_file._BOOST_FS_NATIVE() << std::endl;
		log_os.reset( new ofstream( log_file ) );
	}

	return
		log_os
			? log_os.get()
			: boost::addressof( std::cout );
			// : std::cout;
}

std::string
BioEnvironment::get_log_file_name() const
{
	return log_file_name;
}

/** The prior for TF binding probability. */
float_t
BioEnvironment::get_tf_binding_prior() const
{
	return tf_binding_prior;
}

/** The URI for the KEGG web service. */
std::string
BioEnvironment::get_kegg_wsdl_uri() const
{
	return "c:/data/Kegg/KEGG.wsdl";
}

/** The directory where the serialised files are stored. */
std::string
BioEnvironment::get_serialised_dir() const
{
	return data_dir + DIR_SEP "serialised";
}

/** The directory where the ensembl files are stored. */
std::string
BioEnvironment::get_ensembl_dir() const
{
	return data_dir + DIR_SEP "ensembl";
}

/** The directory where the chromosome files are stored. */
std::string
BioEnvironment::get_chromosome_dir() const
{
	return get_ensembl_dir() + DIR_SEP "chromosomes";
}

/** The directory where the biobase files are stored. */
std::string
BioEnvironment::get_biobase_dir() const
{
	return data_dir + DIR_SEP "biobase";
}

/** The directory where the custom PSSM files are stored. */
std::string
BioEnvironment::get_custom_pssm_dir() const
{
	return data_dir + DIR_SEP "custom-pssms";
}

/** The directory where the transfac files are stored. */
std::string
BioEnvironment::get_transfac_dir() const
{
	return get_biobase_dir() + DIR_SEP "transfac";
}

/** The directory where the transcompel files are stored. */
std::string
BioEnvironment::get_transcompel_dir() const
{
	return get_biobase_dir() + DIR_SEP "transcompel";
}

/** The directory where the transpath files are stored. */
std::string
BioEnvironment::get_transpath_dir() const
{
	return get_biobase_dir() + DIR_SEP "transpath";
}

std::string
BioEnvironment::get_site_test_cases_file() const
{
	return get_serialised_dir() + DIR_SEP "site_test_cases.bin";
}

std::string
BioEnvironment::get_default_remo_archive_file() const
{
	return data_dir + DIR_SEP "ReMos" DIR_SEP "100" DIR_SEP "100.filtered";
}

std::string
BioEnvironment::get_matrix_match_file() const
{
	std::stringstream str;
	str
		<< get_transfac_dir()
		<< DIR_SEP "matrixTFP"
		<< transfac_major_version
		<< transfac_minor_version
		<< ".lib";

	return str.str();
}

std::string
BioEnvironment::get_matrix_min_fp_file() const
{
	std::stringstream str;
	str
		<< get_transfac_dir()
		<< DIR_SEP "minFP"
		<< transfac_major_version
		<< transfac_minor_version
		<< ".prf";

	return str.str();
}

std::string
BioEnvironment::get_matrix_min_fn_file() const
{
	std::stringstream str;
	str
		<< get_transfac_dir()
		<< DIR_SEP "minFN"
		<< transfac_major_version
		<< transfac_minor_version
		<< ".prf";

	return str.str();
}

std::string
BioEnvironment::get_matrix_min_sum_file() const
{
	std::stringstream str;
	str
		<< get_transfac_dir()
		<< DIR_SEP "minSUM"
		<< transfac_major_version
		<< transfac_minor_version
		<< ".prf";

	return str.str();
}

std::string
BioEnvironment::get_likelihoods_cache_file() const
{
	return get_serialised_dir() + DIR_SEP "likelihoods.txt";
}

std::string
BioEnvironment::get_pssm_likelihoods_cache_file() const
{
	return get_serialised_dir() + DIR_SEP "pssm_likelihoods.txt";
}

/** The file where the HMMs for different species are stored. */
std::string
BioEnvironment::get_species_hmm_file() const
{
	return get_serialised_dir() + DIR_SEP "species_hmm.txt";
}

std::string
BioEnvironment::get_svg_script_file() const
{
	return data_dir + DIR_SEP "scripts" DIR_SEP "bifa.js";
}

std::string
BioEnvironment::get_svg_script_file_ver_2() const
{
	return data_dir + DIR_SEP "scripts" DIR_SEP "bifa_ver_2.js";
}

std::string
BioEnvironment::get_serialised_tss_estimates_file() const
{
	return get_serialised_dir() + DIR_SEP "tss_estimates.bin";
}

std::string
BioEnvironment::get_tss_file_new_format() const
{
	return data_dir + DIR_SEP "TSS" DIR_SEP "CurrentTSSData.txt";
}

std::string
BioEnvironment::get_tss_file() const
{
	return data_dir + DIR_SEP "TSS" DIR_SEP "TSS_mouse.txt";
}

std::string
BioEnvironment::get_tss_clones_file() const
{
	return data_dir + DIR_SEP "TSS" DIR_SEP "TSS_mouse_clones.txt";
}

std::string
BioEnvironment::get_factor_synonyms_file() const
{
	return get_serialised_dir() + DIR_SEP "factor_synonyms.txt";
}


BIO_NS_END
