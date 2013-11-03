/* Copyright John Reid 2007
*/

#include "bio-pch.h"


#include "bio/defs.h"


#include "bio/options.h"
#include "bio/environment.h"

#include <boost/filesystem/path.hpp>
#include <boost/filesystem/operations.hpp>
#include <boost/filesystem/fstream.hpp>
namespace po = boost::program_options;
namespace fs = boost::filesystem;


#include <iostream>
using namespace std;


BIO_NS_START


BioOptions::BioOptions(BioEnvironment & env)
: help(false)
, version(false)
{
	config_options.add_options()
        ("help,h", po::bool_switch(&help), "print usage message")
        ("version", po::bool_switch(&version), "print version")
		("data_dir", po::value(&env.data_dir), "where the data is stored")
		("transpath_major_version", po::value(&env.transpath_major_version), "TRANSPATH major version")
		("transpath_minor_version", po::value(&env.transpath_minor_version), "TRANSPATH minor version")
		("transcompel_major_version", po::value(&env.transcompel_major_version), "TRANSCOMPEL major version")
		("transcompel_minor_version", po::value(&env.transcompel_minor_version), "TRANSCOMPEL minor version")
		("transfac_major_version", po::value(&env.transfac_major_version), "TRANSFAC major version")
		("transfac_minor_version", po::value(&env.transfac_minor_version), "TRANSFAC minor version")
		("custom_PSSM_version", po::value(&env.custom_PSSM_version), "custom PSSM version")
		("biobase_url_prefix", po::value(&env.biobase_url_prefix), "prefix for biobase URLs")
		("tf_binding_prior", po::value(&env.tf_binding_prior), "prior for TF binding prob")
		("http_port", po::value(&env.http_port), "the port for the http server")
		("log_file", po::value(&env.log_file_name), "the name of the file to log to")
		;

	hidden_options.add_options()
		("species_prefixes", po::value(&env.species_prefixes), "")
		("num_normalisation_quanta", po::value(&env.num_normalisation_quanta), "")
		("max_pssm_likelihood_map_size", po::value(&env.max_pssm_likelihood_map_size), "")
		;
}

boost::program_options::options_description
BioOptions::get_cmd_line_options()
{
    return po::options_description().add(cmd_line_options).add(config_options).add(hidden_options);
}

boost::program_options::options_description
BioOptions::get_config_file_options()
{
	return po::options_description().add(config_options).add(hidden_options);
}

void
BioOptions::parse_config_file(
	const std::string & filename,
	boost::program_options::variables_map & values,
	boost::program_options::options_description & options)
{
	//Try and find the default config file looking in the current directory and then every parent directory
	//until we find a suitably named file
	fs::path dir(".");
	while (fs::exists(dir))
	{
		fs::path config_file = dir / filename;
		if (fs::exists(config_file))
		{
			cout << "Parsing default config file: " << fs::system_complete( config_file ).normalize()._BOOST_FS_NATIVE() << endl;

			//parse the default config file
			fs::ifstream stream(config_file);
			store(po::parse_config_file(stream, options), values);
			break;
		}
		dir /= "..";
	}
	if (! fs::exists(dir))
	{
		cout << "Could not find default config file: " << filename << endl;
	}

    notify(values);
}

void
BioOptions::parse(
	int argc,
	char * argv [],
	std::basic_istream<char> * config_stream, //0 for default config file
	const po::options_description & additional_options,
	const po::positional_options_description & additional_position_options)
{
    po::options_description combined_cmdline_options = get_cmd_line_options().add(additional_options);

    po::options_description combined_config_options = get_config_file_options().add(additional_options);

	//parse the command line
	po::store(
		po::command_line_parser(argc, argv)
            .options(combined_cmdline_options)
			.positional(additional_position_options)
			.run(),
		values);

	//parse the given config stream
	if (0 != config_stream)
	{
		po::store(po::parse_config_file(*config_stream, combined_config_options), values);
	}

	//Try and find the default config file looking in the current directory and then every parent directory
	//until we find a suitably named file
	parse_config_file("bio_lib.cfg", values, combined_config_options);
    notify(values);
}

po::options_description
BioOptions::get_visible_options() const
{
	po::options_description visible("Bio library options");
	visible.add(cmd_line_options).add(config_options);

	return visible;
}

BIO_NS_END


