

#ifndef BIO_OPTIONS_HPP
#define BIO_OPTIONS_HPP

#include "bio/defs.h"
#include "bio/singleton.h"
#include "bio/environment.h"

#include <boost/program_options.hpp>

#include <string>


BIO_NS_START





struct BioOptions
	: Singleton< BioOptions >
{
	BioOptions( BioEnvironment & environment = BioEnvironment::singleton() );

	/** These options allowed on the command line. */
	boost::program_options::options_description cmd_line_options;

	/** These options allowed on the command line and a config file. */
	boost::program_options::options_description config_options;

	/** Hidden options, will be allowed both on command line and
    in config file, but will not be shown to the user. */
	boost::program_options::options_description hidden_options;

	/** The variables map used in the parsing. */
    boost::program_options::variables_map values;

	/** Parse the command line options and the config stream. */
	void parse(
		int argc,
		char * argv [],
		std::basic_istream<char> * config_file = 0, //0 for no file
		const boost::program_options::options_description & additional_options
			= boost::program_options::options_description(),
		const boost::program_options::positional_options_description & additional_position_options
			= boost::program_options::positional_options_description());

	/** Return the options that are visible to the user. */
	boost::program_options::options_description get_visible_options() const;

	/** True iff --help was on command line. */
	bool help;

	/** True iff --version was on command line. */
	bool version;

	static
	void parse_config_file(
		const std::string & filename,
		boost::program_options::variables_map & values,
		boost::program_options::options_description & options);

	/** Get all the options that can be used on the command line. */
	boost::program_options::options_description get_cmd_line_options();

	/** Get all the options that can be used in the config file. */
	boost::program_options::options_description get_config_file_options();

};





BIO_NS_END

#endif //BIO_OPTIONS_HPP
