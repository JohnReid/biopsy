/* Copyright John Reid 2007
*/

#include "bio-pch.h"


#include "bio/defs.h"


#include "bio/matrix_match.h"
#include "bio/environment.h"
#include "MatrixMatchParser.hpp"

#include <boost/shared_ptr.hpp>
#include <boost/filesystem/fstream.hpp>
namespace fs = boost::filesystem;

#include <string>
using namespace std;

BIO_NS_START


MatrixMatch::MatrixMatch()
: threshold(-1.0)
{
}

void
parse_match_thresholds(
	const fs::path & match_thresholds_file,
	MatrixMatch::map_t & match_map )
{
	if( ! fs::exists( match_thresholds_file ) )
	{
		std::cout << "MatrixMatch: \"" << match_thresholds_file._BOOST_FS_NATIVE() << "\" does not exist" << std::endl;
	}
	else
	{
		fs::ifstream str( match_thresholds_file );
		if (! str)
		{
			std::cout << "MatrixMatch: \"" << match_thresholds_file._BOOST_FS_NATIVE() << "\" cannot open" << std::endl;
		}
		else
		{

			char line[1024];
			str.getline(line, sizeof(line) - 1);
			str.getline(line, sizeof(line) - 1);
			str.getline(line, sizeof(line) - 1);
			str.getline(line, sizeof(line) - 1);

			while (str)
			{
				//did we reach the end
				if ('/' == str.peek())
				{
					break;
				}

				float_t one, three_quarters, min_fp;
				string acc_num, name;
				str >> one >> three_quarters >> min_fp >> acc_num >> name;

				try
				{
					const TableLink link = parse_table_link_accession_number(acc_num);
					match_map[link].threshold = min_fp;
					str.getline(line, sizeof(line) - 1);
				}
				catch (...)
				{
				}
			}

			str.getline(line, sizeof(line) - 1);
			str.get();
			if (! str.eof())
			{
				std::cout << "MatrixMatch: \"" << match_thresholds_file._BOOST_FS_NATIVE() << "\" could not parse everything" << std::endl;
			}
		}
	}
}

const MatrixMatch::map_t & get_min_fp_match_map()
{
	static boost::scoped_ptr< MatrixMatch::map_t > map;
	if( 0 == map )
	{
		map.reset( new MatrixMatch::map_t );
		parse_match_thresholds(
			fs::path(BioEnvironment::singleton().get_matrix_min_fp_file()),
			*map );
	}

	return *map;
}

const MatrixMatch::map_t & get_min_fn_match_map()
{
	static boost::scoped_ptr< MatrixMatch::map_t > map;
	if( 0 == map )
	{
		map.reset( new MatrixMatch::map_t );
		parse_match_thresholds(
			fs::path(BioEnvironment::singleton().get_matrix_min_fn_file()),
			*map );
	}

	return *map;
}

const MatrixMatch::map_t & get_min_sum_match_map()
{
	static boost::scoped_ptr< MatrixMatch::map_t > map;
	if( 0 == map )
	{
		map.reset( new MatrixMatch::map_t );
		parse_match_thresholds(
			fs::path(BioEnvironment::singleton().get_matrix_min_sum_file()),
			*map );
	}

	return *map;
}

BIO_NS_END
