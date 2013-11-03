

#ifndef BIO_TSS_DATA_H_
#define BIO_TSS_DATA_H

#include "bio/defs.h"

#include <boost/filesystem/path.hpp>

#include <map>
#include <string>




BIO_NS_START

enum CloneType
{
	RIKEN_CLONE,
	TIGR_CLONE,
	NO_CLONE,
	NO_CONCLUSION_CLONE
};


struct Clone
{
	std::string transcript;
	CloneType type;
	std::string id;
	std::string chromosome;
	bool strand;
	int rel_tss_pos;
	int tss_pos;

	typedef std::vector<Clone> vec_t;
};



struct TSS
{
	Clone::vec_t clones;

	typedef std::map<std::string, TSS> map_t; /**< Maps gene ids to TSSs. */

	static void parse_files(
		boost::filesystem::path tss_file,
		boost::filesystem::path clone_file,
		map_t & map);
};



BIO_NS_END

#endif //BIO_TSS_DATA_H


