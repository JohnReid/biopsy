
#ifndef BIO_SITE_DATA_H_
#define BIO_SITE_DATA_H_


#include "bio/defs.h"
#include "bio/sequence.h"
#include "bio/common.h"

#include <vector>
#include <iostream>



BIO_NS_START

struct SiteDataGene
{
	typedef std::vector<SiteDataGene> vector_t;

	std::string ensembl_id;
	int index;
	int length;
	std::string region_name;
	int start;
	int end;
	seq_t sequence;
};

struct SiteData
{
	typedef boost::shared_ptr<SiteData> ptr_t;
	typedef std::vector<ptr_t> vector_t;

	TableLink site;
	seq_t site_sequence;
	SiteDataGene::vector_t genes;

	static bool parse(std::istream & stream, vector_t & result);
};


BIO_NS_END


#endif //BIO_SITE_DATA_H_
