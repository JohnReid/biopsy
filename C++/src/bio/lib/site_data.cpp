/* Copyright John Reid 2007
*/

#include "bio-pch.h"


#include "bio/defs.h"


#include "bio/site_data.h"

#include<boost/tokenizer.hpp>

#include <string>
using namespace std;



BIO_NS_START

bool
SiteData::parse(std::istream & stream, SiteData::vector_t & result)
{
	typedef boost::tokenizer< boost::char_separator<char> > tokenizer;
	unsigned line_count = 0;

	//for each entry
	while (stream)
	{
		SiteData::ptr_t site_data(new SiteData);

		{
			//get the line
			char buf[16384];
			stream.getline(buf, sizeof(buf));
			++line_count;
			string line(buf);

			//tokenize the line
			boost::char_separator<char> sep(",");
			tokenizer tok(line, sep);

			//extract the values
			tokenizer::iterator tok_iter = tok.begin();

			site_data->site.table_id = SITE_DATA;
			if (tok.end() == tok_iter)
			{
				return false;
			}
			site_data->site.entry_idx = atoi(tok_iter->c_str());
			++tok_iter;

			if (tok.end() == tok_iter)
			{
				return false;
			}
			++tok_iter;

			if (tok.end() == tok_iter)
			{
				return false;
			}
			site_data->site_sequence = *tok_iter;
			++tok_iter;
		}

		while (stream)
		{
			SiteDataGene gene;

			//get the next line
			char buf[32768];
			stream.getline(buf, sizeof(buf));
			++line_count;
			string line(buf);

			//tokenize the line
			boost::char_separator<char> sep(",");
			tokenizer tok(line, sep);

			//extract the values
			tokenizer::iterator tok_iter = tok.begin();

			//is it a blank line?
			if (tok_iter == tok.end())
			{
				break;
			}

			gene.ensembl_id = *tok_iter;
			++tok_iter;

			if (tok.end() == tok_iter)
			{
				return false;
			}
			gene.index = atoi(tok_iter->c_str());
			++tok_iter;

			if (tok.end() == tok_iter)
			{
				return false;
			}
			gene.length = atoi(tok_iter->c_str());
			++tok_iter;

			if (tok.end() == tok_iter)
			{
				return false;
			}
			gene.region_name = *tok_iter;
			++tok_iter;

			if (tok.end() == tok_iter)
			{
				return false;
			}
			gene.start = atoi(tok_iter->c_str());
			++tok_iter;

			if (tok.end() == tok_iter)
			{
				return false;
			}
			gene.end = atoi(tok_iter->c_str());
			++tok_iter;

			if (tok.end() != tok_iter)
			{
				gene.sequence = *tok_iter;
				++tok_iter;
			}

			//we should be at the end
			if (tok.end() != tok_iter)
			{
				return false;
			}

			site_data->genes.push_back(gene);
		}

		result.push_back(site_data);
	}

	const bool reached_end = stream.eof();

	return reached_end;
}


BIO_NS_END
