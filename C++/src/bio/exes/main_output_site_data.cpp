/* Copyright John Reid 2007
*/

#include "bio-pch.h"



#include <bio/application.h>
#include <bio/biobase_db.h>
USING_BIO_NS

#include <boost/iterator/filter_iterator.hpp>
#include <boost/test/utils/wrap_stringstream.hpp>
namespace po = boost::program_options;

#include <iostream>
using namespace std;


struct IsInterestingSite
{
  bool operator()(Site::map_t::value_type site) const
  {
	  return site.second->id.species_group != "AS";
  }
};


struct OutputSiteData
{
	void operator()(Site::map_t::value_type site) const
	{
		cout
			<< site.first
			<< "," << site.second->id
			<< "," << site.second->sequence
			<< "," << site.second->reference_point
			<< "," << site.second->start_position
			<< "," << site.second->end_position
			<< ",";
		std::copy(
			site.second->database_refs.begin(),
			site.second->database_refs.end(),
			ostream_iterator<DatabaseRef>(cout, ", "));
		cout
			<< "\n";
	}
};


struct OutputSiteDataApp : Application
{
	int
	task()
	{
		std::for_each(
			boost::make_filter_iterator(
				IsInterestingSite(),
				BiobaseDb::singleton().get_sites().begin(),
				BiobaseDb::singleton().get_sites().end()),
			boost::make_filter_iterator(
				IsInterestingSite(),
				BiobaseDb::singleton().get_sites().end(),
				BiobaseDb::singleton().get_sites().end()),
			OutputSiteData());

		return 0;
	}
};

int
main(int argc, char * argv[])
{
	return OutputSiteDataApp().main(argc, argv);
}

