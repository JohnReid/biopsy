/* Copyright John Reid 2007
*/

#include "bio-pch.h"




#include <bio/application.h>
#include <bio/biobase_db.h>
USING_BIO_NS;


#include <boost/program_options.hpp>
using namespace boost;
namespace po = boost::program_options;
namespace fs = boost::filesystem;

#include <iostream>
#include <fstream>
using namespace std;



struct ParseBiobaseApp : Application
{
	ParseBiobaseApp()
	{
		get_options().add_options()
			;
	}

	int task()
	{
		{
			cout << "Loading all of biobase databases\n";
			boost::progress_timer timer;
			BiobaseDb::singleton().load_all();
			cout << "Took ";
		}

		cout << "# matrices:  " << BiobaseDb::singleton().get_matrices().size() << "\n";
		cout << "# sites:     " << BiobaseDb::singleton().get_sites().size() << "\n";
		cout << "# factors:   " << BiobaseDb::singleton().get_factors().size() << "\n";
		cout << "# fragments: " << BiobaseDb::singleton().get_fragments().size() << "\n";
		cout << "# genes:     " << BiobaseDb::singleton().get_genes().size() << "\n";
		cout << "# compels:   " << BiobaseDb::singleton().get_compels().size() << "\n";
		cout << "# evidences: " << BiobaseDb::singleton().get_evidences().size() << "\n";
//		cout << "# pathways:  " << BiobaseDb::singleton().get_pathways().size() << "\n";
//		cout << "# molecules: " << BiobaseDb::singleton().get_molecules().size() << "\n";
		//cout << "finished: type any key to continue\n";
		//char c;
		//std::cin >> c;

		return 0;
	}
};

int
main(int argc, char * argv [])
{
	return ParseBiobaseApp().main(argc, argv);
}
