#include <bio/fasta.h>
USING_BIO_NS

#include <boost/test/unit_test.hpp>
#include <boost/progress.hpp>
using namespace boost;
using boost::unit_test::test_suite;

#include <fstream>
using namespace std;


#define FASTA_FILE "C:\\data\\ensembl\\Homo_sapiens.NCBI35.may.dna.chromosome.21.fa"


void check_parse_fasta()
{
	cout << "******* check_parse_fasta(): Parsing fasta file " FASTA_FILE << endl;

	ifstream file(FASTA_FILE);
	BOOST_CHECK(file.good());

	stringstream sequence;
	progress_timer timer;
	parse_fasta(file, sequence);
}


void check_parse_fasta_2()
{
	const std::string file = "c:\\data\\Helene\\MyoD_DRR.fa";
	cout << "******* check_parse_fasta_2(): Parsing fasta file " << file << "\n";

	ifstream stream(file.c_str());
	BOOST_CHECK(stream.good());

	fasta_file_map_t sequences;
	parse_fasta_2(stream, sequences);
}


void
register_fasta_tests(boost::unit_test::test_suite * test)
{
	//test->add( BOOST_TEST_CASE( &check_parse_fasta ), 0);
	test->add( BOOST_TEST_CASE( &check_parse_fasta_2 ), 0);
}


