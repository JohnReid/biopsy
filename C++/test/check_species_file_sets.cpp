
#include <bio/species_file_sets.h>
USING_BIO_NS

#include <boost/test/unit_test.hpp>

#include <iostream>
using namespace std;

void
check_species_file_sets()
{
	cout << "******* check_species_file_sets()" << endl;

	const unsigned num_seqs = 3;
	const unsigned seq_length = 10;

	SeqList seq_list;
	get_random_sequences(num_seqs, seq_length, seq_list);

	BOOST_CHECK_EQUAL(seq_list.size(), num_seqs);
	for (SeqList::const_iterator i = seq_list.begin();
		seq_list.end() != i;
		++i)
	{
		BOOST_CHECK_EQUAL(seq_length, i->size());
		BOOST_CHECK_EQUAL(unsigned(seq_length), unsigned(std::count_if(i->begin(), i->end(), is_known_nucleotide())));
	}
}


void
register_species_file_sets_tests(boost::unit_test::test_suite * test)
{
	test->add(BOOST_TEST_CASE(&check_species_file_sets), 0);
}

