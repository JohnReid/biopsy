/**
@file

Copyright John Reid 2006, 2007, 2013
*/

#include "bio_test_defs.h"
#include <bio/pssm_match.h>
#include <bio/matrix_test_data.h>
#include <bio/biobase_db.h>
#include <bio/biobase_filter.h>
#include <bio/biobase_data_traits.h>
#include <bio/run_match.h>
#include <bio/remo.h>
#include <bio/svg_match.h>
#include <bio/biobase_db.h>
#include <bio/biobase_data_traits.h>
USING_BIO_NS

#include <boost/foreach.hpp>
#include <boost/io/ios_state.hpp>
#include <boost/filesystem/path.hpp>
#include <boost/test/unit_test.hpp>
#include <boost/test/parameterized_test.hpp>
#include <boost/assign/list_of.hpp>
using namespace boost;
namespace fs = boost::filesystem;

#include <iostream>
#include <vector>
using namespace std;

//#define VERBOSE_CHECKING


/** Prints the hits in a format R can understand. */
template <class OStr>
OStr &
print_hits(OStr & ostr, const string & name, const hit_vec_t & hits)
{
    ostr << name << "<-c(";
    size_t pos = 0;
    for (BIO_NS::float_t score;
        (score = find_best_score_for_position(hits, pos)) >= 0.0;
        ++pos)
    {
        ostr << score << ",";
    }
    ostr << ")" << endl;

    return ostr;
}

void
check_pssm_match_unknowns()
{
    cout << "******* check_pssm_match_unknowns()" << endl;

    const seq_t seq = "actgactgacnacagctagtagctgacnactgatcgatcgatgcnatcgatcgatgcatgcatgcn";
    Scorers scorers(seq.begin(), seq.end());

#ifdef VERBOSE_CHECKING
    //cout << "Real sequence: " << test_data.binding_site->sequence << endl;
    cout << "Test sequence: " << seq << endl;

    print_hits(cout, "biobase", scorers.biobase_scores.get_result(*BiobaseDb::singleton().get_entry<MATRIX_DATA>(TableLink(MATRIX_DATA, 930))));
    print_hits(cout, "ott", scorers.ott_normalised_scores.get_result(*BiobaseDb::singleton().get_entry<MATRIX_DATA>(TableLink(MATRIX_DATA, 930))));
    print_hits(cout, "bayesian", scorers.bayesian_scores.get_result(*BiobaseDb::singleton().get_entry<MATRIX_DATA>(TableLink(MATRIX_DATA, 930))));
    print_hits(cout, "bob", scorers.bayesian_or_better_scores.get_result(*BiobaseDb::singleton().get_entry<MATRIX_DATA>(TableLink(MATRIX_DATA, 930))));
#endif
}

void
check_deaf_match()
{
    cout << "******* check_deaf_match()" << endl;

    vector<TableLink> matrices = boost::assign::list_of
        ( TableLink(MATRIX_DATA, 1001) )
        ;

    vector< seq_t > sequences = boost::assign::list_of
        ( "NNNNNNNNNNNNNNNNNNNNNNNNN" )
        ;

    BOOST_FOREACH( TableLink link, matrices )
    {
        Pssm pssm = make_pssm( link );

        BOOST_FOREACH( const seq_t & s, sequences )
        {
            pssm.score( s.begin(), false );
        }
    }
}


void
check_pssm_match()
{
    cout << "******* check_pssm_match()" << endl;

    boost::io::ios_precision_saver ips(cout);
    cout.precision(3);

    vector<TableLink> matrices;
    matrices.push_back(TableLink(MATRIX_DATA, 930));

    for (vector<TableLink>::const_iterator m = matrices.begin();
        matrices.end() != m;
        ++m)
    {
        MatrixTestData test_data(
            BiobaseDb::singleton().get_entry< MATRIX_DATA >( *m ),
            40 );

        Scorers scorers( test_data.get_input().centre_sequence.begin(), test_data.get_input().centre_sequence.end() );

#ifdef VERBOSE_CHECKING
        //cout << "Real sequence: " << test_data.binding_site->sequence << endl;
        cout << "Real sequence at: " << test_data.binding_site_idx << endl;
        cout << "Test sequence: " << test_data.seq << endl;

        print_hits(cout, "biobase", scorers.biobase_scores.get_result(*test_data.matrix));
        print_hits(cout, "ott", scorers.ott_normalised_scores.get_result(*test_data.matrix));
        print_hits(cout, "bayesian", scorers.bayesian_scores.get_result(*test_data.matrix));
        print_hits(cout, "bob", scorers.bayesian_or_better_scores.get_result(*test_data.matrix));
#endif
    }
}

struct PhyloAdjustParams
{
    typedef std::map<Species, double> ConsMap;

    std::string remo_name;
    Species main_species;
    ConsMap conservations;
    ScoreAlgorithm algorithm;
    BIO_NS::float_t threshold;
};

void
check_phylo_adjust(const PhyloAdjustParams & params)
{
    cout << "******* check_phylo_adjust()" << endl;

    build_test_remos();

    Remo & remo = *(test_remos[params.remo_name]);

    //score the sites and matrices
    match_result_vec_t results;
    BiobasePssmFilter filter;
    pssm_match(
        remo.map[params.main_species].begin(),
        remo.map[params.main_species].end(),
        power(params.threshold, 10),
        params.algorithm,
        false,
        filter,
        filter,
        std::inserter(results, results.begin()));

#ifdef VERBOSE_CHECKING
    cout << results.size() << " matches over threshold" << endl;
#endif



    const bool show_svg =
#ifdef VERBOSE_CHECKING
        true;
#else
        false;
#endif

#if 1
    {
        const std::string name = BIO_MAKE_STRING(params.remo_name << "_" << params.main_species << "_before_adjusting.svg");
        build_svg(
            fs::path(name),
            name,
            remo.map[MOUSE_SPECIES],
            power(params.threshold, 10),
            results,
            16,
            false,
            show_svg);
    }
#endif

    double new_threshold = params.threshold;
    for (PhyloAdjustParams::ConsMap::const_iterator i = params.conservations.begin();
        params.conservations.end() != i;
        ++i)
    {
        adjust_hits(
            results.begin(),
            results.end(),
            remo.map[i->first].begin(),
            remo.map[i->first].end(),
            power(params.threshold, 10),
            params.algorithm);

        match_result_vec_t phylo_results;
        pssm_match(
            remo.map[i->first].begin(),
            remo.map[i->first].end(),
            power(params.threshold, 10),
            params.algorithm,
            false,
            filter,
            filter,
            std::inserter(phylo_results, phylo_results.begin()));

#ifdef VERBOSE_CHECKING
        cout << results.size() << " matches over threshold" << endl;
#endif


#if 1
        const std::string name =
            BIO_MAKE_STRING(params.remo_name << "_" << i->first << "_phylogenetic.svg");
        build_svg(
            fs::path(name),
            name,
            remo.map[i->first],
            power(params.threshold, 10),
            phylo_results,
            16,
            false,
            show_svg);
#endif

        new_threshold *= params.threshold;
    }

    const std::string name =
        BIO_MAKE_STRING(params.remo_name << "_" << params.main_species << "_after_adjusting.svg");
    build_svg(
        fs::path(name),
        name,
        remo.map[MOUSE_SPECIES],
        new_threshold,
        results,
        16,
        false,
        show_svg);
}

void
check_phylo_adjust_multiple()
{
    cout << "******* check_phylo_adjust_multiple()" << endl;

    {
        PhyloAdjustParams params;
        params.remo_name = "dlx5_1";
        params.algorithm = BAYESIAN_SCORE_ALGORITHM;
        params.conservations[HUMAN_SPECIES] = .9;
        params.conservations[XENOPUS_SPECIES] = .65;
        params.main_species = MOUSE_SPECIES;
        params.threshold = 0.999f;

        check_phylo_adjust(params);
    }

    {
        PhyloAdjustParams params;
        params.remo_name = "egr1_1";
        params.algorithm = BAYESIAN_SCORE_ALGORITHM;
        params.main_species = MOUSE_SPECIES;
        params.threshold = 0.9998f;
        params.conservations[HUMAN_SPECIES] = .9;
        params.conservations[DOG_SPECIES] = .65;
        //params.conservations[COW_SPECIES] = .8;
        //params.conservations[CHICKEN_SPECIES] = .8;

        check_phylo_adjust(params);
    }

    {
        PhyloAdjustParams params;
        params.remo_name = "egr1_2";
        params.algorithm = BAYESIAN_SCORE_ALGORITHM;
        params.main_species = MOUSE_SPECIES;
        params.threshold = 0.98f;
        params.conservations[HUMAN_SPECIES] = .9;
        params.conservations[DOG_SPECIES] = .65;
        params.conservations[COW_SPECIES] = .8;
        params.conservations[FUGU_SPECIES] = .8;

        check_phylo_adjust(params);
    }
}

void
check_bayesian_result_calculator(const seq_t & seq)
{
    cout << "******* check_bayesian_result_calculator(): \"" << seq << "\"" << endl;

    //get the pssm
    const Matrix * matrix = BiobaseDb::singleton().get_entry<MATRIX_DATA>(TableLink(MATRIX_DATA, 260));

    Scorers scorers(seq.begin(), seq.end());
    const hit_vec_t results = scorers.bayesian_scores.get_result(*matrix);

#ifdef VERBOSE_CHECKING
    std::copy(results.begin(), results.end(), ostream_iterator<Hit>(cout, "; "));
    cout << endl;
#endif
}

const SeqList bayesian_score_seqs = boost::assign::list_of
    ("GTTACGCAAT") //should match 100%
    ("ATTGCGTAAC") //complementary palindrome of above
    ("GCTACGCAAT") //first mutated in 1 position
    ("GCGACGCAAT") //first mutated in 2 positions
    ("GCGGCGCAAT") //first mutated in 3 positions
    ("CCCTACGTGG") //worst match possible
    ;


void
check_bayesian_score_distributions()
{
    cout << "******* check_bayesian_score_distributions()" << endl;

    static const seq_t mouse_msx1_remo =
        "CCAAGCCACCACAGTGCTAATTCTGCAGGCTGTC"
        "ATCGTTGGGCATCACTGCGTGCTGATGGTCGGCCGAGTGAGTGTAAATAGAGTCAAGATT"
        "TTGCTCTACAAATGCCTTGGCGCACAGCTGCTCCTGGGAACACCCCCCCCTCCCCCCCCC"
        "GTCCATGCCATGTTTGAAGCGACTTCTAATTGGCCTCATCCTGACCCAGCCTTGGAGAGA"
        "GGCCCTAGCTGGTCTCAGCCACCATGAAAGGCAGCAGATGGGGTTCGGAGCACCCCTGAG"
        "CATCACTACTGGGTGTGGAAACTGGCAGAAGCTGATGTTTGGAAGACAGCAATCCAGAGG"
        "CCCAGCTCTGCCACTTTAAAGTTCTGCCTAGTGCTGGGGACAAAACCAAACCAAACCAAA"
        "ACTTTCCCTCTGTCAAATTTCAAGCCCAAAATTTACTGCACATACTTTTGGAGCCTTGTG"
        "TTTGATCCGCTTTTAATTAACTCATTTTGAAGCCTTCATACTTTCCTTCAAAATATTTTA"
        "TTGGCATTTCATTAGCATTGTCTTCCCTGGGCAGGCGGCGATACAGCCAGGGAGGGTCAC"
        "CCAGCAACATGGCTTGGCATTTGGACATGCTTTTGTTCCTGTCTCCCTGCCAGGCCTCTC"
        "TGGGGCCCATGCATTTACTGTTAAAATATTACAGAGGCCCACACCCGAACGACCATAATC"
        "CCTCCTGCTAGACTGTCTCTAATTGATTTAATTATTCATAATTGCTCAGCTGGGTTTGCG"
        "CATACTTGTGTATGATTATGGAGAGAGGGTGGGACCTTGCTGAGTGCTGGGTCCCCATCA"
        "TGATTAGCAGGAGAAT";

    const Matrix * matrix = BiobaseDb::singleton().get_entry<MATRIX_DATA>(TableLink(MATRIX_DATA, 260));

    Scorers scorers(mouse_msx1_remo.begin(), mouse_msx1_remo.end());
    const hit_vec_t results = scorers.bayesian_scores.get_result(*matrix);

    //count scores in a range
    const unsigned num_buckets = 100;
    std::map<unsigned, unsigned> buckets;
    for (hit_vec_t::const_iterator h = results.begin();
        results.end() != h;
        ++h)
    {
        const unsigned bucket = unsigned(h->score * num_buckets);
        if (buckets.end() == buckets.find(bucket))
        {
            buckets[bucket] = 0;
        }
        buckets[bucket]++;
    }

#ifdef VERBOSE_CHECKING
    //output counts
    for (unsigned b = 0; num_buckets != b; ++b)
    {
        cout << b << ": " << buckets[b] << "\n";
    }
#endif
}

void
register_pssm_match_tests(boost::unit_test::test_suite * test)
{
    test->add( BOOST_TEST_CASE( &check_deaf_match ), 0);
    test->add(
        BOOST_PARAM_TEST_CASE(
            &check_bayesian_result_calculator,
            bayesian_score_seqs.begin(),
            bayesian_score_seqs.end()),
        0);
    test->add( BOOST_TEST_CASE( &check_pssm_match ), 0);
    test->add( BOOST_TEST_CASE( &check_bayesian_score_distributions ), 0);
    test->add( BOOST_TEST_CASE( &check_phylo_adjust_multiple ), 0);
    test->add( BOOST_TEST_CASE( &check_pssm_match_unknowns ), 0);
}


//
// Are we going to compile this test into its own executable?
//
#ifdef BIO_STANDALONE_TEST
BIO_DEFINE_STANDALONE_TEST( "pssm match", register_pssm_match_tests )
#endif //BIO_STANDALONE_TEST

