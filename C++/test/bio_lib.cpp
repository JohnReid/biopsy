/**
@file

Copyright John Reid 2013
*/

#include "bio_test_defs.h"
#include "bio_test_data.h"

#include <bio/biobase_match.h>
#include <bio/matrix.h>
#include <bio/score_sequence.h>
USING_BIO_NS;

#include <boost/io/ios_state.hpp>
#include <boost/progress.hpp>
using namespace boost;

#include <string>
using namespace std;


//vector of parameters
vector<BiobaseParams> biobase_params;

//data we use as parameter
const seq_t test_rnd_1000_seq = "tacatcatctgtctgcagtagtctaaccgaccccccccagttttagaagcagactgcatgcggacgggaccgcggatcgcgcggtgcgcctcagtgtacttccgaacgaatgagtcattaatagagcgctatatcgtaactgtctttgacgaagtataccgaaaccgtgcagccagacgtgatccgggcgttgtaaaggcgatcagcgccctaggagtaccatttttgccgtaggcttgcgtctcaaagaccagctggggcgtggtatcactcgtcagtacgatttctgccagatagatagcatagactgaaccttaggcccaatagggacacaattacccgagtgactgactggtctaaggggagtccccccttaaaacgttttacgtaatagcgggctccagaagcaaagcatcggtttgagccccagtactaaacgtttgagtgtttgctctcgtctgataggtaaaccgacaagagaaccaagctcaaggcgcggtaggtgcgccttgcgaactgttgatgccgtgagcgccaccatcccgtgcatcataggcagggagagaagaccacatggccttgcgaccgtatgagctgtttcagattaaatgccaacgggcatggtcggtgtccagcattttttgcagtcagctggtggtacacagtggggacaagaacgcctctggtagatgtcttctgaaggagtaactcatttcgttgaatcgaccttcccttgcgcttgaacgcggacctctagtctctctcgcagactggggtcgaaaatcaaggtagatatggaatgttccgcatgagggtagcgaccggatcgggcgtcaagtatatcctccctgctacgtccccctactagcctcagtccgcctcgaacctaggaagattggccacatcagcttggtggatgcctggtccatacttcagacccgagaatgttagacaggaccccatttggctcctttacgtacgatctatgtagacgcagtga";
const seq_t acgtacgt_seq = "acgtacgt";
const seq_t atcgatcg_seq = "atcgatcg"; //complement is tagctagc - reversed complement is cgatcgat
Pssm::ptr_t acgt_pssm;
Pssm::ptr_t atcg_pssm;
Pssm::ptr_t V$DEAF1_01_pssm;
Pssm::ptr_t V$IPF1_Q4_01_pssm;
Pssm::ptr_t V$EFC_Q6_pssm;
Pssm::ptr_t V$AP1_Q4_pssm;
IupacSeq acgt_iupac;
IupacSeq atcg_iupac;
IupacSeq V$DEAF1_01_iupac;
IupacSeq V$IPF1_Q4_01_iupac;
IupacSeq V$EFC_Q6_iupac;
IupacSeq V$AP1_Q4_iupac;


void ensure_biobase_params_built()
{
    static bool already_done = false;
    if (! already_done)
    {
        //acgt is completely complementary palindromic so we score 1.0 when matching complement
        acgt_pssm.reset(new Pssm);
        acgt_pssm->push_back(PssmEntry(1,0,0,0));
        acgt_pssm->push_back(PssmEntry(0,1,0,0));
        acgt_pssm->push_back(PssmEntry(0,0,1,0));
        acgt_pssm->push_back(PssmEntry(0,0,0,1));
        acgt_iupac.clear();
        for (const char * cp = "acgt"; *cp != '\0'; ++cp) {
            acgt_iupac.push_back(*cp);
        }
        assert(acgt_iupac.size() == acgt_pssm->size());
        biobase_params.push_back(
            BiobaseParams(
                "acgt non-complementary",
                acgt_pssm,
                &acgt_iupac,
                acgtacgt_seq.begin(),
                acgtacgt_seq.begin() + 4,
                false,
                false,
                0.0f,
                Hit(1.0, 0)));
        biobase_params.push_back(
            BiobaseParams(
                "acgt complementary",
                acgt_pssm,
                &acgt_iupac,
                acgtacgt_seq.begin(),
                acgtacgt_seq.begin() + 4,
                true,
                false,
                0.0f,
                Hit(1.0, 0, false)));

        //atcg is not complementary palindromic at all so we score 0.0 when matching complement
        atcg_pssm.reset(new Pssm);
        atcg_pssm->push_back(PssmEntry(1,0,0,0));
        atcg_pssm->push_back(PssmEntry(0,0,0,1));
        atcg_pssm->push_back(PssmEntry(0,1,0,0));
        atcg_pssm->push_back(PssmEntry(0,0,1,0));
        atcg_iupac.clear();
        for (const char * cp = "atcg"; *cp != '\0'; ++cp) {
            atcg_iupac.push_back(*cp);
        }
        assert(atcg_iupac.size() == atcg_pssm->size());
        biobase_params.push_back(
            BiobaseParams(
                "atcg non-complementary",
                atcg_pssm,
                &atcg_iupac,
                atcgatcg_seq.begin(),
                atcgatcg_seq.begin() + 4,
                false,
                false,
                0.0f,
                Hit(1.0, 0)));
        biobase_params.push_back(
            BiobaseParams(
                "atcg complementary",
                atcg_pssm,
                &atcg_iupac,
                atcgatcg_seq.begin(),
                atcgatcg_seq.begin() + 4,
                true,
                false,
                0.0f,
                Hit(0.0, 0, false)));

        //the random sequence was tested on the TRANSFAC web site to get the following scores
        V$DEAF1_01_pssm.reset(new Pssm);
        V$DEAF1_01_pssm->push_back(PssmEntry(4,7,4,1));
        V$DEAF1_01_pssm->push_back(PssmEntry(0,13,3,2));
        V$DEAF1_01_pssm->push_back(PssmEntry(4,2,10,2));
        V$DEAF1_01_pssm->push_back(PssmEntry(1,9,4,4));
        V$DEAF1_01_pssm->push_back(PssmEntry(4,9,4,3));
        V$DEAF1_01_pssm->push_back(PssmEntry(2,10,1,9));
        V$DEAF1_01_pssm->push_back(PssmEntry(0,5,2,16));
        V$DEAF1_01_pssm->push_back(PssmEntry(2,16,1,3));
        V$DEAF1_01_pssm->push_back(PssmEntry(2,2,19,0));
        V$DEAF1_01_pssm->push_back(PssmEntry(2,2,15,4));
        V$DEAF1_01_pssm->push_back(PssmEntry(1,3,17,2));
        V$DEAF1_01_pssm->push_back(PssmEntry(5,1,7,10));
        V$DEAF1_01_pssm->push_back(PssmEntry(10,0,10,3));
        V$DEAF1_01_pssm->push_back(PssmEntry(1,5,2,15));
        V$DEAF1_01_pssm->push_back(PssmEntry(2,3,0,18));
        V$DEAF1_01_pssm->push_back(PssmEntry(0,4,1,18));
        V$DEAF1_01_pssm->push_back(PssmEntry(1,22,0,0));
        V$DEAF1_01_pssm->push_back(PssmEntry(0,22,1,0));
        V$DEAF1_01_pssm->push_back(PssmEntry(0,0,19,1));
        V$DEAF1_01_pssm->push_back(PssmEntry(5,0,8,6));
        V$DEAF1_01_pssm->push_back(PssmEntry(10,3,4,2));
        V$DEAF1_01_pssm->push_back(PssmEntry(7,1,8,3));
        V$DEAF1_01_pssm->push_back(PssmEntry(4,3,5,5));
        V$DEAF1_01_pssm->push_back(PssmEntry(5,4,3,4));
        V$DEAF1_01_pssm->push_back(PssmEntry(2,4,7,3));
        V$DEAF1_01_iupac.clear();
        for (const char * cp = "NCGNNYTCGGGNRTTTCCGDARNNN"; *cp != '\0'; ++cp) {
            V$DEAF1_01_iupac.push_back(*cp);
        }
        assert(V$DEAF1_01_iupac.size() == V$DEAF1_01_pssm->size());
        biobase_params.push_back(
            BiobaseParams(
                "V$DEAF1_01",
                V$DEAF1_01_pssm,
                &V$DEAF1_01_iupac,
                test_rnd_1000_seq.begin(),
                test_rnd_1000_seq.end(),
                true,
                true,
                21.0f,
                Hit(0.857f, 85)));

        V$IPF1_Q4_01_pssm.reset(new Pssm);
        V$IPF1_Q4_01_pssm->push_back(PssmEntry(3,3,1,8));//,T
        V$IPF1_Q4_01_pssm->push_back(PssmEntry(1,6,7,2));//,S
        V$IPF1_Q4_01_pssm->push_back(PssmEntry(3,2,7,4));//,N
        V$IPF1_Q4_01_pssm->push_back(PssmEntry(3,4,9,0));//,G
        V$IPF1_Q4_01_pssm->push_back(PssmEntry(2,5,0,9));//,Y
        V$IPF1_Q4_01_pssm->push_back(PssmEntry(3,12,0,1));//,C
        V$IPF1_Q4_01_pssm->push_back(PssmEntry(16,0,0,0));//,A
        V$IPF1_Q4_01_pssm->push_back(PssmEntry(0,0,0,16));//,T
        V$IPF1_Q4_01_pssm->push_back(PssmEntry(0,1,0,15));//,T
        V$IPF1_Q4_01_pssm->push_back(PssmEntry(15,0,0,1));//,A
        V$IPF1_Q4_01_pssm->push_back(PssmEntry(4,2,7,3));//,N
        V$IPF1_Q4_01_pssm->push_back(PssmEntry(6,4,1,5));//,N
        V$IPF1_Q4_01_pssm->push_back(PssmEntry(3,5,5,3));//,N
        V$IPF1_Q4_01_pssm->push_back(PssmEntry(3,4,2,7));//,N
        V$IPF1_Q4_01_pssm->push_back(PssmEntry(2,8,4,2));//,C
        V$IPF1_Q4_01_iupac.clear();
        for (const char * cp = "TSNGYCATTANNNNC"; *cp != '\0'; ++cp) {
            V$IPF1_Q4_01_iupac.push_back(*cp);
        }
        assert(V$IPF1_Q4_01_iupac.size() == V$IPF1_Q4_01_pssm->size());
        biobase_params.push_back(
            BiobaseParams(
                "V$IPF1_Q4_01",
                V$IPF1_Q4_01_pssm,
                &V$IPF1_Q4_01_iupac,
                test_rnd_1000_seq.begin(),
                test_rnd_1000_seq.end(),
                true,
                true,
                12.0f,
                Hit(0.975f, 110)));

        V$EFC_Q6_pssm.reset(new Pssm);
        V$EFC_Q6_pssm->push_back(PssmEntry(3,3,0,0));//,M
        V$EFC_Q6_pssm->push_back(PssmEntry(2,0,4,0));//,R
        V$EFC_Q6_pssm->push_back(PssmEntry(2,0,0,4));//,W
        V$EFC_Q6_pssm->push_back(PssmEntry(0,2,0,4));//,Y
        V$EFC_Q6_pssm->push_back(PssmEntry(4,0,2,0));//,R
        V$EFC_Q6_pssm->push_back(PssmEntry(1,5,0,0));//,C
        V$EFC_Q6_pssm->push_back(PssmEntry(0,2,1,3));//,Y
        V$EFC_Q6_pssm->push_back(PssmEntry(4,1,0,1));//,A
        V$EFC_Q6_pssm->push_back(PssmEntry(0,0,4,2));//,K
        V$EFC_Q6_pssm->push_back(PssmEntry(1,0,5,0));//,G
        V$EFC_Q6_pssm->push_back(PssmEntry(0,4,2,0));//,S
        V$EFC_Q6_pssm->push_back(PssmEntry(4,0,1,1));//,A
        V$EFC_Q6_pssm->push_back(PssmEntry(5,1,0,0));//,A
        V$EFC_Q6_pssm->push_back(PssmEntry(3,3,0,0));//,M
        V$EFC_Q6_iupac.clear();
        for (const char * cp = "MRWYRCYAKGSAAM"; *cp != '\0'; ++cp) {
            V$EFC_Q6_iupac.push_back(*cp);
        }
        assert(V$EFC_Q6_iupac.size() == V$EFC_Q6_pssm->size());
        biobase_params.push_back(
            BiobaseParams(
                "V$EFC_Q6",
                V$EFC_Q6_pssm,
                &V$EFC_Q6_iupac,
                test_rnd_1000_seq.begin(),
                test_rnd_1000_seq.end(),
                true,
                true,
                6.0f,
                Hit(0.853f, 335, true)));

        V$AP1_Q4_pssm.reset(new Pssm);
        V$AP1_Q4_pssm->push_back(PssmEntry(10,2,10,1));//,R
        V$AP1_Q4_pssm->push_back(PssmEntry(5,6,12,0));//,G
        V$AP1_Q4_pssm->push_back(PssmEntry(0,0,0,23));//,T
        V$AP1_Q4_pssm->push_back(PssmEntry(0,0,23,0));//,G
        V$AP1_Q4_pssm->push_back(PssmEntry(22,0,1,0));//,A
        V$AP1_Q4_pssm->push_back(PssmEntry(0,21,0,2));//,C
        V$AP1_Q4_pssm->push_back(PssmEntry(2,2,2,17));//,T
        V$AP1_Q4_pssm->push_back(PssmEntry(9,9,3,2));//,M
        V$AP1_Q4_pssm->push_back(PssmEntry(16,5,1,1));//,A
        V$AP1_Q4_pssm->push_back(PssmEntry(3,6,10,4));//,N
        V$AP1_Q4_pssm->push_back(PssmEntry(6,4,3,10));//,N
        V$AP1_Q4_iupac.clear();
        for (const char * cp = "RGTGACTMANN"; *cp != '\0'; ++cp) {
            V$AP1_Q4_iupac.push_back(*cp);
        }
        assert(V$AP1_Q4_iupac.size() == V$AP1_Q4_pssm->size());
        biobase_params.push_back(
            BiobaseParams(
                "V$AP1_Q4",
                V$AP1_Q4_pssm,
                &V$AP1_Q4_iupac,
                test_rnd_1000_seq.begin(),
                test_rnd_1000_seq.end(),
                true,
                true,
                3.0f,
                Hit(0.990f, 342)));

        already_done = true;
    }
}



