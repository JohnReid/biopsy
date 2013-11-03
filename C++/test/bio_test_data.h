/**
@file

Copyright John Reid 2006, 2007, 2013
*/

#ifndef BIO_TEST_DATA_H_
#define BIO_TEST_DATA_H_


#include "bio/sequence.h"
#include "bio/biobase_match.h"
#include "bio/matrix.h"
#include <bio/remo.h>

#include <vector>
#include <string>

//data we use as parameter
extern const BIO_NS::seq_t test_rnd_1000_seq;
extern const BIO_NS::seq_t acgtacgt_seq;
extern const BIO_NS::seq_t atcgatcg_seq; //complement is tagctagc - reversed complement is cgatcgat
extern BIO_NS::Pssm::ptr_t acgt_pssm;
extern BIO_NS::Pssm::ptr_t atcg_pssm;
extern BIO_NS::Pssm::ptr_t V$DEAF1_01_pssm;
extern BIO_NS::Pssm::ptr_t V$IPF1_Q4_01_pssm;
extern BIO_NS::Pssm::ptr_t V$EFC_Q6_pssm;
extern BIO_NS::Pssm::ptr_t V$AP1_Q4_pssm;
extern BIO_NS::IupacSeq acgt_iupac;
extern BIO_NS::IupacSeq atcg_iupac;
extern BIO_NS::IupacSeq V$DEAF1_01_iupac;
extern BIO_NS::IupacSeq V$IPF1_Q4_01_iupac;
extern BIO_NS::IupacSeq V$EFC_Q6_iupac;
extern BIO_NS::IupacSeq V$AP1_Q4_iupac;

struct BiobaseParams {
	BiobaseParams(
			const std::string & name,
			BIO_NS::Pssm::ptr_t psmm,					
			BIO_NS::IupacSeq * iupac,					
			std::string::const_iterator seq_begin,
			std::string::const_iterator seq_end,	
			bool match_complement,			
			bool run_over_range,		
			BIO_NS::float_t iupac_tolerance,
			BIO_NS::Hit target_match)					
		: name(name)
		, pssm(psmm)
		, iupac(iupac)
		, seq_begin(seq_begin)
		, seq_end(seq_end)
		, match_complement(match_complement)
		, run_over_range(run_over_range)
		, iupac_tolerance(iupac_tolerance)
		, target_match(target_match)
	{ }
	std::string name;						//name of parameterised test
	BIO_NS::Pssm::ptr_t pssm;				//pssm matrix
	BIO_NS::IupacSeq * iupac;				//iupac sequence
	std::string::const_iterator seq_begin;	//beginning of sequence
	std::string::const_iterator seq_end;	//end of sequence
	bool match_complement;					//match the reversed complement
	bool run_over_range;					//run over the whole sequence or just match the matrix once
	BIO_NS::float_t iupac_tolerance;		//% that the iupac score can differ from the pssm score
	BIO_NS::Hit target_match;		//target score (and position for pssm)
};

//vector of parameters
extern std::vector<BiobaseParams> biobase_params;


void build_biobase_params();
void ensure_biobase_params_built();
void ensure_test_remos_built();





#endif //BIO_TEST_DATA_H_



