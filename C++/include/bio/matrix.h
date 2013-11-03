#ifndef BIOBASE_MATRIX_H_
#define BIOBASE_MATRIX_H_

#include "bio/defs.h"
#include "bio/common.h"
#include "bio/sequence.h"
#include "bio/pssm.h"
#include "bio/pathway_associations.h"

#include <boost/serialization/access.hpp>
#include <boost/shared_ptr.hpp>

#include <string>
#include <list>
#include <vector>
#include <map>
#include <algorithm>

#define ANALOGUE_SITE_EQUIVALENT 50

BIO_NS_START





typedef char ConsensusEntry;
typedef std::vector<ConsensusEntry> ConsensusMatrix;





/** From the TRANSFAC documentation:
BS  Sequence segment (including gaps); Site link; Start; Length; Gaps; Orientation.

Examples:
BS  TGTTTGTCAAT; R08890; 9; 11;; p.             derived from TRANSFAC site    R08890
BS  TATTTACTTTC; C00278,  s00475; 2; 11;; n.    derived from TRANSCompel site s00475
                                                         in composite element C00278

The positions given refer to the sequence as shown in the site table.

For re-construction of alignment from site sequence
1. introduce gaps at the indicated positions
2. determine matrix start
3. determine matrix length
4. orientation: p (positive), i.e. as in site entry / n (negative) = reverse complement


Examples:

BS    AATCCGGAAACGATG; R05296; 1; 17; 1, 2; p
    -----------------
                                                 1111111111
                                        1234567890123456789
Sequence in site entry:                 AATCCGGAAACGATGCG
1. positions of gaps: 1, 2              --AATCCGGAAACGATGCG
2. matrix start: 1                      --AATCCGGAAACGATGCG
3. matrix length: 17                    --AATCCGGAAACGATG
4. orientation: p



BS  CCTGGTCA ATGGGTCATG; R11613; 9; 19; 17; p
    -------------------
                                                 111111111122222222
                                        123456789012345678901234567
Sequence in site entry:                 GGGTCACTCCTGGTCAATGGGTCATG
1. positions of gaps: 17                GGGTCACTCCTGGTCA-ATGGGTCATG
2. matrix start: 9                              CCTGGTCA-ATGGGTCATG
3. matrix length: 19                            CCTGGTCA-ATGGGTCATG
4. orientation: p



BS    AACTACCGG; R08780; 1; 11; 10, 11; n.
    -----------
                                                 11
                                        12345678901
Sequence in site entry:                 CCGGTAGTT
1. positions of gaps: 10, 11            CCGGTAGTT--
2. matrix start: 1                      CCGGTAGTT--
3. matrix length: 11                    CCGGTAGTT--
4. orientation: n                       --AACTACCGG



BS  CATTACAAAATC; R00097; -1; 12;; n.
    ------------
                                                 11
                                      -112345678901
Sequence in site entry:                 ATTTTGTAAT
1. positions of gaps: -                 ATTTTGTAAT
2. matrix start: -1                    GATTTTGTAAT
3. matrix length: 12                   GATTTTGTAATG
4. orientation: n                      CATTACAAAATC
(Flanking positions were retrieved from linked EMBL entry.)
*/
struct AlignDesc {
	std::string sequence;
	TableLink site;
	TableLink secondary_site;
	int start;
	int length;
	std::vector<int> gaps;
	bool positive_orientation;

	AlignDesc();

private:
    friend class boost::serialization::access;
    // When the class Archive corresponds to an output archive, the
    // & operator is defined similar to <<.  Likewise, when the class Archive
    // is a type of input archive the & operator is defined similar to >>.
    template<class Archive>
    void serialize(Archive & ar, const unsigned int version)
    {
        ar & sequence;
        ar & site;
        ar & secondary_site;
		ar & start;
        ar & length;
        ar & gaps;
		ar & positive_orientation;
    }		 
};
typedef boost::shared_ptr<AlignDesc> AlignDescPtr;
typedef std::list<AlignDescPtr> AlignDescList;

/** The information from one entry in the TRANSFAC matrix table.

The MATRIX table contains nucleotide distribution matrices of aligned binding sequences ( Statistics   Statistics). These sequences may have been obtained by in vitro selection studies or may be compiled sites of genes. The source is appropriately indicated. The matrix entries have an identifier that indicates one of six groups of biological species (V$, vertebrates; I$, insects; P$, plants; F$, fungi; N$ nematodes; B$, bacteria), followed by an acronym for the factor the matrix refers to, and a consecutive number discriminating between different matrices for the same factor. Thus, V$OCT1_02 indicates the second matrix for vertebral Oct-1 factor. Instead of a consecutive number, the identifier of those matrices which have been generated from TRANSFAC® SITE entries, end up with an abbreviation of the least quality of the sites used to construct the matrix. For example, V$CREB_Q2 is a matrix constructed of CREB binding sites of quality 2 or better. Finally, a matrix with an identifier like V$AP1_C has been derived from a "consensus description" constructed with the aid of ConsIndex (Frech et al., Nucleic Acids Res. 21:1655-1664, 1993).

The matrix area gives the nucleotide frequencies observed in aligned binding sites of the corresponding transcription factor (or, more general, in aligned sites of the described function); an additional column depicts the IUPAC string consensus derived from the matrix according to the following rules (adapted from Cavener, Nucleic Acids Res. 15:1353-1361, 1987):
Rule 1: A single nucleotide is shown if its frequency is at least 50% and at least twice as high as the second most frequent nucleotide.
Rule 2: A double-degenerate code indicates that the corresponding two nucleotides occur in more than 75% of the underlying sequences and rule 1 does not apply.
Rule 3: Usage of triple-degenerate codes is restricted to those positions where one of the nucleotides did not show up at all in the sequence set and none of the afore mentioned rules applies.
Rule 4: All other frequency distributions are represented by the letter "N".
*/
struct Matrix : BiobaseTablePssmEntry
{
	Matrix() {
		number_of_sites = 0;
	};
	typedef boost::shared_ptr<Matrix> ptr_t;
	typedef std::map<TableLink, ptr_t> map_t;

	TableLink accession_number;

	Identifier id;

	std::string factor_name;
	std::string matrix_basis;
	int number_of_sites;
	std::string description;
	FactorLinkList factor_links;

	Pssm pssm;
	ConsensusMatrix consensus_matrix;

	AlignDescList align_descs;

	TableLink get_link() const { return accession_number; }
	std::string get_name() const { return id.get_text(); }
	std::string get_description() const {
		std::stringstream stream;
		stream << factor_name << ": " << description << ": ";
		stream << get_factors();
		return stream.str();
	}
	size_t get_size() const { return pssm.size(); }
	TableLink get_most_significant_pathway() const
	{
		return PathwayAssociations::singleton().get_most_significant_pathway_for( accession_number );
	};
	
	const FactorLinkList & get_factors() const { return factor_links; }

private:
    friend class boost::serialization::access;
    // When the class Archive corresponds to an output archive, the
    // & operator is defined similar to <<.  Likewise, when the class Archive
    // is a type of input archive the & operator is defined similar to >>.
    template<class Archive>
    void serialize(Archive & ar, const unsigned int version)
    {
        ar & accession_number;
        ar & id;
        ar & factor_name;
		ar & description;
        ar & factor_links;
        ar & pssm;
		ar & consensus_matrix;
		ar & align_descs;
		ar & number_of_sites;
    }		 
};




BIO_NS_END


#endif //BIOBASE_MATRIX_H_



