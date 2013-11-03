/** Copyright John Reid 2012
 *
 * \file
 * \brief Adapts TRANSFAC data structures using boost::fusion.
 */


#ifndef BIO_JR_30AUG2011_TRANSFAC_FUSION_H_
#define BIO_JR_30AUG2011_TRANSFAC_FUSION_H_


#include "bio/defs.h"
#include "bio/common.h"
#include "bio/compel.h"
#include "bio/evidence.h"
#include "bio/factor.h"
#include "bio/fragment.h"
#include "bio/gene.h"
#include "bio/matrix.h"
#include "bio/molecule.h"
#include "bio/pathway.h"
#include "bio/site.h"

#include <boost/fusion/include/adapted.hpp>


// Keys for associative access
BIO_NS_START
namespace spirit {
namespace keys {

// keys for FactorLink
struct link_tag;
struct name_tag;
struct quality_tag;
struct species_tag;
struct cellular_source_tag;
struct sites_included_tag;

// keys for biobase entries
struct accession_number_tag;
struct id_tag;
struct name_tag;
struct matrix_basis_tag;
struct number_of_sites_tag;
struct description_tag;
struct factor_links_tag;
struct pssm_tag;
struct consensus_matrix_tag;
struct align_descs_tag;
struct sequence_tag;
struct db_refs_tag;
struct ref_point_tag;
struct start_tag;
struct end_tag;
struct synonyms_tag;
struct taxonomies_tag;
struct matrices_tag;
struct gene_tag;
struct type_tag;
struct subunits_tag;
struct complexes_tag;
struct subfamilies_tag;
struct superfamilies_tag;

} //namespace keys
} //namespace spirit
BIO_NS_END




// We need to tell fusion about our struct(s)
// to make them first-class fusion citizen(s). This has to
// be in global scope.
BOOST_FUSION_ADAPT_STRUCT(
    BIO_NS::TableLink,
    (BIO_NS::TransData, table_id)
    (int, entry_idx)
)
BOOST_FUSION_ADAPT_ASSOC_STRUCT(
    BIO_NS::FactorLink,
	(BIO_NS::TableLink, link, BIO_NS::spirit::keys::link_tag)
	(std::string, name, BIO_NS::spirit::keys::name_tag)
	(int, quality, BIO_NS::spirit::keys::quality_tag)
	(std::string, cellular_source, BIO_NS::spirit::keys::cellular_source_tag)
	(bool, sites_included, BIO_NS::spirit::keys::sites_included_tag)
	(std::vector< std::string >, species, BIO_NS::spirit::keys::species_tag)
)
BOOST_FUSION_ADAPT_STRUCT(
    BIO_NS::Identifier,
	(std::string, species_group)
	(std::string, factor)
	(std::string, discriminating_extension)
)
BOOST_FUSION_ADAPT_STRUCT(
    BIO_NS::db_ref,
	(BIO_NS::Database, db)
	(std::string, table)
	(int, acc)
)
BOOST_FUSION_ADAPT_STRUCT(
    BIO_NS::AlignDesc,
	(std::string, sequence)
	(BIO_NS::TableLink, site)
	(BIO_NS::TableLink, secondary_site)
	(int, start)
	(int, length)
	(std::vector<int>, gaps)
	(bool, positive_orientation)
)
BOOST_FUSION_ADAPT_ADT(
    BIO_NS::PssmEntry,
    (BIO_NS::float_t, BIO_NS::float_t, obj.get_count('a'), /**/)
    (BIO_NS::float_t, BIO_NS::float_t, obj.get_count('c'), /**/)
    (BIO_NS::float_t, BIO_NS::float_t, obj.get_count('g'), /**/)
    (BIO_NS::float_t, BIO_NS::float_t, obj.get_count('t'), /**/)
)
BOOST_FUSION_DEFINE_STRUCT(
	(BIO_NS),
	PssmAndConsensus,
	(BIO_NS::Pssm, pssm)
	(BIO_NS::ConsensusMatrix, consensus)
)
BOOST_FUSION_ADAPT_ASSOC_STRUCT(
	BIO_NS::Matrix,
	(BIO_NS::TableLink,       accession_number,   BIO_NS::spirit::keys::accession_number_tag)
	(BIO_NS::Identifier,      id,                 BIO_NS::spirit::keys::id_tag)
	(std::string,             factor_name,        BIO_NS::spirit::keys::name_tag)
	(std::string,             matrix_basis,       BIO_NS::spirit::keys::matrix_basis_tag)
	(int,                     number_of_sites,    BIO_NS::spirit::keys::number_of_sites_tag)
	(std::string,             description,        BIO_NS::spirit::keys::description_tag)
	(BIO_NS::FactorLinkList,  factor_links,       BIO_NS::spirit::keys::factor_links_tag)
	(BIO_NS::AlignDescList,   align_descs,        BIO_NS::spirit::keys::align_descs_tag)
	(BIO_NS::Pssm,            pssm,               BIO_NS::spirit::keys::pssm_tag)
	(BIO_NS::ConsensusMatrix, consensus_matrix,   BIO_NS::spirit::keys::consensus_matrix_tag)
)
BOOST_FUSION_ADAPT_ASSOC_STRUCT(
	BIO_NS::Site,
	(BIO_NS::TableLink,       accession_number,   BIO_NS::spirit::keys::accession_number_tag)
	(BIO_NS::Identifier,      id,                 BIO_NS::spirit::keys::id_tag)
	(BIO_NS::seq_t,           sequence,           BIO_NS::spirit::keys::sequence_tag)
	(std::string,             description,        BIO_NS::spirit::keys::description_tag)
	(BIO_NS::FactorLinkList,  factor_links,       BIO_NS::spirit::keys::factor_links_tag)
	(BIO_NS::DatabaseRefVec,  database_refs,      BIO_NS::spirit::keys::db_refs_tag)
	(std::string,             reference_point,    BIO_NS::spirit::keys::ref_point_tag)
	(int,                     start_position,     BIO_NS::spirit::keys::start_tag)
	(int,                     end_position,       BIO_NS::spirit::keys::end_tag)
)
BOOST_FUSION_ADAPT_ASSOC_STRUCT(
    BIO_NS::Factor,
    (BIO_NS::TableLink,               accession_number,   BIO_NS::spirit::keys::accession_number_tag)
    (std::string,                     name,               BIO_NS::spirit::keys::name_tag)
    (BIO_NS::Factor::type,            _type,              BIO_NS::spirit::keys::type_tag)
    (BIO_NS::TableLink,               gene,               BIO_NS::spirit::keys::gene_tag)
    (BIO_NS::TableLinkVec,            matrices,           BIO_NS::spirit::keys::matrices_tag)
    (BIO_NS::Factor::synonym_set_t,   synonyms,           BIO_NS::spirit::keys::synonyms_tag)
    (BIO_NS::DatabaseRefVec,          database_refs,      BIO_NS::spirit::keys::db_refs_tag)
    (BIO_NS::Factor::taxonomies_t,    taxonomies,         BIO_NS::spirit::keys::taxonomies_tag)
    (BIO_NS::TableLinkVec,            subunits,           BIO_NS::spirit::keys::subunits_tag)
    (BIO_NS::TableLinkVec,            complexes,          BIO_NS::spirit::keys::complexes_tag)
    (BIO_NS::TableLinkVec,            sub_families,       BIO_NS::spirit::keys::subfamilies_tag)
    (BIO_NS::TableLinkVec,            super_families,     BIO_NS::spirit::keys::superfamilies_tag)
)
BOOST_FUSION_ADAPT_STRUCT(
    BIO_NS::Fragment,
    (BIO_NS::TableLink,               accession_number)
    (BIO_NS::FactorLinkList,          factor_links)
    (BIO_NS::TableLinkVec,            genes)
    (BIO_NS::seq_t,                   sequence)
)
BOOST_FUSION_ADAPT_STRUCT(
    BIO_NS::Gene,
    (BIO_NS::TableLink,               accession_number)
    (std::string,                     name)
    (std::string,                     species)
    (BIO_NS::DatabaseRefVec,          database_refs)
)
BOOST_FUSION_ADAPT_STRUCT(
    BIO_NS::Compel,
    (BIO_NS::TableLink,               accession_number)
    (BIO_NS::Identifier,              id)
    (BIO_NS::TableLink,               gene)
    (BIO_NS::seq_t,                   sequence)
    (int,                             begin)
    (int,                             end)
    (BIO_NS::CompelBindingSite::vec,  binding_sites)
    (BIO_NS::CompelType,              type)
    (BIO_NS::DatabaseRefVec,          database_refs)
    (std::string,                     comment)
    (BIO_NS::TableLinkVec,            evidences)
)
BOOST_FUSION_ADAPT_STRUCT(
    BIO_NS::Evidence,
    (BIO_NS::TableLink,               accession_number)
    (BIO_NS::TableLink,               composite_element)
    (BIO_NS::DatabaseRefVec,          database_refs)
)
BOOST_FUSION_ADAPT_STRUCT(
    BIO_NS::Molecule,
    (BIO_NS::TableLink,               accession_number)
    (BIO_NS::TableLinkVec,            pathways)
    (BIO_NS::TableLinkVec,            super_families)
    (BIO_NS::DatabaseRefVec,          database_refs)
)
BOOST_FUSION_ADAPT_STRUCT(
    BIO_NS::Pathway,
    (BIO_NS::TableLink,               accession_number)
    (BIO_NS::PathwayType,             pathway_type)
    (BIO_NS::TableLinkVec,            super_families)
    (std::string,                     name)
)




#endif //BIO_JR_30AUG2011_TRANSFAC_FUSION_H_
