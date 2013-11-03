/* Copyright John Reid 2007
*/

#include "bio-pch.h"


#include "bio/defs.h"


#include "bio/biobase_data_traits.h"
#include "bio/environment.h"
#include "bio/biobase_db.h"


#ifdef WIN32
# define DIR_SEP "\\"
#else //WIN32
# define DIR_SEP "/"
#endif //WIN32



BIO_NS_START


std::string DataTraits<MATRIX_DATA>::get_name() { return "matrix"; }
std::string DataTraits<MATRIX_DATA>::get_biobase_file() { return BioEnvironment::singleton().get_transfac_dir() + DIR_SEP + get_name() + ".dat"; }
std::string DataTraits<MATRIX_DATA>::get_serialised_file() { return BioEnvironment::singleton().get_serialised_dir() + DIR_SEP + "biobase_" + get_name() + ".txt"; }
std::string DataTraits<MATRIX_DATA>::get_serialised_binary_file() { return BioEnvironment::singleton().get_serialised_dir() + DIR_SEP + "biobase_" + get_name() + ".bin"; }
size_t DataTraits<MATRIX_DATA>::get_num_data() { return 795; }
DataTraits<MATRIX_DATA>::entry_t::map_t & DataTraits<MATRIX_DATA>::get_map_in(BiobaseDb & db) { return db.get_matrices(); }
const DataTraits<MATRIX_DATA>::entry_t::map_t & DataTraits<MATRIX_DATA>::get_map_in(const BiobaseDb & db) { return db.get_matrices(); }

std::string DataTraits<SITE_DATA>::get_name() { return "site"; }
std::string DataTraits<SITE_DATA>::get_biobase_file() { return BioEnvironment::singleton().get_transfac_dir() + DIR_SEP + get_name() + ".dat"; }
std::string DataTraits<SITE_DATA>::get_serialised_file() { return BioEnvironment::singleton().get_serialised_dir() + DIR_SEP + "biobase_" + get_name() + ".txt"; }
std::string DataTraits<SITE_DATA>::get_serialised_binary_file() { return BioEnvironment::singleton().get_serialised_dir() + DIR_SEP + "biobase_" + get_name() + ".bin"; }
size_t DataTraits<SITE_DATA>::get_num_data() { return 16992; }
DataTraits<SITE_DATA>::entry_t::map_t & DataTraits<SITE_DATA>::get_map_in(BiobaseDb & db) { return db.get_sites(); }
const DataTraits<SITE_DATA>::entry_t::map_t & DataTraits<SITE_DATA>::get_map_in(const BiobaseDb & db) { return db.get_sites(); }

std::string DataTraits<FACTOR_DATA>::get_name() { return "factor"; }
std::string DataTraits<FACTOR_DATA>::get_biobase_file() { return BioEnvironment::singleton().get_transfac_dir() + DIR_SEP + get_name() + ".dat"; }
std::string DataTraits<FACTOR_DATA>::get_serialised_file() { return BioEnvironment::singleton().get_serialised_dir() + DIR_SEP + "biobase_" + get_name() + ".txt"; }
std::string DataTraits<FACTOR_DATA>::get_serialised_binary_file() { return BioEnvironment::singleton().get_serialised_dir() + DIR_SEP + "biobase_" + get_name() + ".bin"; }
size_t DataTraits<FACTOR_DATA>::get_num_data() { return 3314; }
DataTraits<FACTOR_DATA>::entry_t::map_t & DataTraits<FACTOR_DATA>::get_map_in(BiobaseDb & db) { return db.get_factors(); }
const DataTraits<FACTOR_DATA>::entry_t::map_t & DataTraits<FACTOR_DATA>::get_map_in(const BiobaseDb & db) { return db.get_factors(); }

std::string DataTraits<FRAGMENT_DATA>::get_name() { return "fragment"; }
std::string DataTraits<FRAGMENT_DATA>::get_biobase_file() { return BioEnvironment::singleton().get_transfac_dir() + DIR_SEP + get_name() + ".dat"; }
std::string DataTraits<FRAGMENT_DATA>::get_serialised_file() { return BioEnvironment::singleton().get_serialised_dir() + DIR_SEP + "biobase_" + get_name() + ".txt"; }
std::string DataTraits<FRAGMENT_DATA>::get_serialised_binary_file() { return BioEnvironment::singleton().get_serialised_dir() + DIR_SEP + "biobase_" + get_name() + ".bin"; }
size_t DataTraits<FRAGMENT_DATA>::get_num_data() { return 3841; }
DataTraits<FRAGMENT_DATA>::entry_t::map_t & DataTraits<FRAGMENT_DATA>::get_map_in(BiobaseDb & db) { return db.get_fragments(); }
const DataTraits<FRAGMENT_DATA>::entry_t::map_t & DataTraits<FRAGMENT_DATA>::get_map_in(const BiobaseDb & db) { return db.get_fragments(); }

std::string DataTraits<GENE_DATA>::get_name() { return "gene"; }
std::string DataTraits<GENE_DATA>::get_biobase_file() { return BioEnvironment::singleton().get_transfac_dir() + DIR_SEP + get_name() + ".dat"; }
std::string DataTraits<GENE_DATA>::get_serialised_file() { return BioEnvironment::singleton().get_serialised_dir() + DIR_SEP + "biobase_" + get_name() + ".txt"; }
std::string DataTraits<GENE_DATA>::get_serialised_binary_file() { return BioEnvironment::singleton().get_serialised_dir() + DIR_SEP + "biobase_" + get_name() + ".bin"; }
size_t DataTraits<GENE_DATA>::get_num_data() { return 1533; }
DataTraits<GENE_DATA>::entry_t::map_t & DataTraits<GENE_DATA>::get_map_in(BiobaseDb & db) { return db.get_genes(); }
const DataTraits<GENE_DATA>::entry_t::map_t & DataTraits<GENE_DATA>::get_map_in(const BiobaseDb & db) { return db.get_genes(); }

std::string DataTraits<PATHWAY_DATA>::get_name() { return "pathway"; }
std::string DataTraits<PATHWAY_DATA>::get_biobase_file() { return BioEnvironment::singleton().get_transpath_dir() + DIR_SEP + get_name() + ".dat"; }
std::string DataTraits<PATHWAY_DATA>::get_serialised_file() { return BioEnvironment::singleton().get_serialised_dir() + DIR_SEP + "biobase_" + get_name() + ".txt"; }
std::string DataTraits<PATHWAY_DATA>::get_serialised_binary_file() { return BioEnvironment::singleton().get_serialised_dir() + DIR_SEP + "biobase_" + get_name() + ".bin"; }
size_t DataTraits<PATHWAY_DATA>::get_num_data() { return 1002; }
DataTraits<PATHWAY_DATA>::entry_t::map_t & DataTraits<PATHWAY_DATA>::get_map_in(BiobaseDb & db) { return db.get_pathways(); }
const DataTraits<PATHWAY_DATA>::entry_t::map_t & DataTraits<PATHWAY_DATA>::get_map_in(const BiobaseDb & db) { return db.get_pathways(); }

std::string DataTraits<MOLECULE_DATA>::get_name() { return "molecule"; }
std::string DataTraits<MOLECULE_DATA>::get_biobase_file() { return BioEnvironment::singleton().get_transpath_dir() + DIR_SEP + get_name() + ".dat"; }
std::string DataTraits<MOLECULE_DATA>::get_serialised_file() { return BioEnvironment::singleton().get_serialised_dir() + DIR_SEP + "biobase_" + get_name() + ".txt"; }
std::string DataTraits<MOLECULE_DATA>::get_serialised_binary_file() { return BioEnvironment::singleton().get_serialised_dir() + DIR_SEP + "biobase_" + get_name() + ".bin"; }
size_t DataTraits<MOLECULE_DATA>::get_num_data() { return 52040; }
DataTraits<MOLECULE_DATA>::entry_t::map_t & DataTraits<MOLECULE_DATA>::get_map_in(BiobaseDb & db) { return db.get_molecules(); }
const DataTraits<MOLECULE_DATA>::entry_t::map_t & DataTraits<MOLECULE_DATA>::get_map_in(const BiobaseDb & db) { return db.get_molecules(); }

std::string DataTraits<COMPEL_DATA>::get_name() { return "compel"; }
std::string DataTraits<COMPEL_DATA>::get_biobase_file() { return BioEnvironment::singleton().get_transcompel_dir() + DIR_SEP + get_name() + ".dat"; }
std::string DataTraits<COMPEL_DATA>::get_serialised_file() { return BioEnvironment::singleton().get_serialised_dir() + DIR_SEP + "biobase_" + get_name() + ".txt"; }
std::string DataTraits<COMPEL_DATA>::get_serialised_binary_file() { return BioEnvironment::singleton().get_serialised_dir() + DIR_SEP + "biobase_" + get_name() + ".bin"; }
size_t DataTraits<COMPEL_DATA>::get_num_data() { return 421; }
DataTraits<COMPEL_DATA>::entry_t::map_t & DataTraits<COMPEL_DATA>::get_map_in(BiobaseDb & db) { return db.get_compels(); }
const DataTraits<COMPEL_DATA>::entry_t::map_t & DataTraits<COMPEL_DATA>::get_map_in(const BiobaseDb & db) { return db.get_compels(); }

std::string DataTraits<EVIDENCE_DATA>::get_name() { return "evidence"; }
std::string DataTraits<EVIDENCE_DATA>::get_biobase_file() { return BioEnvironment::singleton().get_transcompel_dir() + DIR_SEP + get_name() + ".dat"; }
std::string DataTraits<EVIDENCE_DATA>::get_serialised_file() { return BioEnvironment::singleton().get_serialised_dir() + DIR_SEP + "biobase_" + get_name() + ".txt"; }
std::string DataTraits<EVIDENCE_DATA>::get_serialised_binary_file() { return BioEnvironment::singleton().get_serialised_dir() + DIR_SEP + "biobase_" + get_name() + ".bin"; }
size_t DataTraits<EVIDENCE_DATA>::get_num_data() { return 1460; }
DataTraits<EVIDENCE_DATA>::entry_t::map_t & DataTraits<EVIDENCE_DATA>::get_map_in(BiobaseDb & db) { return db.get_evidences(); }
const DataTraits<EVIDENCE_DATA>::entry_t::map_t & DataTraits<EVIDENCE_DATA>::get_map_in(const BiobaseDb & db) { return db.get_evidences(); }



BIO_NS_END


