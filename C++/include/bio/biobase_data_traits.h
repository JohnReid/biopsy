#ifndef BIO_BIOBASE_DATA_TRAITS_H_
#define BIO_BIOBASE_DATA_TRAITS_H_

#include "bio/defs.h"
#include "bio/common.h"
#include "bio/matrix.h"
#include "bio/site.h"
#include "bio/factor.h"
#include "bio/fragment.h"
#include "bio/gene.h"
#include "bio/pathway.h"
#include "bio/molecule.h"
#include "bio/compel.h"
#include "bio/evidence.h"


BIO_NS_START


class MatrixParser;
class SiteParser;
class FactorParser;
class FragmentParser;
class GeneParser;
class PathwayParser;
class MoleculeParser;
class CompelParser;
class EvidenceParser;

struct BiobaseDb;

template <TransData type>
struct DataTraits ;

template<>
struct DataTraits<MATRIX_DATA>
{
	typedef Matrix entry_t;
	typedef MatrixParser parser_t;
	static std::string get_name();
	static std::string get_biobase_file();
	static std::string get_serialised_file();
	static std::string get_serialised_binary_file();
	static size_t get_num_data();
	static entry_t::map_t & get_map_in(BiobaseDb & db);
	static const entry_t::map_t & get_map_in(const BiobaseDb & db);
};

template<>
struct DataTraits<SITE_DATA>
{
	typedef Site entry_t;
	typedef SiteParser parser_t;
	static std::string get_name();
	static std::string get_biobase_file();
	static std::string get_serialised_file();
	static std::string get_serialised_binary_file();
	static size_t get_num_data();
	static entry_t::map_t & get_map_in(BiobaseDb & db);
	static const entry_t::map_t & get_map_in(const BiobaseDb & db);
};

template<>
struct DataTraits<FACTOR_DATA>
{
	typedef Factor entry_t;
	typedef FactorParser parser_t;
	static std::string get_name();
	static std::string get_biobase_file();
	static std::string get_serialised_file();
	static std::string get_serialised_binary_file();
	static size_t get_num_data();
	static entry_t::map_t & get_map_in(BiobaseDb & db);
	static const entry_t::map_t & get_map_in(const BiobaseDb & db);
};

template<>
struct DataTraits<FRAGMENT_DATA>
{
	typedef Fragment entry_t;
	typedef FragmentParser parser_t;
	static std::string get_name();
	static std::string get_biobase_file();
	static std::string get_serialised_file();
	static std::string get_serialised_binary_file();
	static size_t get_num_data();
	static entry_t::map_t & get_map_in(BiobaseDb & db);
	static const entry_t::map_t & get_map_in(const BiobaseDb & db);
};

template<>
struct DataTraits<GENE_DATA>
{
	typedef Gene entry_t;
	typedef GeneParser parser_t;
	static std::string get_name();
	static std::string get_biobase_file();
	static std::string get_serialised_file();
	static std::string get_serialised_binary_file();
	static size_t get_num_data();
	static entry_t::map_t & get_map_in(BiobaseDb & db);
	static const entry_t::map_t & get_map_in(const BiobaseDb & db);
};

template<>
struct DataTraits<COMPEL_DATA>
{
	typedef Compel entry_t;
	typedef CompelParser parser_t;
	static std::string get_name();
	static std::string get_biobase_file();
	static std::string get_serialised_file();
	static std::string get_serialised_binary_file();
	static size_t get_num_data();
	static entry_t::map_t & get_map_in(BiobaseDb & db);
	static const entry_t::map_t & get_map_in(const BiobaseDb & db);
};

template<>
struct DataTraits<EVIDENCE_DATA>
{
	typedef Evidence entry_t;
	typedef EvidenceParser parser_t;
	static std::string get_name();
	static std::string get_biobase_file();
	static std::string get_serialised_file();
	static std::string get_serialised_binary_file();
	static size_t get_num_data();
	static entry_t::map_t & get_map_in(BiobaseDb & db);
	static const entry_t::map_t & get_map_in(const BiobaseDb & db);
};

template<>
struct DataTraits<PATHWAY_DATA>
{
	typedef Pathway entry_t;
	typedef PathwayParser parser_t;
	static std::string get_name();
	static std::string get_biobase_file();
	static std::string get_serialised_file();
	static std::string get_serialised_binary_file();
	static size_t get_num_data();
	static entry_t::map_t & get_map_in(BiobaseDb & db);
	static const entry_t::map_t & get_map_in(const BiobaseDb & db);
};

template<>
struct DataTraits<MOLECULE_DATA>
{
	typedef Molecule entry_t;
	typedef MoleculeParser parser_t;
	static std::string get_name();
	static std::string get_biobase_file();
	static std::string get_serialised_file();
	static std::string get_serialised_binary_file();
	static size_t get_num_data();
	static entry_t::map_t & get_map_in(BiobaseDb & db);
	static const entry_t::map_t & get_map_in(const BiobaseDb & db);
};

BIO_NS_END

#endif //BIO_BIOBASE_DATA_TRAITS_H_
