#ifndef BIO_BIOBASE_DB_H_
#define BIO_BIOBASE_DB_H_

#include "bio/defs.h"
#include "bio/matrix.h"
#include "bio/site.h"
#include "bio/factor.h"
#include "bio/compel.h"
#include "bio/evidence.h"
#include "bio/pathway.h"
#include "bio/fragment.h"
#include "bio/gene.h"
#include "bio/molecule.h"
#include "bio/singleton.h"

#include <boost/shared_ptr.hpp>

BIO_NS_START


/** This is specialised for the various types. */
template <TransData type>
struct DataTraits {
};



/** Contains all the tables in biobase. */
struct BiobaseDb
	: Singleton< BiobaseDb >
{
	typedef boost::shared_ptr<BiobaseDb> ptr_t;

protected:
	//transfac
	Matrix::map_t matrices;
	Site::map_t sites;
	Factor::map_t factors;
	Fragment::map_t fragments;
	Gene::map_t genes;

	//transcompel
	Compel::map_t compels;
	Evidence::map_t evidences;

	//transpath
	Pathway::map_t pathways;
	Molecule::map_t molecules;

public:
	BiobaseTableEntry * get_entry(const TableLink & link) const;
	BiobaseTablePssmEntry * get_pssm_entry(const TableLink & link) const;

	Matrix::map_t & get_matrices();
	const Matrix::map_t & get_matrices() const;

	Site::map_t & get_sites();
	const Site::map_t & get_sites() const;

	Factor::map_t & get_factors();
	const Factor::map_t & get_factors() const;

	Fragment::map_t & get_fragments();
	const Fragment::map_t & get_fragments() const;

	Gene::map_t & get_genes();
	const Gene::map_t & get_genes() const;

	Compel::map_t & get_compels();
	const Compel::map_t & get_compels() const;

	Evidence::map_t & get_evidences();
	const Evidence::map_t & get_evidences() const;
	Pathway::map_t & get_pathways();
	const Pathway::map_t & get_pathways() const;

	Molecule::map_t & get_molecules();
	const Molecule::map_t & get_molecules() const;


	/** Make sure all tables are loaded. */
	void load_all() const;


	template <TransData data_type>
	typename DataTraits<data_type>::entry_t *
	get_entry(const TableLink & link) const
	{
		if (link.table_id != data_type)
		{
			throw std::logic_error( "Wrong table link data type" );
		}

		typedef DataTraits<data_type> traits_t;

		typename traits_t::entry_t::map_t::const_iterator i = traits_t::get_map_in(*this).find(link);
		if (traits_t::get_map_in(*this).end() == i)
		{
			return 0; //null pointer for not in biobase_db
		}
		return i->second.get();
	}

	template < TransData data_type >
	typename DataTraits< data_type >::entry_t *
	get_entry(unsigned accession_number) const
	{
		typedef DataTraits< data_type > traits_t;

		typename traits_t::entry_t::map_t::const_iterator i = traits_t::get_map_in(*this).find(TableLink(data_type, accession_number));
		if (traits_t::get_map_in(*this).end() == i)
		{
			return 0; //null pointer for not in biobase_db
		}
		return i->second.get();
	}

	//sometimes we need the shared ptr
	template <TransData data_type>
	typename DataTraits<data_type>::entry_t::ptr_t
	get_entry_ptr(const TableLink & link) const
	{
		if (link.table_id != data_type)
		{
			throw std::logic_error( "Wrong table link data type" );
		}

		typedef DataTraits<data_type> traits_t;

		typename traits_t::entry_t::map_t::const_iterator i = traits_t::get_map_in(*this).find(link);
		if (traits_t::get_map_in(*this).end() == i)
		{
			return typename DataTraits<data_type>::entry_t::ptr_t(); //null pointer for not in biobase_db
		}
		return i->second;
	}
};





BIO_NS_END


#endif //BIO_BIOBASE_DB_H_

