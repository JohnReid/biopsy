
#ifndef BIO_TSS_ESTIMATES_H_
#define BIO_TSS_ESTIMATES_H_

#include "bio/defs.h"
#include "bio/singleton.h"

#ifdef _MSC_VER
#pragma warning(push)
#pragma warning(disable : 4312)
#pragma warning(disable : 4311)
#endif // _MSC_VER

#include <boost/multi_index_container.hpp>
#include <boost/multi_index/member.hpp>
#include <boost/multi_index/ordered_index.hpp>

#ifdef _MSC_VER
#pragma warning(pop)
#endif //_MSC_VER

#include <iostream>

BIO_NS_START



struct TssEstimate
{
	enum Status
	{
		UNKNOWN_STATUS,		/**< Unknown status. */
		NO_CLONES,			/**< Couldn't find clones. */
		NO_CONCLUSION,		/**< Found clones but not interesting. */
		CONFIRMATION,		/**< Only confirms ensembl. */
		ALTERNATIVE,		/**< Found inward TSS. */
		EXONS_ONLY,			/**< No prediction for new TSS, just new exons. */
		PROOF				/**< New TSS. */
	};

	enum CoordSystem
	{
		UNKNOWN_COORDS,			/**< Unknown coord system. */
		NULL_COORDS,			/**< Null coord system. */
		CHROMOSOME_COORDS,		/**< Chromosome coord system. */
		SCAFFOLD_COORDS,		/**< Scaffold coord system. */
		SUPERCONTIG_COORDS,		/**< Supercontig coord system. */
	};

	enum CloneType
	{
		UNKNOWN_CLONE_TYPE,					/**< Unknown clone type. */
		NULL_CLONE_TYPE,					/**< Null clone type. */
		STANDARD_CLONE_TYPE,				/**< Standard clone type. */
		RIKEN_CLONE_TYPE,					/**< RIKEN clone type. */
		RIKEN_FRAGMENT_CLONE_TYPE,			/**< RIKEN fragment clone type. */
		ENSEMBL_TRANSCRIPT_CLONE_TYPE,		/**< Ensembl transcript clone type. */
	};

	struct Exon
	{
		typedef std::vector< Exon > vec_t;

		int start;
		int end;

		Exon( int start = 0, int end = 0 );

		Exon & parse( const std::string & line );

private:
		friend class boost::serialization::access;

		template< typename Archive >
		void serialize( Archive & ar, const unsigned int version )
		{
			ar
				& start
				& end;
		}
	};

	std::string gene;
	std::string transcript;
	std::string database;
	bool five_prime;
	Status status;
	CoordSystem coord_system;
	std::string seq_region;
	int abs_tss_pos;
	bool positive_strand;
	int rel_tss_pos;
	std::string seq_1st_fifty;
	CloneType clone_type;
	std::string clone_set;
	std::string clone_id;
	std::string gene_species;
	std::string clone_species;
	Exon::vec_t exons;

	TssEstimate();

	TssEstimate & parse( const std::string & line );

	static Status parse_status( const std::string & str );
	static CoordSystem parse_coord_system( const std::string & str );
	static CloneType parse_clone_type( const std::string & str );

private:
    friend class boost::serialization::access;

	template< typename Archive >
	void serialize( Archive & ar, const unsigned int version )
	{
		ar
			& gene
			& transcript
			& database
			& five_prime
			& status
			& coord_system
			& seq_region
			& abs_tss_pos
			& positive_strand
			& rel_tss_pos
			& seq_1st_fifty
			& clone_type
			& clone_set
			& clone_id
			& gene_species
			& clone_species
			& exons;
	}
};

std::ostream &
operator<<( std::ostream & os, TssEstimate::Status status );

std::ostream &
operator<<( std::ostream & os, TssEstimate::CoordSystem coord_system );

std::ostream &
operator<<( std::ostream & os, TssEstimate::CloneType clone_type );

std::ostream &
operator<<( std::ostream & os, const TssEstimate::Exon & exon );

std::ostream &
operator<<( std::ostream & os, const TssEstimate tss_estimate );




struct TssEstimates
	: Singleton< TssEstimates >
{
	void init_singleton();

	struct by_gene { };

	typedef boost::multi_index::multi_index_container<
		TssEstimate,
		boost::multi_index::indexed_by<
			// sort by less< string > on gene
			boost::multi_index::ordered_non_unique<
				boost::multi_index::tag< by_gene >,
				boost::multi_index::member<
					TssEstimate,
					std::string,
					&TssEstimate::gene
				>
			>
		> 
	> estimate_set_t;

	estimate_set_t estimates;

	estimate_set_t & operator()();
	const estimate_set_t & operator()() const;

	void parse();

private:
    friend class boost::serialization::access;

	template< typename Archive >
	void save( Archive & ar, const unsigned int version ) const
	{
		unsigned num_entries = estimates.size();
		ar & num_entries;
		for( estimate_set_t::const_iterator e = estimates.begin();
			estimates.end() != e;
			++e )
		{
			ar & *e;
		}
	}

	template< typename Archive >
	void load( Archive & ar, const unsigned int version )
	{
		unsigned num_entries;
		ar & num_entries;
		while( 0 != num_entries )
		{
			TssEstimate estimate;
			ar & estimate;
			estimates.insert( estimate );

			--num_entries;
		}
	}
    BOOST_SERIALIZATION_SPLIT_MEMBER()
};


BIO_NS_END

#endif //BIO_TSS_ESTIMATES_H_
