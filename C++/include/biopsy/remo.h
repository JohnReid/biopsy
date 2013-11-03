/**
@file

Copyright John Reid 2006

*/

#ifndef BIOPSY_REMO_H_
#define BIOPSY_REMO_H_

#ifdef _MSC_VER
# pragma once
#endif //_MSC_VER

#include "biopsy/defs.h"



namespace biopsy {
namespace remo {

	
typedef std::string gene;
typedef std::string species;
typedef int position;


/**
Where w.r.t. a gene a remo is located.
*/
enum region
{
	region_upstream,
	region_downstream,
	region_gene,
	region_undefined
};
std::ostream &
operator<<( std::ostream &, region r );


/**
Defines a location w.r.t. a gene.
*/
struct location
	: boost::equality_comparable< location >
	, boost::less_than_comparable< location >

{
	typedef std::vector< location > list;
	typedef boost::shared_ptr< list > list_ptr;

	int _start;
	int _end;

	location( );
	location( int start, int end );

	bool operator==( const location & rhs ) const;
	bool operator<( const location & rhs ) const;

	std::string str() const;
};
std::ostream &
operator<<( std::ostream &, const location & l );

/**
An ensembl id.
*/
struct ensembl_id
	: boost::equality_comparable< ensembl_id >
	, boost::less_than_comparable< ensembl_id >
{
	typedef std::vector< ensembl_id > list;
	typedef boost::shared_ptr< list > list_ptr;

	std::string _prefix;
	unsigned _num;

	ensembl_id( const std::string & prefix = "", unsigned num = 0 );

	bool operator<( const ensembl_id & rhs ) const;
	bool operator==( const ensembl_id & rhs ) const;

	std::string str() const;

	/** Does the text look parseable into an ensembl id?. */
	static bool looks_like( const std::string & text );
	static ensembl_id parse( const std::string & text );
};
std::ostream &
operator<<( std::ostream &, const ensembl_id & );



/**
Defines an exon.
*/
struct exon
	: boost::equality_comparable< exon >
	, boost::less_than_comparable< exon >
{
	typedef std::vector< exon > list;
	typedef boost::shared_ptr< list > list_ptr;

	location _location;
	ensembl_id _id;

	exon( );
	exon( location l, const ensembl_id & id );

	bool operator==( const exon & rhs ) const;
	bool operator<( const exon & rhs ) const;
};


/**
Identifies a particular ensembl database.
*/
struct ensembl_database_id
	: boost::equality_comparable< ensembl_database_id >
	, boost::less_than_comparable< ensembl_database_id >
{
	species _species;
	int _software_version;
	int _ncbi_build;
	std::string _build_version;		/**< Normally a lowercase character, 'a', 'b', 'c', etc... */

	ensembl_database_id( );

	ensembl_database_id(
		const species & s,
		int software_version,
		int ncbi_build,
		const std::string & build_version );

	ensembl_database_id( const std::string & id );

	bool operator<( const ensembl_database_id & rhs ) const;
	bool operator==( const ensembl_database_id & rhs ) const;

	std::string str() const;
	static bool looks_like( const std::string & );
	static ensembl_database_id parse( const std::string & );
};
std::ostream &
operator<<( std::ostream &, const ensembl_database_id & );


/**
Identifies a sequence by gene id, transcript id, ensembl database and version.
*/
struct alignment_sequence_id
	: boost::equality_comparable< alignment_sequence_id >
	, boost::less_than_comparable< alignment_sequence_id >
	, boost::noncopyable
{
	typedef boost::shared_ptr< alignment_sequence_id > ptr;
	typedef std::vector< ptr > list;
	typedef boost::shared_ptr< list > list_ptr;

	ensembl_id _gene_id;
	ensembl_id _transcript_id;
	ensembl_database_id _db_id;
	int _version;					/**< A counter to differentiate different alignment runs. */

	alignment_sequence_id( );

	alignment_sequence_id( 
		const ensembl_id & gene_id,
		const ensembl_id & transcript_id,
		const ensembl_database_id & db_id,
		int version );

	alignment_sequence_id( const std::string & id );

	bool operator<( const alignment_sequence_id & rhs ) const;
	bool operator==( const alignment_sequence_id & rhs ) const;

	std::string str() const;
	static bool looks_like( const std::string & );
	static alignment_sequence_id::ptr parse( const std::string & );
};
std::ostream &
operator<<( std::ostream &, const alignment_sequence_id & );


/** 
Holds information about a sequence, such as position, the exons, ...
*/
struct alignment_sequence_info
	: boost::noncopyable
{
	typedef boost::shared_ptr< alignment_sequence_info > ptr;

	int _length;
	bool _has_position;
	int _position;
	exon::list_ptr _exons;
	region _region;

	alignment_sequence_info( );

	alignment_sequence_info(
		int length,
		bool has_position,
		int position,
		exon::list_ptr exons,
		region r );
};



/**
A sequence that makes up part of a remo in one species.
*/
struct remo_sequence
	: boost::equality_comparable< remo_sequence >
	, boost::less_than_comparable< remo_sequence >
	, boost::noncopyable
{
	typedef boost::shared_ptr< remo_sequence > ptr;
	typedef std::vector< ptr > list;
	typedef boost::shared_ptr< list > list_ptr;

	sequence _masked_sequence;
	sequence _unmasked_sequence;
	location _location;
	location _target_location;
	unsigned _conservation;
	unsigned _repeat_ratio;
	double _belief;

	remo_sequence( );

	remo_sequence(
		const sequence & masked_sequence,
		const sequence & unmasked_sequence,
		location l,
		location target_l,
		unsigned conservation,
		unsigned repeat_ratio,
		double belief );

	const sequence & get_sequence( bool masked = true ) const;

	bool operator==( const remo_sequence & rhs ) const;
	bool operator<( const remo_sequence & rhs ) const;
};


struct smart_ptr_less_than
{
	template< typename T >
	bool operator()( T p1, T p2 ) const
	{
		return p1.get() < p2.get();
	}
};

struct smart_ptr_less_than_value
{
	template< typename T >
	bool operator()( T p1, T p2 ) const
	{
		return *( p1.get() ) < *( p2.get() );
	}
};

/**
A regulatory module.
*/
struct module
	: boost::equality_comparable< module >
	, boost::less_than_comparable< module >
	, boost::noncopyable
{
	typedef boost::shared_ptr< module > ptr;
	typedef std::vector< ptr > list;
	typedef boost::shared_ptr< list > list_ptr;
	typedef std::map< alignment_sequence_id::ptr, remo_sequence::list_ptr, smart_ptr_less_than_value > sequence_remo_map;

	sequence_remo_map _sequences;

	alignment_sequence_id::list_ptr get_sequence_ids() const;
	remo_sequence::list_ptr get_sequences( alignment_sequence_id::ptr ) const;
	sequence get_sequence_for( alignment_sequence_id::ptr, bool masked = true ) const;

	bool operator==( const module & rhs ) const;
	bool operator<( const module & rhs ) const;
};



/**
Defines the sequences that have been aligned to search for remos. 
*/
struct aligned_sequence_set
	: boost::equality_comparable< aligned_sequence_set >
	, boost::less_than_comparable< aligned_sequence_set >
	, boost::noncopyable
{
	typedef boost::shared_ptr< aligned_sequence_set > ptr;
	typedef std::vector< ptr > list;
	typedef boost::shared_ptr< list > list_ptr;
	typedef std::map< alignment_sequence_id::ptr, alignment_sequence_info::ptr, smart_ptr_less_than_value > alignment_sequence_info_map;

	alignment_sequence_info_map _sequences;			/**< The sequences that were aligned with info about them. */
	alignment_sequence_id::ptr _centre_sequence;	/**< The sequence used in the centre of the alignment. */			

	alignment_sequence_id::list_ptr get_sequence_ids() const;
	alignment_sequence_info::ptr get_sequence_info( alignment_sequence_id::ptr seq_id ) const;

	bool operator==( const aligned_sequence_set & rhs ) const;
	bool operator<( const aligned_sequence_set & rhs ) const;
};



/**
Maps ensembl ids (genes probably) to aligned sequences.
*/
struct ensembl_id_alignment_map
	: boost::noncopyable
{
	typedef boost::shared_ptr< ensembl_id_alignment_map > ptr;
	typedef std::map< ensembl_id, aligned_sequence_set::list_ptr > ensembl_aligned_sequences_map;

	ensembl_aligned_sequences_map _map;

	ensembl_id::list_ptr get_genes( const std::string & prefix = "" ) const;
	aligned_sequence_set::list_ptr get_alignments_for( const ensembl_id & id ) const;
};

/**
Takes a list of aligned sequences and maps gene ids to them.
*/
ensembl_id_alignment_map::ptr
make_gene_alignment_map( 
	aligned_sequence_set::list_ptr aligned_sequence_set
);



/**
A space of ReMos. Maps aligned sequences to lists of remos.
*/
struct remome
	: boost::noncopyable
{
	typedef boost::shared_ptr< remome > ptr;
	typedef std::map< aligned_sequence_set::ptr, module::list_ptr, smart_ptr_less_than_value > remo_map;
	typedef boost::shared_ptr< remo_map > remo_map_ptr;

	remo_map_ptr _remos;

	aligned_sequence_set::list_ptr get_aligned_sequences() const;
	module::list_ptr get_remos_for( aligned_sequence_set::ptr ) const;
	void remove_remos_for( aligned_sequence_set::ptr );

	static remome::ptr deserialise( const std::string & path );
	void serialise( const std::string & path ) const;
};

remome::ptr load_remome_from_file( const std::string & filename );

remome::ptr parse_remome_from_file( const std::string & filename );

std::string get_remo_id( aligned_sequence_set::ptr aligned_seqs, module::ptr _remo );

typedef boost::tuple< module::ptr, aligned_sequence_set::ptr > remo_locator;
remo_locator get_remo_from_id( remome::ptr _remome, const std::string & id );

sequence_vec_ptr get_sequences_for_remo( aligned_sequence_set::ptr aligned_seqs, module::ptr _remo );



} //namespace remo
} //namespace biopsy

#endif //BIOPSY_REMO_H_
