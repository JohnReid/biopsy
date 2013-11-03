#ifndef BIOPSY_CUSTOM_PSSM_H_
#define BIOPSY_CUSTOM_PSSM_H_

#ifdef _MSC_VER
# pragma once
#endif //_MSC_VER

/**
Code to parse custom PSSM files.

E.g.

NA  P$ABF1.01
WI  12
PO  01 02 03 04 05 06 07 08 09 10 11 12
CA  00 02 20 00 20 00 00 00 00 00 00 08
CC  01 00 00 20 00 20 00 00 00 00 19 04
CG  14 17 00 00 00 00 20 00 20 16 00 02
CT  05 01 00 00 00 00 00 20 00 03 00 03
IU  G  G  A  C  A  C  G  T  G  G  C  N 
UR  http://www.genomatix.de/

*/



#include "biopsy/defs.h"
#include "biopsy/pssm.h"

namespace biopsy {

/** Info about a custom PSSM. */
struct custom_pssm
{
	typedef boost::shared_ptr< custom_pssm > ptr;

	std::string name;
	std::string id;
	std::string url;
	nucleo_dist::vec counts;
};

typedef std::set< std::string > pssm_set;
typedef boost::shared_ptr< pssm_set > pssm_set_ptr;

/** The filename for the custom pssm with this id. */
std::string custom_pssm_filename( const std::string & pssm_id );

/** Parse the custom pssm file. */
custom_pssm::ptr parse_custom_pssm_file( const std::string & custom_pssm_filename );

/** Get the custom pssm by id. */
custom_pssm::ptr get_custom_pssm( const std::string & pssm_id );

/** The names of the installed custom pssm sets. */
std::vector< std::string > get_installed_custom_pssm_set_names();
string_vec_ptr get_custom_pssm_set_names();

/** Get the custom pssm set by name. */
pssm_set_ptr get_pssm_set( const std::string & name );
string_vec_ptr get_custom_pssms(const std::string & pssm_set_name);

/** Get all the custom pssms. */
pssm_set_ptr all_custom_pssms();
string_vec_ptr get_all_custom_pssms();

/** Add pssms from named set to the pssm set given. */
void add_pssms_from_pssm_set(
	pssm_set & set,
	const std::string & pssm_set_name,
	bool use_transfac_consensus_sequences,
	const std::string & matrix_species,
	const std::string & matrix_name_match );

/** Build a pssm set from a range of pssm set names (I.e. the union). */
template< typename PssmSetNameRange >
pssm_set_ptr
build_pssm_set_from_names(
	const PssmSetNameRange & names,
	bool use_transfac_consensus_sequences,
	const std::string & matrix_species,
	const std::string & matrix_name_match )
{
	using namespace boost;
	pssm_set_ptr result( new pssm_set );
	//use transfac as default if range is empty
	add_pssms_from_pssm_set( *result, "transfac", use_transfac_consensus_sequences, matrix_species, matrix_name_match );
	BOOST_FOREACH( const typename range_value< PssmSetNameRange >::type & name, names )
	{
		add_pssms_from_pssm_set( *result, name, use_transfac_consensus_sequences, matrix_species, matrix_name_match );
	}
	return result;
}




} //namespace biopsy

#endif //BIOPSY_CUSTOM_PSSM_H_

