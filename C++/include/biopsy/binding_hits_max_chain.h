/**
@file

Copyright John Reid 2006

*/

#ifndef BIOPSY_BINDING_HITS_MAX_CHAIN_H_
#define BIOPSY_BINDING_HITS_MAX_CHAIN_H_

#ifdef _MSC_VER
# pragma once
#endif //_MSC_VER


#include "biopsy/defs.h"
#include "biopsy/binding_hits.h"

#include <bio/max_chain.h>



namespace biopsy {






/** 
Describes some traits of binding hits for use in max chain algorithms. 
*/
struct binding_hit_traits
{
	typedef binding_hit hit;
	typedef const hit * data;
	typedef std::string character;
	typedef double weight;
	typedef int coord;

	static data get_data( const hit & h );
	static const character & get_char( const hit & h );
	static weight get_weight( const hit & h );
	static coord get_start( const hit & h );
	static coord get_end( const hit & h );
};




} //namespace biopsy

#endif //BIOPSY_BINDING_HITS_MAX_CHAIN_H_

