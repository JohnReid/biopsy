/**
@file

Copyright John Reid 2006

*/

#include "biopsy/gapped_pssm_hmm.h"
#include "biopsy/python.h"
#include "biopsy/python_convert.h"

using namespace boost;
using namespace boost::python;
using namespace std;


namespace biopsy {

numpy_converter converter;

namespace gapped_pssm {

	
object sequence_to_python( const sequence & seq )
{
	boost::scoped_array< char > s( new char[ seq.size() + 1 ] );

	unsigned i = 0;
	BOOST_FOREACH( dna_base b, seq )
	{
		switch( b )
		{
		case 0: s[ i ] = 'a'; break;
		case 1: s[ i ] = 'c'; break;
		case 2: s[ i ] = 'g'; break;
		case 3: s[ i ] = 't'; break;
		case 4: s[ i ] = 'n'; break;
		default: throw std::logic_error( "Bad base in sequence" );
		}
		++i;
	}
	s[ i ] = '\0';
	return python::str( ( const char * ) s.get() );
}

void python_to_sequence( object python_seq, sequence & seq )
{
	const std::string str = extract< const char * >( python_seq.attr( "__str__" )() );
	seq.clear();
	BOOST_FOREACH( char c, str )
	{
		switch( c )
		{
		case 'a': case 'A': seq.push_back( 0 ); break;
		case 'c': case 'C': seq.push_back( 1 ); break;
		case 'g': case 'G': seq.push_back( 2 ); break;
		case 't': case 'T': seq.push_back( 3 ); break;
		case 'n': case 'N': seq.push_back( 4 ); break;
		default: throw std::logic_error( "Bad base in sequence" );
		}
	}
}

object sequences_to_python( const sequence_vec & seqs )
{
	python::list result;
	BOOST_FOREACH( const sequence & seq, seqs )
	{
		result.append( sequence_to_python( seq ) );
	}
	return result;
}

void python_to_sequences( object python_seqs, sequence_vec & seqs )
{
	seqs.clear();
	for ( unsigned i = 0; int( i ) < len( python_seqs ); ++i )
	{
		seqs.push_back( sequence() );
		python_to_sequence( python_seqs[i], seqs[ i ] );
	}
}


namespace hmm {

object p_r_given_predecessor( const model & m )
{
	double_array p_r;
	calculate_p_r_given_predecessor( m.trans_params, *( m.var_dist ), p_r );
	return convert_to_python( p_r );
}

object p_predecessor_given_r( const model & m )
{
	double_array p_r;
	calculate_p_r_given_predecessor( m.trans_params, *( m.var_dist ), p_r );
	double_array p_pre;
	calculate_p_predecessor_given_r( p_r, p_pre );
	return convert_to_python( p_pre );
}

object r_mode( const variational_distribution & var_dist )
{
	unsigned_vector_vec r;
	var_dist.r_mode( r );
	return convert_to_python( r );
}

object draw_sequence( const hidden_data & data, unsigned I )
{
	sequence seq;
	unsigned_vector r;
	data.draw_sequence( I, r, seq );
	return boost::python::make_tuple( sequence_to_python( seq ), convert_to_python( r ) );
}


observed_sequences::ptr 
create_sequences( object python_seqs )
{
	observed_sequences::ptr result( new observed_sequences );
	python_to_sequences( python_seqs, result->_sequences );
	return result;
}

observed_data::ptr
create_data(
	unsigned K,
	observed_sequences::ptr seqs,
	object psi,
	object theta,
	object phi,
	object upsilon )
{
	double_vector Psi;
	double_vector Theta;
	double_vector Phi;
	double_vector Upsilon;

	converter.from_numpy( psi, Psi );
	converter.from_numpy( theta, Theta );
	converter.from_numpy( phi, Phi );
	converter.from_numpy( upsilon, Upsilon );

	observed_data::ptr result(
		new observed_data(
			K,
			*seqs,
			Psi,
			Theta,
			Phi,
			Upsilon ) );
	return result;
}

object get_sequences( const observed_sequences & sequences )
{
	return sequences_to_python( sequences._sequences );
}

void test_state_map( const state_map_uncached & map, bool verbose = false )
{
	for( unsigned s = 0; map.num_states() != s; ++s )
	{
		const bool g = map.g( s );
		const bool c = map.c( s );
		const unsigned k = map.k( s );
		const unsigned calc_s = map.s( g, c, k );
		if( verbose ) std::cout << s << "," << g << "," << c << "," << k << "," << calc_s << "\n";
		if( calc_s != s ) throw std::logic_error( "Problem with s" );
	}

	for( unsigned k = 0; map.K + 1 != k; ++k ) //for each base
	{
		for( unsigned g = 0; 2 != g; ++g ) //whether gap or not
		{
			if( 0 == k && 1 == g ) continue; //no gap for background base
			if( map.K == k && 1 == g ) continue; //no gap for last base

			for( unsigned c = 0; 2 != c; ++c ) //whether rev-comp or not
			{
				if( 0 == k && 1 == c ) continue; //no rev-comp for background

				const unsigned s = map.s( bool( g ), bool( c ), k );
				if( bool( g ) != map.g( s ) ) throw std::logic_error( "Bad g in test_state_map" );
				if( bool( c ) != map.c( s ) ) throw std::logic_error( "Bad c in test_state_map" );
				if( k != map.k( s ) ) throw std::logic_error( "Bad k in test_state_map" );
			}
		}
	}
}

object predecessor_states( const state_map_uncached & map )
{
	unsigned_vector_vec predecessor_states;
	calculate_predecessor_states( map.K, predecessor_states );
	return convert_to_python( predecessor_states );
}

void export_gapped_pssms_hmm()
{
	using boost::python::arg;

	class_< 
		observed_sequences
	>( 
		"ObservedSequences",
		"Observed sequences",
		no_init
	)
		.def(
			"__init__",
			make_constructor( create_sequences ) 
		)
		.add_property(
			"sequences",
			get_sequences,
			"The sequences" 
		)
		.add_property(
			"N",
			&observed_sequences::N,
			"# sequences" 
		)
		.def(
			"I",
			&observed_sequences::I,
			"# bases in sequence n" 
		)
		.def(
			"x",
			&observed_sequences::x,
			"Base at position i in sequence n" 
		)
		;
	register_ptr_to_python< observed_sequences::ptr >();

	class_< 
		state_map_uncached
	>( 
		"StateMap",
		"Maps from state indices to their attributes and vice versa",
		init< unsigned >()
	)
	.def_readonly( "K", &state_map_uncached::K, "pssm size (excluding gap bases)" )
	.add_property( "S", &state_map_uncached::num_states, "# states" )
	.def( "g", &state_map_uncached::g, "Does the state have a gap?" )
	.def( "c", &state_map_uncached::c, "Is the state reverse complemented?" )
	.def( "b", &state_map_uncached::b, "Is it the background state?" )
	.def( "k", &state_map_uncached::k, "The index into the pssm" )
	.def( "m", &state_map_uncached::m, "The index into the emission probabilities" )
	.def( "s", &state_map_uncached::s, "The state index given the state attributes" )
	.def( "predecessor_states", predecessor_states, "Maps states to their predecessors" );
	;

	class_< 
		observed_data
	>( 
		"ObservedData",
		"Observed data",
		no_init
	)
		.def(
			"__init__",
			make_constructor( create_data ) 
		)
		.add_property(
			"X",
			&observed_data::X,
			"The sequences"
		)
		.def_readonly(
			"K",
			&observed_data::K,
			"# pssm bases (excluding gaps)" 
		)
		.add_property(
			"S",
			&observed_data::S,
			"# states"
		)
		.add_property(
			"E",
			&observed_data::E,
			"# different emission distributions"
		)
		.add_property(
			"Psi",
			access_and_convert< observed_data, double_vector, &observed_data::Psi >,
			"Prior for pssm emission distribution"
		)
		.add_property(
			"Theta",
			access_and_convert< observed_data, double_vector, &observed_data::Theta >,
			"Prior for background emission distribution"
		)
		.add_property(
			"Phi",
			access_and_convert< observed_data, boost::array< double, 2 >, &observed_data::Phi >,
			"Prior for pssm transition distribution"
		)
		.add_property(
			"Upsilon",
			access_and_convert< observed_data, boost::array< double, 2 >, &observed_data::Upsilon >,
			"Prior for background transition distribution"
		)
		;
	register_ptr_to_python< observed_data::ptr >();

	class_< 
		variational_distribution
	>( 
		"VariationalDistribution",
		"Variational distribution over hidden variables",
		init< observed_data >()
	)
	.add_property(
		"rho",
		access_and_convert< variational_distribution, double_vector_vec_vec, &variational_distribution::rho >,
		convert_and_set< variational_distribution, double_vector_vec_vec, &variational_distribution::rho >,
		"variational dist over hidden variable r : the state for each base"
	)
	.add_property(
		"eta",
		access_and_convert< variational_distribution, double_array, &variational_distribution::eta >,
		convert_and_set< variational_distribution, double_array, &variational_distribution::eta >,
		"variational dist over hidden variable e : the emission probabilities"
	)
	.add_property(
		"tau",
		access_and_convert< variational_distribution, double_array, &variational_distribution::tau >,
		convert_and_set< variational_distribution, double_array, &variational_distribution::tau >,
		"variational dist over hidden variable t : the transition probabilities"
	)
	.add_property(
		"r_mode",
		r_mode,
		"The mode of r under the variational distribution"
	)
	;
	register_ptr_to_python< variational_distribution::ptr >();

	class_< 
		hidden_data
	>( 
		"HiddenData",
		"The model's hidden data",
		init< const observed_data & >()
	)
	.add_property(
		"e",
		access_and_convert< hidden_data, double_array, &hidden_data::e >,
		"Emission probabilities"
	)
	.add_property(
		"t",
		access_and_convert< hidden_data, double_vector, &hidden_data::t >,
		"Transition probabilities"
	)
	.def(
		"draw_sequence",
		draw_sequence,
		"Draw a sequence from the hidden data.",
		( 
			arg( "I" )
		)
	)
	;

	class_< 
		model
	>( 
		"Model",
		"Model to do variational inference in",
		init< observed_data::ptr >()
	)
	.add_property(
		"data",
		&model::data,
		"The observed data"
	)
	.add_property(
		"var_dist",
		&model::var_dist,
		"The variational distribution over the model's hidden variables"
	)
	.add_property(
		"predecessor_states",
		access_and_convert< model, unsigned_vector_vec, &model::predecessor_states >,
		"Those states that can preceed each indexed state"
	)
	.add_property(
		"p_r_given_predecessor",
		p_r_given_predecessor,
		"numpy.array containing p(state|predecessor)"
	)
	.add_property(
		"p_predecessor_given_r",
		p_predecessor_given_r,
		"numpy.array containing p(predecessor|state)"
	)
	.def(
		"update",
		&model::update,
		"Update the model.",
		( 
			arg( "update_rho" ) = true,
			arg( "update_eta" ) = true,
			arg( "update_tau" ) = true
		)
	)
	;

	def( 
		"test_state_map",
		test_state_map,
		"Tests transformations between states and bases",
		( arg( "map" ), arg( "verbose" ) = false ) );
}










} //namespace hmm
} //namespace gapped_pssm
} //namespace biopsy


BOOST_PYTHON_MODULE( _gapped_pssms_hmm )
{
	biopsy::gapped_pssm::hmm::export_gapped_pssms_hmm();
}

