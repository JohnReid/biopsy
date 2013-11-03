/**
@file

Copyright John Reid 2006

*/

#include "biopsy/defs.h"
#include "biopsy/gapped_pssm.h"
#include "biopsy/gapped_pssm_2.h"
#include "biopsy/python.h"
#include "biopsy/python_convert.h"

using namespace boost;
using namespace boost::python;
using namespace std;


namespace biopsy {

numpy_converter converter;

/**
Tests how quick the compressed dna data structure is compared with vector< char or unsigned >.

Not very quick apparently!
*/
struct test_compressed_dna
{
	std::vector< unsigned char > char_dna;
	std::vector< unsigned > unsigned_dna;
	compressed_dna_seq compressed_dna;
	test_compressed_dna( unsigned length )
		: char_dna( length )
		, unsigned_dna( length )
		, compressed_dna( length )
	{
		for( unsigned i = 0; length != i; ++i )
		{
			char_dna[i] = (unsigned char)( i % 4 );
			unsigned_dna[i] = i % 4;
			compressed_dna.set( i, i % 4 );
		}
	}

	unsigned test_char() const { 
		unsigned count = 0;
		for (unsigned i = 0; char_dna.size() != i; ++i )
		{
			if( char_dna[i] == 1 ) ++count;
		}
		return count;
	}

	unsigned test_unsigned() const {
		unsigned count = 0;
		for (unsigned i = 0; char_dna.size() != i; ++i )
		{
			if( unsigned_dna[i] == 1 ) ++count;
		}
		return count;
	}

	unsigned test_compressed() const {
		unsigned count = 0;
		for (unsigned i = 0; char_dna.size() != i; ++i )
		{
			if( compressed_dna.get< unsigned >(i) == 1 ) ++count;
		}
		return count;
	}
};

object dna_seq_list_to_python( const dna_seq_list & seqs )
{
	python::list result;
	BOOST_FOREACH( const compressed_dna_seq & seq, seqs )
	{
		std::string s;
		compressed_dna_to_string( seq, s );
		result.append( python::str( s.c_str() ) );
	}
	return result;
}

void python_to_dna_seq_list( object python_seqs, dna_seq_list & dna_seqs )
{
	dna_seqs.clear();
	for ( unsigned i = 0; i < python_seqs.attr( "__len__" )(); ++i )
	{
		const std::string str = 
			extract< const char * >( 
				python_seqs[i].attr( "__str__" )() );
		compressed_dna_seq seq;
		string_to_compressed_dna( str, seq );
		dna_seqs.push_back( seq );
	}
}

object dna_vec_list_to_python( const dna_vec_list & seqs )
{
	python::list result;
	BOOST_FOREACH( const dna_vec & seq, seqs )
	{
		std::string s;
		dna_vec_to_string( seq, s );
		result.append( python::str( s.c_str() ) );
	}
	return result;
}

void python_to_dna_vec_list( object python_seqs, dna_vec_list & dna_seqs )
{
	dna_seqs.clear();
	for ( unsigned i = 0; i < python_seqs.attr( "__len__" )(); ++i )
	{
		const std::string str = 
			extract< const char * >( 
				python_seqs[i].attr( "__str__" )() );
		dna_vec seq;
		string_to_dna_vec( str, seq );
		dna_seqs.push_back( seq );
	}
}



namespace gapped_pssm {

object get_seqs( const variational_model & model )
{
	return dna_vec_list_to_python( model.seqs );
}

double alpha( const variational_model & model, unsigned i )
{
	if( model.alpha.size() <= i ) throw std::logic_error( "Index out of range" );
	return model.alpha[i];
}

double varphi( const variational_model & model, unsigned i )
{
	if( model.varphi.size() <= i ) throw std::logic_error( "Index out of range" );
	return model.varphi[i];
}

double phi( const variational_model & model, unsigned i )
{
	if( model.phi.size() <= i ) throw std::logic_error( "Index out of range" );
	return model.phi[i];
}

double lambda( const variational_model & model, unsigned i )
{
	if( model.lambda.size() <= i ) throw std::logic_error( "Index out of range" );
	return model.lambda[i];
}

double eta( const variational_model & model, unsigned i )
{
	if( model.eta.size() <= i ) throw std::logic_error( "Index out of range" );
	return model.eta[i];
}

double mu( const variational_model & model, unsigned i )
{
	if( model.mu.size() <= i ) throw std::logic_error( "Index out of range" );
	return model.mu[i];
}

double omega( const variational_model & model, unsigned r, unsigned x )
{
	if( model.K+1 <= r ) throw std::logic_error( "Index out of range" );
	if( 4 <= x ) throw std::logic_error( "Index out of range" );
	return model.omega[r][x];
}

double nu( const variational_model & model, unsigned n, unsigned i )
{
	if( model.seqs.size() <= n ) throw std::logic_error( "Index out of range" );
	if( model.nu[n].size() <= i ) throw std::logic_error( "Index out of range" );
	return model.nu[n][i];
}

unsigned num_seqs( const variational_model & model )
{
	return model.seqs.size();
}

unsigned seq_length( const variational_model & model, unsigned n )
{
	if( model.seqs.size() <= n ) throw std::logic_error( "Index out of range" );
	return model.seqs[n].size();
}

variational_model create( 
	unsigned K,
	object seqs )
{
	dna_vec_list dna_seqs;
	python_to_dna_vec_list( seqs, dna_seqs );
	return variational_model( 
		K,
		dna_seqs,
		std::vector< double >( 2, 1.0 ),
		std::vector< double >( 4, 10.0 ),
		std::vector< double >( 4, 0.1 ) );
}

variational_model clone( const variational_model & model )
{
	return model;
}

struct compressed_dna_wrapper
{
	compressed_dna_seq s;
	compressed_dna_wrapper( const std::string & str ) { string_to_compressed_dna( str, s ); }
	compressed_dna_wrapper( const compressed_dna_seq & s ) : s( s ) { }
	operator std::string() const {
		std::string str; 
		compressed_dna_to_string( s, str ); 
		return str; 
	} 
	const char * get( unsigned i ) const {
		switch( s.get< unsigned >( i ) )
		{
		case 0: return "a";
		case 1: return "c";
		case 2: return "g";
		case 3: return "t";
		default: throw std::logic_error( "Base out of range - bad sequence" );
		}
	}
	void set( unsigned i, const std::string & bases )
	{
		for( unsigned j = 0; bases.size() != j; ++j )
		{
			switch( bases[j] )
			{
			case 'a': case 'A': s.set( i, 0 ); break;
			case 'c': case 'C': s.set( i, 1 ); break;
			case 'g': case 'G': s.set( i, 2 ); break;
			case 't': case 'T': s.set( i, 3 ); break;
			default: throw std::logic_error( "Unknown base" );
			}
		}
	}
};

void export_gapped_pssms()
{
	using boost::python::arg;

	class_< 
		test_compressed_dna
	>( 
		"TestCompressedDnaSeq",
		"Tests speed of access to compressed dna sequences",
		init< unsigned >()
	)
		.def(
			"test_char",
			&test_compressed_dna::test_char,
			"Tests char storage" )
		.def(
			"test_unsigned",
			&test_compressed_dna::test_unsigned,
			"Tests unsigned storage" )
		.def(
			"test_compressed",
			&test_compressed_dna::test_compressed,
			"Tests compressed storage" )
		;

	class_< 
		compressed_dna_wrapper
	>( 
		"CompressedDnaSeq",
		"Compresses unambiguous dna to 2 bits per base",
		init< std::string >()
	)
		.def(
			"__repr__",
			&compressed_dna_wrapper::operator std::string,
			"String representation" )
		.def(
			"get",
			&compressed_dna_wrapper::get,
			"The base at nucleotide i" )
		.def(
			"set",
			&compressed_dna_wrapper::set,
			"Sets the bases starting at given nucleotide index to the bases in the string argument" )
		;

	class_< 
		variational_model
	>( 
		"VariationalModel_C",
		"C implementation of gapped pssm variation model",
		no_init
	)
		.def(
			"create",
			create,
			( arg("seqs") ) )
		.staticmethod( "create" )
		.def(
			"clone",
			clone )
		.def_readonly(
			"K",
			&variational_model::K, 
			"Pssm length" )
		.add_property(
			"N",
			&num_seqs, 
			"Number of sequences in the model" )
		.add_property(
			"sequences",
			get_seqs,
			"The sequences" )
		.def(
			"sequence_length",
			seq_length,
			( arg("n") ),
			"The length of sequence n" )
		.def(
			"alpha",
			alpha,
			( arg("i") ),
			"Beta prior parameters for gamma (p(has gap))" )
		.def(
			"varphi",
			varphi,
			( arg("i") ),
			"Dirichlet prior parameters for background distribution" )
		.def(
			"phi",
			phi,
			( arg("i") ),
			"Dirichlet prior parameters for pssm distribution" )
		.def(
			"lambda_",
			lambda,
			( arg("i") ),
			"Variational parameter for gamma" )
		.def(
			"eta",
			eta,
			( arg("i") ),
			"Variational parameter for location of the gap" )
		.def(
			"mu",
			mu,
			( arg("i") ),
			"Variational parameter for g: has_gap variable" )
		.def(
			"nu",
			nu,
			( arg("n"), arg("i") ),
			"Variational parameters for start positions of sites" )
		.def(
			"omega",
			omega,
			( arg("r"), arg("x") ), 
			"Variational parameters for background and pss distributions" )
		.def( 
			"update",
			&variational_model::update,
			"Perform one update on the model" )
		.def( 
			"log_likelihood",
			&variational_model::log_likelihood, 
			"Likelihood of the sequences given the model" )
		.def( 
			"shift",
			&variational_model::shift, 
			( arg("offset") ), 
			"Shift the pssm and starting offsets by the given offset" )
		.def( 
			"blank_sites",
			&variational_model::blank_sites, 
			( arg("likelihood") = 0.8 ), 
			"Blank those sites that explain at least likelihood of sites. Returns number of sites blanked. (Blanked = set to 'n's)" )
		.def( 
			"initialise_variational_params",
			&variational_model::initialise_variational_params, 
			"Initialise variational parameters" )
		;
}

} //namespace gapped_pssm















namespace gapped_pssm_2 {
namespace impl {

observed_data::ptr 
create_observed_data( 
	object sequences,
	unsigned K,
	double expected_num_sites_per_seq,
	double strength_of_belief_in_num_sites_per_seq,
	double strength_of_gap_dist,
	double strength_of_pssm_dist )
{
	dna_vec_list dna_seqs;
	python_to_dna_vec_list( sequences, dna_seqs );
	observed_data::ptr result( 
		new observed_data( 
			dna_seqs,
			K,
			expected_num_sites_per_seq,
			strength_of_belief_in_num_sites_per_seq,
			strength_of_gap_dist,
			strength_of_pssm_dist ) );
	return result;
}

template< 
	typename Class,
	typename Container,
	Container Class:: * Member, 
	typename Result 
>
Result
get_element_1d(
	const Class & obj,
	unsigned index_0 )
{
	if( ( obj.*Member ).size() <= index_0 ) throw std::logic_error( BIOPSY_MAKE_STRING( "Index out of range: " << index_0 ) );
	return ( obj.*Member )[ index_0 ];
}

template< 
	typename Class,
	typename Container,
	Container Class:: * Member, 
	typename Result 
>
Result
get_element_2d(
	const Class & obj,
	unsigned index_0,
	unsigned index_1 )
{
	if( ( obj.*Member ).size() <= index_0 ) throw std::logic_error( BIOPSY_MAKE_STRING( "Index 0 out of range: " << index_0 ) );
	if( ( obj.*Member )[ index_0 ].size() <= index_1 ) throw std::logic_error( BIOPSY_MAKE_STRING( "Index 1 out of range: " << index_1 ) );
	return ( obj.*Member )[ index_0 ][ index_1 ];
}

object
expected_w( const variational_distribution & var_dist )
{
	double_array w;
	var_dist.calc_expected_w( w );
	return biopsy::converter.to_numpy( w );
}

} //namespace impl

void export_gapped_pssms_2()
{
	using namespace impl;
	using boost::python::arg;

	def(
		"r",
		calc_r,
		"Converts an offset into a gapped pssm binding site to r" );

	def(
		"offset",
		calc_offset,
		"Converts r an index into a gapped pssm to an offset into a binding site for the pssm" );

	class_< 
		observed_data
	>( 
		"ObservedData",
		"Observed data for gapped pssm model",
		no_init
	)
		.def(
			"__init__",
			make_constructor( create_observed_data ) 
		)
		.add_property(
			"A",
			access_and_convert< observed_data, double_vector, &observed_data::A >,
			"Prior for a : does site start at given base" 
		)
		.add_property(
			"B",
			access_and_convert< observed_data, double_vector, &observed_data::B >,
			"Prior for t : does given site have gap" 
		)
		.add_property(
			"V",
			access_and_convert< observed_data, double_vector, &observed_data::V >,
			"Prior for part of w : the gap base distribution" 
		)
		.add_property(
			"W",
			access_and_convert< observed_data, double_vector, &observed_data::W >,
			"Prior for part of w : the pssm base distribution" 
		)
		.def(
			"X",
			get_element_2d< observed_data, dna_vec_list, &observed_data::X, unsigned char >,
			"The observed sequences" 
		)
		.add_property(
			"N",
			&observed_data::N,
			"# of sequences" 
		)
		.def_readonly(
			"K",
			&observed_data::K,
			"# of bases in pssm (excluding gap)" 
		)
		.def(
			"L",
			&observed_data::L,
			"# bases in sequence n"
		)
		;



	class_< 
		variational_distribution
	>( 
		"VariationalDistribution",
		"Distribution over hidden variables in the model",
		init< const observed_data & >()
	)
	.add_property(
		"alpha",
		access_and_convert< variational_distribution, double_vector_vec, &variational_distribution::alpha >,
		"Variational distribution over a : does site start at base i in sequence n?" 
	)
	.add_property(
		"beta",
		access_and_convert< variational_distribution, double_vector_vec, &variational_distribution::beta >,
		"Variational distribution over b : does site at base i in sequence n have gap?" 
	)
	.add_property(
		"gamma",
		access_and_convert< variational_distribution, double_vector_vec, &variational_distribution::gamma >,
		"Variational distribution over c : is site at base i in sequence n a reverse complement?"
	)
	.add_property(
		"eta",
		access_and_convert< variational_distribution, double_vector, &variational_distribution::eta >,
		"Variational distribution over w : the base distribution." 
	)
	.add_property(
		"tau",
		access_and_convert< variational_distribution, double_vector, &variational_distribution::tau >,
		"Variational distribution over t : prob. any given site has gap." 
	)
	.add_property(
		"omega",
		access_and_convert< variational_distribution, double_array, &variational_distribution::omega >,
		"Variational distribution over w : the base distribution." 
	)
	.add_property(
		"expected_w",
		expected_w,
		"The expected w under the variational distribution"
	)
	;

	class_< 
		variational_model
	>( 
		"VariationalModel2",
		"Variational model",
		init< const observed_data & >()
	)
		.add_property(
			"data",
			&variational_model::data,
			"The observed data" )
		.add_property(
			"var_dist",
			&variational_model::var_dist,
			"The variational distribution" )
		.def(
			"update",
			&variational_model::update,
			"Update the variational distribution parameters." )
	;
}

} //namespace gapped_pssm
} //namespace biopsy


BOOST_PYTHON_MODULE( _gapped_pssms )
{
	biopsy::gapped_pssm::export_gapped_pssms();
	biopsy::gapped_pssm_2::export_gapped_pssms_2();
}

