/**
@file

Copyright John Reid 2006-2010

*/

#ifndef BIOPSY_BUILD_SVG_H_
#define BIOPSY_BUILD_SVG_H_

#ifdef _MSC_VER
# pragma once
#endif //_MSC_VER

#include "biopsy/defs.h"

#include "bio/xml_builder.h"

namespace biopsy
{

extern const std::string svg_ns;
extern const std::string xlink_ns;
extern const std::string ev_ns;

double p_binding_min( const binding_hit::vec & hits );
double p_binding_max( const binding_hit::vec & hits );

struct BiFaDetails
{
	struct factor {
		typedef std::vector< factor > vec;

		std::string id;
		std::string name;
		std::set< size_t > hit_indices;
	};

	typedef std::vector< std::string > string_vec;
	typedef std::set< std::string > string_set;
	typedef std::map< std::string, std::string > string_map;
	typedef std::map< std::string, string_set > string_set_map;

	std::string info; //displayed behind the hits
	std::string notes; //tooltip
	const binding_hit::vec & hits;
	const binding_hit::vec & maximal_chain;
	factor::vec factors;
	string_map name_for_binder;
	string_map url_for_binder;
	string_set_map factors_for_binder;
	string_map pathway_name_for_binder;
	string_map url_for_pathway;
	string_map pathway_colour;

	BiFaDetails( const binding_hit::vec & hits, const binding_hit::vec & maximal_chain ) : hits( hits ), maximal_chain( maximal_chain ) { }
};

void build_bifa_details( BiFaDetails & details );

/// Namespace for implementation details.
namespace details {

/// An interval.
template< typename T >
struct interval {
	T _min;
	T _max;
	interval( T _min, T _max ) : _min( _min ), _max( _max ) {
		if( _max < _min ) throw std::logic_error( BIOPSY_MAKE_STRING( "Bad arguments for interval constructor: "<<_min<<" > "<<_max ) );
	}
	bool inside( T x ) const { return _min <= x && x <= _max; }
	T size() const { return _max - _min; }
	T normalise( T x ) const { return (x - _min) / size(); }
};



/// A bijection can map forwards and backwards to and from its codomain and domain.
template< typename Forward, typename Backward >
struct bijection {
	Forward _forward;
	Backward _backward;
};


/// Map from p(binding) to odds.
struct p_binding_to_odds
{
	double operator()( double p_binding ) const {
		return p_binding / (1.0 - p_binding) / pssm_parameters::singleton().binding_background_odds_prior;
	}
};

/// Map from odds to p(binding).
struct odds_to_p_binding
{
	double operator()( double odds ) const {
		odds *= pssm_parameters::singleton().binding_background_odds_prior;
		return odds / (1. + odds);
	}
};


/// A bijection type from p(binding) to odds.
typedef bijection< p_binding_to_odds, odds_to_p_binding > bifa_p_binding_to_odds_bijection;


/// A bijection from p(binding) to odds.
inline
bifa_p_binding_to_odds_bijection
p_binding_odds_mapping() {
	return bijection< p_binding_to_odds, odds_to_p_binding >();
}



/// Calculates where markers should go in an interval.
template< typename OutputIt, typename T >
void
calculate_marker_positions(
	const interval< T > & range,
	OutputIt output_it )
{
	const double size = range.size();
	const int scale = int( log10( size ) );
	double spacing = pow( 10., scale-1 );
	const int num_markers = int( size / spacing );
	if( num_markers > 40 ) spacing *= 10.;
	else if( num_markers > 20 ) spacing *= 5.;
	else if( num_markers > 10 ) spacing *= 2.;
	const double begin = ceil( range._min / spacing ) * spacing;
	const double end = floor( range._max / spacing ) * spacing;
	for( double m = begin; m <= range._max; m += spacing ) {
		*output_it = m;
		++output_it;
	}
}


/// Interface for a mapped range.
template< typename T >
struct mapped_range
{
	typedef boost::shared_ptr< mapped_range > ptr;

	virtual ~mapped_range() { }

	virtual T forward( T x ) const = 0;
	virtual T backward( T y ) const = 0;
	virtual interval< T > get_mapped_range() const = 0;
	virtual interval< T > get_range() const = 0;
};




/// A concrete mapped range that depends on the templated Bijection.
template< typename Bijection, typename T >
struct concrete_mapped_range : mapped_range< T >
{
	Bijection _bijection;
	interval< T > _range;
	interval< T > _mapped_range;

	concrete_mapped_range( Bijection _bijection, interval< T > _range )
		: _bijection( _bijection )
		, _range( _range )
		, _mapped_range( _bijection._forward( _range._min ), _bijection._forward( _range._max ) )
	{ }

	virtual ~concrete_mapped_range() { }

	virtual T forward( T x ) const {
		return _bijection._forward( x );
	}

	virtual T backward( T x ) const {
		return _bijection._backward( x );
	}

	virtual interval< T > get_mapped_range() const {
		return _mapped_range;
	}

	virtual interval< T > get_range() const {
		return _range;
	}
};


/// Make a mapped range from a bijection and a range.
template< typename Bijection, typename T >
concrete_mapped_range< Bijection, T >
make_mapped_range( Bijection _bijection, interval< T > _range )
{
	return concrete_mapped_range< Bijection, T >( _bijection, _range );
}


/// A mapped range type from p(binding) to odds
typedef concrete_mapped_range< bifa_p_binding_to_odds_bijection, double > bifa_odds_mapped_range;


/// Make a mapped range from p(binding) to odds with the given limits.
inline
bifa_odds_mapped_range
bifa_p_binding_to_odds_mapping( double _min_p_binding, double _max_p_binding )
{
	return make_mapped_range( p_binding_odds_mapping(), interval< double >( _min_p_binding, _max_p_binding ) );
}


} //namespace details





/**
sets pathway info based on pssm details
*/
void set_pathways( BiFaDetails & details );




class BiFaSvgBuilder
{
public:
	typedef XERCES_CPP_NAMESPACE::DOMElement element;

	enum show_sequence {
		show_sequence_none,
		show_sequence_once,
		show_sequence_multiple,
	};


	BiFaDetails details;
	const binding_hit::vec & hits;
	sequence seq;
	details::mapped_range< double >::ptr            p_binding_to_odds;       ///< Map p(binding) to odds ratio.
	std::vector< int >                              odds_to_show;            ///< Those odds we show on y-axis
	unsigned vertical_guide_separation;
	unsigned logo_height;
	unsigned doc_x_margin;
	unsigned graph_border;
	unsigned graph_canvas_height;
	unsigned graph_canvas_width;
	unsigned graph_height;
	unsigned graph_width;
	unsigned info_box_height;
	unsigned doc_height;
	unsigned doc_width;
	XERCES_CPP_NAMESPACE::DOMImplementation * impl;
	XERCES_CPP_NAMESPACE::DOMDocumentType * doc_type;
	XERCES_CPP_NAMESPACE::DOMDocument * doc;
	element * doc_root;
	element * defs_el;
	element * script_el;
	element * sequence_el;
	element * info_box_el;
	element * notes_panel_el;
	element * factors_text_el;
	element * graph_el;
	element * graph_canvas_el;
	element * hits_el;
	element * underlines_el;
	element * maximal_chain_el;
	element * right_y_axis_el;
	element * x_axis_el;
	bool show_factor_list;
	bool show_labels;
	bool show_underlines;
	bool y_axis_on_left;
	bool label_y_axis_limits;
	double stroke_width;
	show_sequence _show_sequence;
	unsigned max_factors_in_hit_info;

	BiFaSvgBuilder(
		const sequence & seq,
		const BiFaDetails & details,
		double min_p_binding,
		double max_p_binding
	);

	~BiFaSvgBuilder();

	BiFaSvgBuilder & build();
	void initialise();
	void create_doc();
	void add_script(bool ref = false);
	void create_info_box();
	void add_notes();
	void add_notes_label();
	void add_factors();
	void create_graph();
	void add_info();
	void create_sequence_elements();
	void label_axes();
	double hit_y(double p_binding) const;
	void add_hit_as_logo(const binding_hit & hit, unsigned idx);
	void add_hit_as_shape(const binding_hit & hit, unsigned idx);
	void add_hits();
	element * append_linked_text(element * parent, const std::string & text, const std::string & url);

	template< typename Value >
	BIO_NS::XmlSetAttribute attr( const std::string & name, const Value & value, const std::string & ns_uri = svg_ns )
	{
		return BIO_NS::XmlSetAttribute( name, value, ns_uri );
	}



};

} //namespace biopsy

#endif //BIOPSY_BUILD_SVG_H_
