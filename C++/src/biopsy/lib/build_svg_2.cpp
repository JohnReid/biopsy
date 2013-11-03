/**
@file

Copyright John Reid 2007-2010

*/

#include "biopsy/defs.h"
#include "biopsy/binding_hits.h"
#include "biopsy/analyse.h"
#include "biopsy/transfac.h"
#include "biopsy/build_svg.h"

#include "bio/open_file.h"
#include "bio/environment.h"
#include "bio/biobase_db.h"
#include "bio/svg.h"
#include "bio/svg_pssm_logo.h"
#include "bio/xml_builder.h"
#include "bio/biobase_filter.h"
#include "bio/pathway_associations.h"

#include <boost/thread/mutex.hpp>
#include <boost/thread/locks.hpp>
#include <boost/assign/list_inserter.hpp>
USING_BIO_NS;

using namespace std;

XERCES_CPP_NAMESPACE_USE
using namespace boost;
using namespace std;

namespace biopsy {

using XERCES_CPP_NAMESPACE::DOMDocument;

const string svg_ns = "http://www.w3.org/2000/svg";
const string xlink_ns = "http://www.w3.org/1999/xlink";
const string ev_ns = "http://www.w3.org/2001/xml-events";

static const std::vector< std::string > svg_colours //an array of colours to be used in SVG
    = assign::list_of
/*
        ("white")
        ("aquamarine")
        ("blueviolet")
        ("brown")
        ("cadetblue")
        ("chartreuse")
        ("chocolate")
        ("crimson")
        ("cyan")
        ("darkgreen")
        ("darkkhaki")
        ("gold")
        ("lightpink")
        ("lightsalmon")
        ("lightskyblue")
        ("olive")
        ("orchid")
        ("peru")
        ("powderblue")
        ("red")
        ("steelblue")
        ("yellow")
        ("yellowgreen")
*/
        ("white")

        ("#3366FF")
        ("#6633FF")
        ("#FF3366")
        ("#FF6633")
        ("#33FF66")
        ("#66FF33")

        ("#33CCFF")
        ("#CC33FF")
        ("#33FFCC")
        ("#CCFF33")
        ("#FFCC33")
        ("#FF33CC")

        ("#003DF5")
        ("#3D00F5")
        ("#00F53D")
        ("#3DF500")
        ("#F5003D")
        ("#F53D00")

        ("#002EB8")
        ("#2E00B8")
        ("#00B82E")
        ("#2EB800")
        ("#B8002E")
        ("#B82E00")

        ;


double p_binding_min( const binding_hit::vec & hits )
{
    double _min = 1.0;
    BOOST_FOREACH( const binding_hit & hit, hits ) _min = std::min( _min, hit._p_binding );
    return _min;
}

double p_binding_max( const binding_hit::vec & hits )
{
    double _max = 0.0;
    BOOST_FOREACH( const binding_hit & hit, hits ) _max = std::max( _max, hit._p_binding );
    return _max;
}

DOMElement *
create_star( DOMDocument * doc, DOMElement * parent, double poly_x, double poly_y, double poly_size )
{
    const double star_outer_min_x = poly_x-poly_size*1.5;
    const double star_inner_min_x = poly_x-poly_size/2.0;
    const double star_outer_max_x = poly_x+poly_size*1.5;
    const double star_inner_max_x = poly_x+poly_size/2.0;
    const double star_outer_min_y = poly_y-poly_size*1.5;
    const double star_inner_min_y = poly_y-poly_size/2.0;
    const double star_outer_max_y = poly_y+poly_size*1.5;
    const double star_inner_max_y = poly_y+poly_size/2.0;
    return XmlBuilder(doc, parent)
        << XmlStartElement("polygon")
            << XmlSetAttribute(
                "points",
                BIOPSY_MAKE_STRING(
                       poly_x           << "," << star_outer_max_y << " "
                    << star_inner_max_x << "," << star_inner_max_y << " "
                    << star_outer_max_x << "," << poly_y           << " "
                    << star_inner_max_x << "," << star_inner_min_y << " "
                    << poly_x           << "," << star_outer_min_y << " "
                    << star_inner_min_x << "," << star_inner_min_y << " "
                    << star_outer_min_x << "," << poly_y           << " "
                    << star_inner_min_x << "," << star_inner_max_y << " " ));
}

const string & get_pathway_colour_for(const TableLink & pathway_link)
{
    assert( svg_colours.size() > interesting_pathways.size() ); //make sure we have enough (plus one for not identified pathways)

    //find the index of the pathway in interesting_pathways
    size_t index = 0;

    TableLinkVec::const_iterator interesting_pathway = find(interesting_pathways.begin(), interesting_pathways.end(), pathway_link);
    if (interesting_pathways.end() != interesting_pathway)
        index = interesting_pathway - interesting_pathways.begin() + 1;

    return svg_colours[index];
}




const std::string & bifa_javascript()
{
    static std::string _bifa_javascript;
    if( ! _bifa_javascript.size() )
    {
        const std::string script_file = BioEnvironment::singleton().get_svg_script_file_ver_2();
        ifstream is( script_file.c_str() );
        if( ! is ) throw std::logic_error( BIOPSY_MAKE_STRING( "Could not open: " << script_file ) );
        std::string line;
        while( ! is.eof() )
        {
            getline( is, line );
            if( ! is.eof() && ! is ) throw std::logic_error( BIOPSY_MAKE_STRING( "Problem reading: " << script_file ) );
            _bifa_javascript += line;
            _bifa_javascript += "\n";
        }
    }
    return _bifa_javascript;
}




void set_names( BiFaDetails & details )
{
    BOOST_FOREACH( const binding_hit & hit, details.hits )
        details.name_for_binder[ hit._binder_name ] = get_pssm_name( hit._binder_name );
};

void set_urls( BiFaDetails & details )
{
    BOOST_FOREACH( const binding_hit & hit, details.hits )
        details.url_for_binder[ hit._binder_name ] = get_pssm_url( hit._binder_name );
};

void set_pathways( BiFaDetails & details )
{
    BOOST_FOREACH( const binding_hit & hit, details.hits )
    {
        details.pathway_name_for_binder[ hit._binder_name ] = "Pathway not known";
        details.pathway_colour[ hit._binder_name ] = svg_colours[0];
        try {
            TableLink link( hit._binder_name );
            TableLink pathway_link = PathwayAssociations::singleton().get_most_significant_pathway_for( link );
            details.pathway_colour[ hit._binder_name ] = get_pathway_colour_for( pathway_link );
            if( PATHWAY_DATA == pathway_link.table_id ) {
                details.pathway_name_for_binder[ hit._binder_name ] = pathway_link.get_name();
                details.url_for_pathway[ pathway_link.get_name() ] = pathway_link.get_url();
            }
        } catch( ... ) {
        }
    }
}

/** Take some hits and output factors associated with them ordered in some way by strength of hits. */
void set_factors( BiFaDetails & details, unsigned max_num_factors=30 )
{
    unsigned i = 0;
    factor_scores_map_t factor_scores;
    BOOST_FOREACH( const binding_hit & hit, details.hits )
    {
        BiobaseTablePssmEntry * entry = 0;
        try {
            entry = BiobaseDb::singleton().get_pssm_entry(TableLink(hit._binder_name));
        } catch ( ... ) {
        }
        if( 0 != entry )
        {
            const FactorLinkList & factors = entry->get_factors();

            BOOST_FOREACH( const FactorLinkPtr f, factors )
                details.factors_for_binder[ hit._binder_name ].insert( f->name );

            //we give each factor hit by this matrix a fraction of the hits total value (1.0)
            const BIO_NS::float_t score_per_factor =
                1.0 == hit._p_binding
                    ? std::numeric_limits<BIO_NS::float_t>::max()
                    : BIO_NS::float_t(1.0) / BIO_NS::float_t(factors.size()) / (BIO_NS::float_t(1.0) - BIO_NS::float_t(hit._p_binding));

            for (FactorLinkList::const_iterator f = factors.begin();
                factors.end() != f;
                ++f)
            {
                FactorInfo & factor_info = factor_scores[(*f)->link];
                factor_info.score += score_per_factor;
                factor_info.hits.insert(i);
            }
        }
        ++i;
    }

    //for each factor - in score order
    while (! factor_scores.empty() && 0 != max_num_factors--)
    {
        //find the factor with the highest score
        factor_scores_map_t::iterator best = factor_scores.end();
        for (factor_scores_map_t::iterator f = factor_scores.begin();
            factor_scores.end() != f;
            ++f)
        {
            if (factor_scores.end() == best || f->second.score > best->second.score)
            {
                best = f;
            }
        }
        assert(best != factor_scores.end());

        //output
        BiFaDetails::factor f;
        f.id = best->first.get_text();
        f.name = best->first.get_name();
        f.hit_indices = best->second.hits;
        details.factors.push_back(f);

        //remove from factors
        factor_scores.erase(best);
    }
}

void build_bifa_details( BiFaDetails & details )
{
    set_names( details );
    set_urls( details );
    set_pathways( details );
    set_factors( details );
}



BiFaSvgBuilder::BiFaSvgBuilder(
    const sequence & seq,
    const BiFaDetails & details,
    double min_p_binding,
    double max_p_binding)
    : details(details)
    , hits(details.hits)
    , seq(seq)
    , p_binding_to_odds( new details::bifa_odds_mapped_range( details::bifa_p_binding_to_odds_mapping( 0, 1. ) ) )
    , vertical_guide_separation(50)
    , logo_height(10)
    , doc_x_margin(10)
    , graph_border(10)
    , graph_canvas_height(100)
    , graph_canvas_width(seq.size())
    , graph_height(graph_canvas_height + graph_border)
    , graph_width(graph_canvas_width + graph_border)
    , info_box_height(20)
    , doc_height(graph_height+info_box_height)
    , doc_width(graph_width)
    , impl(0)
    , doc_type(0)
    , doc(0)
    , doc_root(0)
    , defs_el(0)
    , script_el(0)
    , sequence_el(0)
    , info_box_el(0)
    , notes_panel_el(0)
    , factors_text_el(0)
    , graph_el(0)
    , graph_canvas_el(0)
    , hits_el(0)
    , underlines_el(0)
    , maximal_chain_el(0)
    , right_y_axis_el(0)
    , x_axis_el(0)
    , show_factor_list(true)
    , show_labels(false)
    , show_underlines(true)
    , y_axis_on_left(false)
    , label_y_axis_limits(true)
    , stroke_width(.3)
    , _show_sequence(show_sequence_multiple)
    , max_factors_in_hit_info(4)
{
    using namespace boost::assign;
    using namespace std;
    odds_to_show = { 1000, 10000, 20000, 50000, 100000, 200000, 500000 };
    //odds_to_show = list_of< int >( 1000 )( 10000 )( 20000 )( 50000 )( 100000 )( 200000 )( 500000 );
    boost::algorithm::to_upper(this->seq);
}

BiFaSvgBuilder & BiFaSvgBuilder::build()
{
    initialise();
    create_doc();
    create_sequence_elements();
    add_hits();
    return *this;
}

BiFaSvgBuilder::~BiFaSvgBuilder()
{
    delete doc;
    delete doc_type;
}

void BiFaSvgBuilder::initialise()
{
    // get DOM implementation
    impl = DOMImplementation::getImplementation();
    if( ! impl ) throw std::logic_error("Could not get DOM implementation");
}

void BiFaSvgBuilder::create_doc()
{
    //create document type and document
    doc_type = impl->createDocumentType(
        XStr("svg"),
        XStr("-//W3C//DTD SVG 1.1//EN"),
        XStr("http://www.w3.org/Graphics/SVG/1.1/DTD/svg11.dtd") );
    doc = impl->createDocument(
        XStr("http://www.w3.org/2000/svg"),        // root element namespace URI.
        XStr("svg"),        // root element name
        0 );            // document type object (DTD).
        //doc_type );            // document type object (DTD).
    doc->setEncoding(XStr("UTF-8"));

    //set attributes for document root
    doc_root = doc->getDocumentElement();
    XmlBuilder(doc, doc_root)
        << attr("version", "1.1")
        << attr("baseProfile", "full")
        << attr("viewBox", BIOPSY_MAKE_STRING(-int(doc_x_margin) << " 0 " << doc_width+2*doc_x_margin << " " << doc_height))
        << attr("preserveAspectRatio", "xMidYMid")
        << attr("onload", "init(evt);")
        << attr("pointer-events", "none")
        ;
    XmlBuilder(doc, doc_root)
        << attr("xmlns:xlink", xlink_ns, "")
        << attr("xmlns:ev", ev_ns, "")
        ;

    add_script();

    //add logo definitions
    defs_el = add_logo_defs(doc);

    create_info_box();
    if( "" != details.notes ) add_notes_label();
    create_graph();
    if( "" != details.notes ) add_notes();
}

void BiFaSvgBuilder::add_script(bool ref)
{
    script_el = XmlBuilder(doc, doc_root)
        << XmlStartElement("script")
            << attr("type", "text/javascript")
            ;
    if( ref ) XmlBuilder(doc, script_el) << attr("xlink:href", "bifa_ver_2.js", xlink_ns); //refer to script
    else script_el->appendChild( doc->createCDATASection(XStr(bifa_javascript())) ); //include script in svg
}

void BiFaSvgBuilder::create_info_box()
{
    //create info box element
    info_box_el = XmlBuilder(doc, doc_root)
        << XmlStartElement("g")
            << XmlStartElement("text")
                << attr("font-family", "Arial")
                << attr("font-size", "3px")
                << attr("id", "info_box_text")
                << attr("x", 0)
                << attr("y", info_box_height-5)
            << XmlEndElement()
        ;
    if( show_factor_list ) add_factors();

}

void BiFaSvgBuilder::add_notes_label()
{
    //create text to hover over to display notes
    XmlBuilder(doc, doc_root)
        << XmlStartElement("g")
            << XmlStartElement("text")
                << attr("font-family", "Arial")
                << attr("font-size", "3px")
                << attr("pointer-events", "visiblePainted")
                << attr("x", graph_canvas_width)
                << attr("y", info_box_height-5)
                << attr("text-anchor", "end")
                << XmlSetTextContent("notes")
                << attr("onmouseover", "show_notes();")
                << attr("onmouseout", "hide_notes();")
            << XmlEndElement()
        ;
}

void BiFaSvgBuilder::add_notes()
{
    notes_panel_el = XmlBuilder( doc, doc_root )
        << XmlStartElement("g")
            << attr( "id", "notes_panel" )
            << attr( "display", "none" )
            << attr( "transform", BIOPSY_MAKE_STRING("translate("<<graph_border<<" "<<2*info_box_height<<")"))
            << XmlStartElement( "rect" )
                << attr( "id", "notes_background" )
                << attr( "stroke", "black" )
                << attr( "stroke-width", stroke_width )
                << attr( "fill", "white" )
                << attr( "rx", stroke_width )
            << XmlEndElement()
        ;

    DOMElement * notes_text_el = XmlBuilder( doc, notes_panel_el )
        << XmlStartElement( "text" )
            << attr("font-family", "Arial")
            << attr("font-size", "3px")
        ;

    //separate the notes into individual lines
    typedef tokenizer< char_separator< char > > tokenizer;
    char_separator< char > sep( "\n" );
    tokenizer tokens( details.notes, sep );
    std::vector< std::string > lines;
    for( tokenizer::iterator tok_iter = tokens.begin();
        tok_iter != tokens.end();
        ++tok_iter )
    {
        XmlBuilder( doc, notes_text_el )
            << XmlStartElement( "tspan" )
                << attr( "x", "0px" )
                << attr( "dy", "4px" )
                << XmlSetTextContent( *tok_iter )
            << XmlEndElement()
            ;
    }
}

void BiFaSvgBuilder::add_factors()
{
    factors_text_el = XmlBuilder(doc, info_box_el)
        << XmlStartElement("text")
            << attr("font-family", "Arial")
            << attr("font-size", "3px")
            << attr("x", 0)
            << attr("y", info_box_height-10)
            << attr("pointer-events", "visiblePainted")
            ;
    double estimated_length = 0.;
    BOOST_FOREACH( const BiFaDetails::factor & f, details.factors )
    {
        stringstream os;
        std::copy( f.hit_indices.begin(), f.hit_indices.end(), ostream_iterator< unsigned >( os, "," ) );
        os << "null";
        const string args = os.str();
        XmlBuilder(doc, factors_text_el)
            << XmlStartElement("tspan")
                << XmlSetTextContent(f.name)
                << attr("id", BIOPSY_MAKE_STRING("factor_"<<f.id))
                << attr("onmouseover", BIOPSY_MAKE_STRING("mouse_over_factor("<<args<<");"))
                << attr("onmouseout", BIOPSY_MAKE_STRING("mouse_out_factor("<<args<<");"))
            ;
        estimated_length += 1.9*(1 + f.name.size());
        if( estimated_length > seq.size() ) break;
    }
}

void BiFaSvgBuilder::create_graph()
{
    //create graph element
    graph_el =
        XmlBuilder(doc, doc_root)
            << XmlStartElement("g")
                << attr("transform", BIOPSY_MAKE_STRING("translate(0 "<<info_box_height<<")"))
            ;

    //crate graph canvas
    graph_canvas_el =
        XmlBuilder(doc, graph_el)
            << XmlStartElement("g")
                << attr("transform", BIOPSY_MAKE_STRING("translate("<<(y_axis_on_left ? graph_border : 0)<<")"))
                << XmlStartElement("rect")
                    << attr("width", graph_canvas_width)
                    << attr("height", graph_canvas_height)
                    << attr("stroke", "black")
                    << attr("fill", "lightgray")
                    << attr("opacity", .1)
                << XmlEndElement()
            ;

    //place a vertical line on the canvas every so many bases
    for(unsigned i = 1; vertical_guide_separation*i < seq.size(); ++i)
        XmlBuilder(doc, graph_canvas_el)
            << XmlStartElement("line")
                << attr("x1", vertical_guide_separation*i)
                << attr("x2", vertical_guide_separation*i)
                << attr("y1", 0)
                << attr("y2", graph_canvas_height)
                << attr("stroke-width", stroke_width)
                << attr("stroke", "lightgray")
                ;

    label_axes();

    add_info();
}

void BiFaSvgBuilder::add_info()
{
    DOMElement * info_text = XmlBuilder(doc, graph_canvas_el)
        << XmlStartElement("g")
            << attr("transform", BIOPSY_MAKE_STRING("translate("<<graph_canvas_width/2<<" 13)"))
            << XmlStartElement("text")
                << attr("text-anchor", "middle")
                << attr("font-family", "Arial")
                << attr("font-size", "12px")
                << attr("fill", "lightgray")
                ;

    //separate the notes into individual lines
    typedef tokenizer< char_separator< char > > tokenizer;
    char_separator< char > sep( "\n" );
    tokenizer tokens( details.info, sep );
    std::vector< std::string > lines;
    for( tokenizer::iterator tok_iter = tokens.begin();
        tok_iter != tokens.end();
        ++tok_iter )
    {
        XmlBuilder( doc, info_text )
            << XmlStartElement( "tspan" )
                << attr( "x", "0px" )
                << attr( "dy", "13px" )
                << XmlSetTextContent( *tok_iter )
            << XmlEndElement()
            ;
    }
}

void BiFaSvgBuilder::create_sequence_elements()
{
    //add a sequence definition
    sequence_el =
        XmlBuilder(doc, defs_el)
            << XmlStartElement("text")
                << attr("font-family", "Arial")
                << attr("id", "sequence")
                << attr("font-size", "1.8px")
                << attr("fill", "lightgray")
                << attr("pointer-events", "visiblePainted")
                ;
    for (size_t i = 0; i < seq.size(); ++i) //add each base
        XmlBuilder(doc, sequence_el)
            << XmlStartElement("tspan")
                << attr("x", i)
                << XmlSetTextContent(seq.substr(i, 1))
                ;

//    if( _show_sequence == show_sequence_once || _show_sequence == show_sequence_multiple ) {
//        //put one copy of the sequence at the top
//        XmlBuilder(doc, graph_canvas_el)
//            << XmlStartElement("use")
//                << attr("xlink:href", "#sequence", xlink_ns)
//                << attr("x", 0)
//                << attr("y", -.5)
//                ;
//    }

    if( _show_sequence == show_sequence_multiple ) {
        BOOST_FOREACH( int odds, odds_to_show ) {
            const double y = p_binding_to_odds->backward( odds );
            XmlBuilder(doc, graph_canvas_el)
                << XmlStartElement("use")
                    << attr("xlink:href", "#sequence", xlink_ns)
                    << attr("x", 0)
                    << attr("y", (1.-y)*100.-.5*stroke_width)
                    ;
        }
    }
}

void BiFaSvgBuilder::label_axes()
{
    BOOST_FOREACH( int odds, odds_to_show ) {
        const double p_binding = p_binding_to_odds->backward( odds );
        XmlBuilder(doc, graph_canvas_el)
            << XmlStartElement("line")
                << attr("x1", 0)
                << attr("y1", (1.-p_binding)*100.)
                << attr("x2", graph_canvas_width)
                << attr("y2", (1.-p_binding)*100.)
                << attr("stroke-width", stroke_width)
                << attr("stroke", "lightgray")
                ;
    }

    const double label_y_shift = 2.;
    const int y_axis_x_offset = y_axis_on_left ? 0 : seq.size()+1;
    //label the y axis
    right_y_axis_el =
        XmlBuilder(doc, graph_canvas_el)
            << XmlStartElement("text")
                << attr("font-family", "Arial")
                << attr("font-size", "6px")
                << attr("fill", "gray")
                //<< attr("baseline-shift", "-35%")
                << attr("text-anchor", y_axis_on_left ? "end" : "start")
                ;
    XmlBuilder(doc, right_y_axis_el)
        << XmlStartElement("tspan")
            << attr("x", y_axis_x_offset)
            << attr("y", label_y_shift)
            << XmlSetTextContent( "Perfect" )
            ;
    BOOST_FOREACH( int odds, odds_to_show ) {
        const double p_binding = p_binding_to_odds->backward( odds );
        XmlBuilder(doc, right_y_axis_el)
            << XmlStartElement("tspan")
                << attr("x", y_axis_x_offset)
                << attr("y", (1.-p_binding)*100+label_y_shift)
                << XmlSetTextContent( BIOPSY_MAKE_STRING( odds/1000. << "Kb" ) )
                ;
    }
//    if( label_y_axis_limits ) {
//        XmlBuilder(doc, right_y_axis_el)
//            << XmlStartElement("tspan")
//                << attr("x", y_axis_x_offset)
//                << attr("y", 100+label_y_shift)
//                << XmlSetTextContent( BIOPSY_MAKE_STRING(int(y_axis->get_mapped_range()._min)) )
//            << XmlEndElement()
//            << XmlStartElement("tspan")
//                << attr("x", y_axis_x_offset)
//                << attr("y", label_y_shift)
//                << XmlSetTextContent( BIOPSY_MAKE_STRING(int(y_axis->get_mapped_range()._max)) )
//                ;
//    }
    if( y_axis_on_left ) {
        XmlBuilder(doc, right_y_axis_el) //a bit of a hack to ensure correct alignment when y axis is placed on left of figure
            << XmlStartElement("tspan")
                << attr("visibility", "hidden")
                << attr("x", y_axis_x_offset)
                << attr("y", -10+label_y_shift)
                << XmlSetTextContent( "0" )
                ;
    }

    //label the x axis
    x_axis_el =
        XmlBuilder(doc, graph_canvas_el)
            << XmlStartElement("text")
                << attr("font-family", "Arial")
                << attr("font-size", "6px")
                << attr("fill", "gray")
                << attr("text-anchor", "middle")
                ;
    for(unsigned i = 0; 100*i < seq.size(); ++i)
        XmlBuilder(doc, x_axis_el)
            << XmlStartElement("tspan")
                << attr("x", 100*i)
                << attr("y", graph_canvas_height + graph_border - 1)
                << XmlSetTextContent(BIOPSY_MAKE_STRING(100*i))
                ;
}

double BiFaSvgBuilder::hit_y( double p_binding ) const
{
    return ( 1. - p_binding ) * graph_canvas_height;
}

void BiFaSvgBuilder::add_hit_as_logo(const binding_hit & hit, unsigned idx)
{
    const unsigned num_bases = hit._location._length;
    const unsigned hit_width = logo_height*num_bases/2;
    const unsigned hit_height = logo_height;
    DOMElement * hit_el = XmlBuilder(doc, graph_canvas_el)
        << XmlStartElement("g")
            << attr(
                "transform",
                BIOPSY_MAKE_STRING(
                    "translate("
                    <<hit._location._position+num_bases/2.
                    <<" "
                    <<hit_y(hit._p_binding)-logo_height
                    <<")"))
            << XmlStartElement("rect")
                << attr("fill", "white")
                << attr("x", 0)
                << attr("y", 0)
                << attr("width", hit_width)
                << attr("height", hit_height)
                << attr("rx", hit_height/10.)
                << attr("opacity", .8)
            << XmlEndElement()
            ;

    const pssm_info & pssm_info = get_pssm(hit._binder_name);
    sequence hit_seq;
    if( hit._location._positive_strand ) {
        hit_seq = seq.substr(hit._location._position, hit._location._length);
    } else {
        reverse_complement(seq.substr(hit._location._position, hit._location._length), back_inserter(hit_seq));
    }

    const Pssm pssm = make_transfac_pssm( pssm_info._dists );
    DOMElement * svg_logo_el = create_svg_pssm_logo(
        pssm,
        doc,
        hit_seq);
    XmlBuilder(doc, svg_logo_el)
        << attr("x", 0)
        << attr("y", 0)
        << attr("width", hit_width)
        << attr("height", hit_height)
        //<< attr("preserveAspectRatio", "xMinYMin")
        ;
    hit_el->appendChild(svg_logo_el);
}

/**
    <g id="hit_29" transform="translate(527 98.3875)">
      <g display="none" id="hit_29_underlines">
        <path d="M -6 -97.8875 L 6 -97.8875" opacity="0.5" stroke="orange" stroke-width="1"/>
        <path d="M -6 -77.8875 L 6 -77.8875" opacity="0.5" stroke="orange" stroke-width="1"/>
        <path d="M -6 -57.8875 L 6 -57.8875" opacity="0.5" stroke="orange" stroke-width="1"/>
        <path d="M -6 -37.8875 L 6 -37.8875" opacity="0.5" stroke="orange" stroke-width="1"/>
        <path d="M -6 -17.8875 L 6 -17.8875" opacity="0.5" stroke="orange" stroke-width="1"/>
        <path d="M -6 2.11254 L 6 2.11254" opacity="0.5" stroke="orange" stroke-width="1"/>
      </g>
      <circle cx="0.0000" cy="0.0000" display="none" fill="white" id="hit_29_bg" opacity="0.9" r="7.5000" stroke="black" stroke-width="0.1000"/>
      <polygon fill="white" points="0,2.25 0.75,0.75 2.25,0 0.75,-0.75 0,-2.25 -0.75,-0.75 -2.25,0 -0.75,0.75 " stroke="black" stroke-width="0.2000">
        <desc>Pathway not known</desc>
        <title>V$RP58_01</title>
      </polygon>
      <g display="none" id="hit_29_info" opacity="0.7" transform="translate(3 0) scale(1.5)">
        <rect fill="white" id="hit_29_info_back" rx="1"/>
        <text font-family="Arial" style="font-size:10px" transform="scale(.2)">
          <a xlink:href="http://www.biobase-international.com/cgi-bin/biobase/transfac/current/bin/getTFProf.cgi?M00532">
            <tspan font-weight="bold" text-decoration="underline">V$RP58_01</tspan>
          </a>
          <tspan dx="1.2em" fill="blue">Pathway not known</tspan>
          <tspan dy="1.2em" fill="blue" x="0">Consensus: NNAACATCTGGA</tspan>
          <tspan dy="1.2em" fill="blue" x="0">Sequence: TTCACACCTGGA</tspan>
          <tspan dy="1.2em" fill="blue" x="0">Position: [521,533]</tspan>
          <tspan dy="1.2em" fill="blue" x="0">Factors: RP58, </tspan>
        </text>
        <svg> <!-- motif goes here -->
        </svg>
      </g>
    </g>
*/

void BiFaSvgBuilder::add_hit_as_shape(const binding_hit & hit, unsigned idx)
{
    const unsigned start = hit._location._position;
    const unsigned num_bases = hit._location._length;
    const pssm_info & pssm_info = get_pssm(hit._binder_name);
    const sequence hit_seq = seq.substr(start, num_bases);
    sequence pssm_seq;
    if( hit._location._positive_strand ) {
        pssm_seq = hit_seq;
    } else {
        reverse_complement(hit_seq, back_inserter(pssm_seq));
    }

    const Pssm pssm = make_transfac_pssm( pssm_info._dists );
    DOMElement * svg_logo_el = create_svg_pssm_logo(pssm, doc, pssm_seq);
    XmlBuilder(doc, svg_logo_el)
        << attr("height", logo_height)
        << attr("width", logo_height*num_bases/2)
        ;

    const string name = details.name_for_binder[ hit._binder_name ];
    const string url = details.url_for_binder[ hit._binder_name ];
    const string id = BIOPSY_MAKE_STRING("hit_"<<idx);
    const string pathway = details.pathway_name_for_binder[ hit._binder_name ];
    const string pathway_url = details.url_for_pathway[ pathway ];
    const string pathway_colour = details.pathway_colour[ hit._binder_name ];
    stringstream factors_stream;
    factors_stream << "Factors:";
    unsigned max_factors = max_factors_in_hit_info+1;
    BOOST_FOREACH( const string & f, details.factors_for_binder[ hit._binder_name ] ) {
        if( ! --max_factors ) {
            factors_stream << " ...";
            break;
        }
        factors_stream << " " << f;
    }
    const string factors = factors_stream.str();
    const bool in_max_chain = details.maximal_chain.end() != std::find( details.maximal_chain.begin(), details.maximal_chain.end(), hit );
    const double line_sep = 5.;
    const string text_colour = "black";
    const string position_string =
        hit._location._positive_strand
            ? BIOPSY_MAKE_STRING("Position: ["<<start<<","<<start+num_bases<<"]")
            : BIOPSY_MAKE_STRING("Position: ["<<start+num_bases<<","<<start<<"]");
    if( show_underlines ) {
        DOMElement * underlines_el = XmlBuilder(doc, hits_el)
            << XmlStartElement("g")
                << attr("id", BIOPSY_MAKE_STRING(id<<"_underlines"))
                << attr("display", "none")
                << attr("opacity", .5)
                ;
        BOOST_FOREACH( int odds, odds_to_show ) {
            const double y = p_binding_to_odds->backward( odds );
            XmlBuilder(doc, underlines_el)
                << XmlStartElement("line")
                    << attr("x1", start)
                    << attr("x2", start+num_bases)
                    << attr("y1", (1.-y)*100.+.5+.5*stroke_width)
                    << attr("y2", (1.-y)*100.+.5+.5*stroke_width)
                    << attr("stroke", "orange")
                << XmlEndElement()
                ;
        }
    }
    DOMElement * hit_el = XmlBuilder(doc, hits_el)
        << XmlStartElement("g")
            << attr(
                "transform",
                BIOPSY_MAKE_STRING(
                    "translate("
                    <<start+num_bases/2.
                    <<" "
                    <<hit_y(hit._p_binding)
                    <<")"))
            << attr("id", id)
            << XmlStartElement("rect")
                << attr("fill", "gray")
                << attr("display", "none")
                << attr("id", BIOPSY_MAKE_STRING(id<<"_bg"))
                << attr("opacity", .7)
                << attr("x", -int(num_bases)/2.-1.)
                << attr("y", -2)
                << attr("rx", .5)
                << attr("width", num_bases+2)
                << attr("height", 4)
                << attr("stroke", "black")
                << attr("stroke-width", stroke_width)
                << attr("stroke-dasharray", BIOPSY_MAKE_STRING(stroke_width<<","<<stroke_width))
            << XmlEndElement()
            << XmlStartElement("rect")
                << attr("fill", pathway_colour)
                << attr("x", -int(num_bases)/2.)
                << attr("y", -1)
                << attr("width", num_bases)
                << attr("height", 2)
                << attr("rx", .5)
                << attr("stroke", "black")
                << attr("stroke-width", stroke_width)
                << attr("pointer-events", "visiblePainted")
                << attr("onclick", BIOPSY_MAKE_STRING("click_hit('"<<id<<"')"))
                << attr("onmouseover", BIOPSY_MAKE_STRING("info('"<<name<<"')"))
                << attr("onmouseout", "info('')")
            << XmlEndElement()
            ;
    if( show_labels ) {
        XmlBuilder(doc, hit_el)
            << XmlStartElement("text")
                << attr("id", BIOPSY_MAKE_STRING(id<<"_label"))
                << attr("display", "inline")
                << attr("text-anchor", "middle")
                << attr("x", 0)
                << attr("y", 5)
                << attr("font-family", "Arial")
                << attr("font-size", "3px")
                << attr("fill", "gray")
                << XmlSetTextContent(name)
            << XmlEndElement()
            ;
    }
    if( in_max_chain )
        XmlBuilder(doc, create_star(doc, maximal_chain_el, start+num_bases/2., hit_y(hit._p_binding), 6))
            << attr("fill", pathway_colour)
            << attr("opacity", .5)
            << attr("stroke", "black")
            << attr("stroke-width", stroke_width)
            ;
    DOMElement * info_el = XmlBuilder(doc, hit_el)
        << XmlStartElement("g")
            << attr("id", BIOPSY_MAKE_STRING(id<<"_info"))
            << attr("pointer-events", "visiblePainted")
            << attr("display", "none")
            //<< attr("opacity", .9)
            << attr("transform", BIOPSY_MAKE_STRING("translate("<<num_bases/2.+1<<" 0) scale(.5) translate(0 -17)"))
            << XmlStartElement("rect")
                << attr("id", BIOPSY_MAKE_STRING(id<<"_info_back"))
                << attr("fill", "white")
                << attr("x", 0)
                << attr("y", 0)
                << attr("height", 0)
                << attr("width", 0)
                << attr("stroke", "lightgray")
                << attr("stroke-width", stroke_width)
                << attr("rx", stroke_width)
            << XmlEndElement()
            ;
    info_el->appendChild(svg_logo_el);
    XmlBuilder(doc, append_linked_text(info_el, name, url))
        << attr("y", logo_height+5)
        << attr("font-size", "5px")
        << attr("fill", text_colour)
        ;
    XmlBuilder(doc, append_linked_text(info_el, pathway, pathway_url))
        << attr("font-size", "5px")
        << attr("y", logo_height+5+line_sep)
        << attr("fill", pathway_colour == "white" ? text_colour : pathway_colour)
        << XmlSetTextContent(pathway)
        ;
    XmlBuilder(doc, info_el)
        << XmlStartElement("text")
            << attr("font-family", "Arial")
            << attr("font-size", "5px")
            << attr("fill", text_colour)
            << attr("y", logo_height+2*line_sep)
            << XmlStartElement("tspan")
                << attr("x", 0)
                << attr("dy", line_sep)
                << XmlSetTextContent(BIOPSY_MAKE_STRING("Sequence: "<<hit_seq))
            << XmlEndElement()
            << XmlStartElement("tspan")
                << attr("x", 0)
                << attr("dy", line_sep)
                << XmlSetTextContent(position_string)
            << XmlEndElement()
            << XmlStartElement("tspan")
                << attr("x", 0)
                << attr("dy", line_sep)
                << XmlSetTextContent(factors)
            << XmlEndElement()
            ;
}

DOMElement * BiFaSvgBuilder::append_linked_text(DOMElement * parent, const std::string & text, const std::string & url)
{
    if( "" == url ) { //no url
        return XmlBuilder(doc, parent)
            << XmlStartElement("text")
                << attr("font-family", "Arial")
                << XmlSetTextContent(text)
                ;
    } else {
        return XmlBuilder(doc, parent)
            << XmlStartElement("a")
                << XmlSetAttribute("xlink:href", url)
                << XmlStartElement("text")
                    << attr("font-family", "Arial")
                    << attr("font-weight", "bold")
                    << attr("text-decoration", "underline")
                    << XmlSetTextContent(text)
                ;
    }
}

void BiFaSvgBuilder::add_hits()
{
    maximal_chain_el = XmlBuilder(doc, graph_canvas_el)
        << XmlStartElement("g")
        ;

    underlines_el = XmlBuilder(doc, graph_canvas_el)
        << XmlStartElement("g")
        ;

    hits_el = XmlBuilder(doc, graph_canvas_el)
        << XmlStartElement("g")
        ;

    //put the hits on the canvas
    for(unsigned idx = 0; hits.size() != idx; ++idx) add_hit_as_shape(hits[idx], idx);
}

void
build_svg(
    const boost::filesystem::path & file,
    const string & title,
    const seq_t & seq,
    double min_threshold,
    const binding_hit::vec & hits,
    size_t max_num_factors,
    bool show_labels,
    bool open,
    const binding_hit::vec * max_chain,
    const std::string & notes,
    double max_threshold )
{
    //    xerces is not thread safe.
    //    Do one svg file at a time
    static boost::mutex  _m;
    boost::lock_guard<boost::mutex> lock(_m);
    {

        XMLPlatformUtils::Initialize();

        //make sure we have the best hit
        max_threshold = std::max( p_binding_max( hits ), max_threshold );
        if( max_threshold < min_threshold ) throw std::logic_error("No hits above minimum threshold");

        binding_hit::vec empty_chain;
        BiFaDetails details( hits, max_chain ? *max_chain : empty_chain );
        build_bifa_details( details );
        details.notes = notes;
        details.info = title;

        BiFaSvgBuilder bifa_builder(seq, details, min_threshold, max_threshold);
        bifa_builder.show_labels = show_labels;
        //bifa_builder.show_factor_list = false;
        bifa_builder.build();
        dom_print( bifa_builder.doc, file.string().c_str() );

        if( open ) open_file(file);
    }

    //    Terminate the XML after BiFaSvgBuilder has been deleted.
    XMLPlatformUtils::Terminate();
}

} //namespace biopsy
