/* Copyright John Reid 2007
*/

#include "bio-pch.h"


#include "bio/defs.h"


#include "bio/xml_builder.h"
#include "bio/svg_match.h"
#include "bio/svg_pssm_logo.h"
#include "bio/biobase_db.h"
#include "bio/biobase_match.h"
#include "bio/biobase_data_traits.h"
#include "bio/pathway_associations.h"

XERCES_CPP_NAMESPACE_USE

using namespace boost;

#include <sstream>
#include <string>
#include <vector>
#include <set>
#include <iomanip>
using namespace std;

//the area containing the graph and associated axes
#define GRAPH_HEIGHT 300
#define GRAPH_BORDER 40


#define DEFAULT_TEXT_STYLE "fill: #000000; font-weight: normal; "

BIO_NS_START

using XERCES_CPP_NAMESPACE::DOMDocument;

DOMElement *
create_text_element(DOMDocument * doc, const string & name, const string & value)
{
    DOMElement * text_element = doc->createElement(XStr(name));
    text_element->setTextContent(XStr(value));
    return text_element;
}

DOMElement * create_svg(
    DOMDocument * doc,
    size_t x,
    size_t y,
    size_t width,
    size_t height,
    size_t view_box_width,
    size_t view_box_height,
    const std::string stroke_colour,
    const std::string fill_colour,
    float_t opacity)
{
    stringstream view_box_stream;
    view_box_stream << "0 0 " << view_box_width << " " << view_box_height << ends;

    DOMElement * result = doc->createElement(XStr("svg"));
    XmlBuilder builder(doc, result);
    builder
        << XmlSetAttribute("x", x)
        << XmlSetAttribute("y", y)
        << XmlSetAttribute("width", width)
        << XmlSetAttribute("height", height)
        << XmlSetAttribute("viewBox", view_box_stream.str())
        << XmlSetAttribute("preserveAspectRatio", "xMidYMid")
        << XmlStartElement("rect")
            << XmlSetAttribute("height", view_box_height)
            << XmlSetAttribute("width", view_box_width)
            << XmlSetAttribute("stroke", stroke_colour)
            << XmlSetAttribute("fill", fill_colour)
            << XmlSetAttribute("opacity", opacity)
        << XmlEndElement()
        ;

    return result;

}


std::string
build_transform_string(
    const std::string & type,
    float_t x,
    float_t y)
{
    stringstream transform_value;
    transform_value
        << "translate("
        << x
        << ","
        << y
        << ") ";
    return transform_value.str();
}


//one colour for each pathway
struct SvgColours : std::vector<std::string>
{
    SvgColours()
    {
        push_back("aquamarine");
        push_back("blueviolet");
        push_back("brown");
        push_back("cadetblue");
        push_back("chartreuse");
        push_back("chocolate");
        push_back("crimson");
        push_back("cyan");
        push_back("darkgreen");
        push_back("darkkhaki");
        push_back("gold");
        push_back("lightpink");
        push_back("lightsalmon");
        push_back("lightskyblue");
        push_back("olive");
        push_back("orchid");
        push_back("peru");
        push_back("powderblue");
        push_back("red");
        push_back("steelblue");
        push_back("yellow");
        push_back("yellowgreen");
    }
};
static const SvgColours colours; //an array of colours to be used in SVG


SvgDomDocument::SvgDomDocument(
    size_t num_bases,
    float_t min_threshold,
    float_t max_threshold,
    const std::string & title,
    const seq_t & seq,
    bool show_titles)
    : doc(0)
    , doc_root(0)
    , graph_graphic(0)
    , graph_area_graphic(0)
    , points_graphic(0)
    , defs_el(0)
    , factor_el(0)
    , notes_el(0)
    , height(GRAPH_HEIGHT)
    , width(float_t( num_bases + 2 * GRAPH_BORDER ))
    , num_bases(num_bases)
    , min_threshold(min_threshold)
    , max_threshold(max_threshold)
    , seq(seq)
    , show_titles(show_titles)
{
    DOMImplementation* impl
        = DOMImplementation::getImplementation();

    if (impl != NULL)
    {
        //DOMDocumentType* doc_type =
            impl->createDocumentType(
                XStr("svg"),
                XStr("-//W3C//DTD SVG 1.1//EN"),
                XStr("http://www.w3.org/Graphics/SVG/1.1/DTD/svg11.dtd") );
        doc = impl->createDocument(
                    XStr("http://www.w3.org/2000/svg"),        // root element namespace URI.
                    XStr("svg"),        // root element name
//                    doc_type);            // document type object (DTD).
                    0);            // document type object (DTD).

        doc->setEncoding(XStr("UTF-8"));

        //set attributes for document root
        {
            doc_root = doc->getDocumentElement();
            stringstream view_box_stream;
            view_box_stream << "0 0 " << width << " " << height << ends;
            set_attribute(doc_root, "viewBox", view_box_stream.str());
            set_attribute(doc_root, "preserveAspectRatio", "xMidYMid");
            set_attribute(doc_root, "xmlns:xlink", "http://www.w3.org/1999/xlink");
        }

        //add the svg logo definitions
        add_logo_defs(doc);

        //add a defs element
        defs_el = doc->createElement(XStr("defs"));
        doc_root->appendChild(defs_el);

        //add a sequence definition
        {
            DOMElement * seq_def_el = doc->createElement(XStr("text"));
            set_attribute(seq_def_el, "font-size", "1.8px");

            //add each base
            for (size_t i = 0; i < num_bases; ++i)
            {
                DOMElement * base_el = doc->createElement(XStr("tspan"));
                set_attribute(base_el, "x", (width - 2 * GRAPH_BORDER) * i / num_bases);
                base_el->setTextContent(XStr(seq.substr(i, 1)));
                seq_def_el->appendChild(base_el);
            }

            //identify and add to the defs element
            set_attribute(seq_def_el, "id", "sequence");
            defs_el->appendChild(seq_def_el);
        }

        //the graphic we put graph elements in
        //    <svg x="0" y="0" width="1000" height="400" viewBox="0 0 100 60" preserveAspectRatio="none">
        graph_graphic = create_svg(
            doc,
            0,
            0,
            size_t(width),
            GRAPH_HEIGHT,
            size_t(width),
            GRAPH_HEIGHT,
            "none",
            "none",
            0.1f);
        doc_root->appendChild(graph_graphic);

        build_graph(title);

        //the graphic transformed to graph coordinates
        {
            graph_area_graphic = doc->createElement(XStr("g"));
            stringstream transform_stream;
            transform_stream
                << "translate(" << GRAPH_BORDER << " " << GRAPH_BORDER << ") "
                << "scale(1 "
                << float_t(GRAPH_HEIGHT - 2 * GRAPH_BORDER) / float_t(100)
                << ")"
                << ends;
            set_attribute(graph_area_graphic, "transform", transform_stream.str());

            //
            // create the graphic which we put points in
            //
            points_graphic = doc->createElement(XStr("g"));
            graph_area_graphic->appendChild(points_graphic);
            
            graph_graphic->appendChild(graph_area_graphic);
        }

        //build_pathways();
    }
}

void SvgDomDocument::build_graph(const std::string & title)
{
    const size_t graph_area_height = GRAPH_HEIGHT - 2 * GRAPH_BORDER;
    const size_t graph_area_width = size_t(width) - 2 * GRAPH_BORDER;

    //
    // put a background in the graph
    //
    DOMElement * background_rect = doc->createElement(XStr("rect"));
    set_attribute(background_rect, "y", GRAPH_BORDER);
    set_attribute(background_rect, "x", GRAPH_BORDER);
    set_attribute(background_rect, "height", float_t( graph_area_height ));
    set_attribute(background_rect, "width", float_t( graph_area_width ));
    set_attribute(background_rect, "stroke", "none");
    set_attribute(background_rect, "fill", "gray");
    set_attribute(background_rect, "opacity", "0.2");
    graph_graphic->appendChild(background_rect);

    //
    // put text in the graph
    //
    {
        typedef tokenizer<char_separator<char> > tokenizer;
        char_separator<char> sep("\n");
        tokenizer tokens(title, sep);
        size_t spacing = 0;
        for (tokenizer::iterator tok_iter = tokens.begin();
            tok_iter != tokens.end();
            ++tok_iter, spacing += 30)
        {
            DOMElement * graph_text = doc->createElement(XStr("text"));
            set_attribute(graph_text, coord_t(float_t( width / 2 ), float_t( GRAPH_HEIGHT / 2 + spacing )));
            set_attribute(graph_text, "style", "text-anchor:middle; fill: #000000; font-size: 24px; font-weight: normal");
            graph_text->setTextContent(XStr(*tok_iter));
            set_attribute(graph_text, "opacity", "0.5");
            graph_graphic->appendChild(graph_text);
        }
    }

    //
    //axes text elements
    //
    DOMElement * axes_text = doc->createElement(XStr("text"));
    set_attribute(axes_text, "style", DEFAULT_TEXT_STYLE "font-size: 12px");
    set_attribute(axes_text, "opacity", "0.5");

    { //y axis
        for (size_t i = 0; i != 6; ++i)
        {
            const float_t t = min_threshold + get_range() * float_t(i) / float_t(5);
            const float_t y = GRAPH_BORDER + graph_area_height * (1.0f - (t - min_threshold) / get_range());

            axes_text->appendChild(
                create_label( //label for start of x axis
                    coord_t(GRAPH_BORDER, y),
                    num2str(t),
                    "text-anchor:end"));
            axes_text->appendChild(
                create_label( //label for end of x axis
                    coord_t(float_t( GRAPH_BORDER + graph_area_width ), y),
                    num2str(t),
                    "text-anchor:start"));

            //line
            graph_graphic->appendChild(
                create_line(
                    coord_t(float_t( GRAPH_BORDER ), y),
                    coord_t(float_t( GRAPH_BORDER + graph_area_width ), y)));

            DOMElement * seq_el = doc->createElement(XStr("use"));
            set_attribute(seq_el, "xlink:href", "#sequence");
            set_attribute(seq_el, coord_t(GRAPH_BORDER, y));
            graph_graphic->appendChild(seq_el);
        }
    }

    { //x axis
        for (size_t i = 0; i < num_bases; i += 100)
        {
            const float_t x = GRAPH_BORDER + float_t( graph_area_width * i ) / num_bases;

            //label every 100'th bases
            axes_text->appendChild(
                create_label(
                    coord_t(x, GRAPH_BORDER + graph_area_height),
                    num2str(i),
                    "text-anchor:middle; baseline-shift:-16px")); 

            //line every 100'th base
            graph_graphic->appendChild(
                create_line(
                    coord_t(x, GRAPH_BORDER),
                    coord_t(x, GRAPH_BORDER + graph_area_height)));
        }

        //label for end of x axis
        axes_text->appendChild(
            create_label(
                coord_t(GRAPH_BORDER + float_t( graph_area_width ), GRAPH_BORDER + float_t( graph_area_height )),
                num2str(num_bases),
                "text-anchor:middle; baseline-shift:-16px")); 

        //line at end of x axis
        graph_graphic->appendChild(
            create_line(
                coord_t(GRAPH_BORDER + float_t( graph_area_width ), float_t( GRAPH_BORDER )),
                coord_t(GRAPH_BORDER + float_t( graph_area_width ), GRAPH_BORDER + float_t( graph_area_height ))));

        //x axis label
        axes_text->appendChild(
            create_label(
                coord_t(GRAPH_BORDER + float_t( graph_area_width ) / 2, GRAPH_BORDER + float_t( graph_area_height )),
                "Position in sequence",
                "font-size: 16px; text-anchor:middle; baseline-shift:-30px"));
    }
    graph_graphic->appendChild(axes_text);

    // factor element
    {
        XmlBuilder builder(doc, graph_graphic);
        builder << XmlStartElement("text");
        factor_el = builder.current_element;
        builder
            << XmlSetAttribute("y", 10)
            << XmlSetAttribute("font-size", 10)
            << XmlSetTextContent("Factors: ");
        builder << XmlEndElement();
    }

    // notes element
    {
        XmlBuilder builder( doc, graph_graphic );
        builder << XmlStartElement( "text" );
        notes_el = builder.current_element;
        builder << XmlEndElement();
    }
}


void
SvgDomDocument::add_notes( const std::string & notes )
{
    // notes element
    XmlBuilder( doc, notes_el )
        << XmlSetTextContent( "notes" )
        << XmlSetAttribute( "opacity", 0.5 )
        << XmlSetAttribute( "y", 25 )
        << XmlSetAttribute( "font-size", 10 )
        << XmlSetAttribute( "onmouseover", "NotesMouseOver(true)" )
        << XmlSetAttribute( "onmouseout", "NotesMouseOver(false)" )
        ;

    XmlBuilder builder( doc, graph_graphic );
    builder << XmlStartElement( "g" );
    XERCES_CPP_NAMESPACE::DOMElement * notes_panel_el = builder.current_element;
    builder << XmlEndElement();

    XmlBuilder( doc, notes_panel_el )
        << XmlSetAttribute( "id", "notes_panel" )
        << XmlSetAttribute( "display", "none" )
        << XmlStartElement( "rect" )
            << XmlSetAttribute( "id", "notes_background" )
            << XmlSetAttribute( "stroke", "black" )
            << XmlSetAttribute( "fill", "white" )
            << XmlSetAttribute( "rx", "3" )
        << XmlEndElement()
        ;

    //separate the notes into individual lines
    typedef tokenizer< char_separator< char > > tokenizer;
    char_separator< char > sep( "\n" );
    tokenizer tokens( notes, sep );
    size_t spacing = 0;
    std::vector< std::string > lines;
    for( tokenizer::iterator tok_iter = tokens.begin();
        tok_iter != tokens.end();
        ++tok_iter, ++spacing )
    {
        XmlBuilder( doc, notes_panel_el )
            << XmlStartElement( "text" )
                << XmlSetTextContent( *tok_iter )
                << XmlSetAttribute( "x", 50 )
                << XmlSetAttribute( "y", 60 + spacing * 15 )
            << XmlEndElement()
            ;
    }
}


void
SvgDomDocument::add_factor(Factor * factor, const std::set<size_t> & hits)
{
    const std::string factor_id = BIO_MAKE_STRING("factor_" << factor->get_name());

    //add the text to the description
    {
        XmlBuilder builder(doc, factor_el);
        builder
            << XmlStartElement("a")
                << XmlSetAttribute("xlink:href", factor->get_link().get_url())
                << XmlSetAttribute("onmouseover", BIO_MAKE_STRING("FactorMouseOver('" << factor_id << "', true)"))
                << XmlSetAttribute("onmouseout", BIO_MAKE_STRING("FactorMouseOver('" << factor_id << "', false)"))
                << XmlStartElement("tspan")
                    << XmlSetAttribute("text-decoration", "underline")
                    << XmlSetTextContent(factor->get_name())
                << XmlEndElement()
            << XmlEndElement()
        ;
    }

    //add the definition
    {
        XmlBuilder defs_builder(doc, defs_el);
        defs_builder
            << XmlStartElement("bifa:factor_hits")
                << XmlSetAttribute("xmlns:bifa", "http://bifa.org/bifa")
                << XmlSetAttribute("id", factor_id)
                ;

        for (std::set<size_t>::const_iterator i = hits.begin();
            hits.end() != i;
            ++i)
        {
            defs_builder
                << XmlStartElement("bifa:hit")
                    << XmlSetAttribute("index", *i)
                << XmlEndElement()
                ;
        }

        defs_builder
            << XmlEndElement()
            ;
    }
}


void SvgDomDocument::add_result(
    const Hit & details,
    BiobaseTablePssmEntry & entry,
    size_t hit_idx,
    bool is_in_max_chain )
{
    /** this is what we're trying for
        <g id="hit_1">
          <polygon fill="white" points="48,97.8536 48.75,96.3536 50.25,95.6036 48.75,94.8536 48,93.3536 47.25,94.8536 45.75,95.6036 47.25,96.3536 " stroke="black" stroke-width="0.2000">
          </polygon>
          <rect id="hit_1_info_back" fill="white"/>
          <g id="hit_1_info" display="inline" opacity="0.8000" transform="translate(50 100) scale(1.5)">
            <!-- <rect fill="white" height="100%" rx="2" ry="2" stroke="black" transform="scale(1,1)" width="100%"/> -->
            <text style="font-size:10px" transform="scale(.2)">
              <a xlink:href="http://www.biobase.de/cgi-bin/biobase/transfac/9.2/bin/getTFProf.cgi?R02116">
                <tspan font-weight="bold" text-decoration="underline">V$CLOX_01</tspan>
              </a>
              <tspan dx="1.2em" fill="blue">Pathway not known</tspan>
              <tspan dy="1.2em" x="0" fill="blue" font-style="italic">Consensus: NNTATCGATTANYNW</tspan>
              <tspan dy="1.2em" x="0" fill="blue">Sequence: tgcatcgatcacann</tspan>
              <tspan dy="1.2em" x="0" fill="blue">Factors: Cutl, </tspan>
            </text>
              <svg height="6.0000" preserveAspectRatio="xMinYMin" viewBox="0 0 1100 200" width="33" x="0.0000" y="-9.0000">
            ... SVG logo ...
              </svg>
          </g>
        </g>
    */

    const BiobaseTablePssmEntry * pssm_entry = BiobaseDb::singleton().get_pssm_entry(entry.get_link());

    TableLink pathway_link =
        PathwayAssociations::singleton().get_most_significant_pathway_for(entry.get_link());

    const string pathway_description =
        PATHWAY_DATA != pathway_link.table_id
            ? "Pathway not known"
            : BiobaseDb::singleton().get_entry<PATHWAY_DATA>(pathway_link)->get_name();

    const Shape shape =
        MATRIX_DATA == entry.get_link().table_id
            ? STAR
            : CIRCLE;

    Pssm pssm;
    string iupac;
    if (MATRIX_DATA == entry.get_link().table_id)
    {
        Matrix * matrix = BiobaseDb::singleton().get_entry<MATRIX_DATA>(entry.get_link());
        copy(matrix->consensus_matrix.begin(), matrix->consensus_matrix.end(), inserter(iupac, iupac.end()));
        pssm = make_pssm(matrix);
    }
    else
    {
        Site * site = BiobaseDb::singleton().get_entry<SITE_DATA>(entry.get_link());
        iupac += site->sequence;
        pssm = make_pssm(site);
    }

    add_hit(
        entry.get_name(),
        pssm,
        details.position,
        details.complement,
        details.score,
        hit_idx,
        is_in_max_chain,
        shape,
        iupac,
        pssm_entry->get_factors(),
        pathway_link,
        pathway_description,
        entry.get_link().get_url() );
}

void SvgDomDocument::add_hit(
    const std::string & pssm_name,
    const Pssm & pssm,
    unsigned position,
    bool complement,
    float_t p_binding,
    size_t hit_idx,
    bool is_in_max_chain,
    Shape shape,
    const std::string & iupac,
    const FactorLinkList & factors,
    TableLink pathway_link,
    const std::string & pathway_description,
    const std::string & url )
{
    const size_t mid_position = position + pssm.size() / 2;

    DOMElement * hit_el = doc->createElement(XStr("g"));
    set_attribute(hit_el, "id", BIO_MAKE_STRING("hit_" << hit_idx));
    set_attribute(
        hit_el,
        "transform",
        BIO_MAKE_STRING(
            "translate(" << mid_position
            << " "
            << float_t(100) - (p_binding - min_threshold) / get_range() * float_t(100)
            << ")"));

    //add the underline graphic
    {
        /**
          <g id="hit_36_underlines">
            <path d="M -5 -86.0656 L 5 -86.0656" stroke="orange" stroke-width="1" opacity="0.4"/>
          </g>
          */
        DOMElement * underline_el = doc->createElement(XStr("g"));
        set_attribute(underline_el, "id", BIO_MAKE_STRING("hit_" << hit_idx << "_underlines"));
        set_attribute(underline_el, "display", "none");
        for (size_t i = 0; i != 6; ++i)
        {
            //const float_t t = min_threshold + get_range() * float_t(i) / float_t(5);
            const float_t y =
                ((p_binding - min_threshold) / get_range() * float_t(100) - float_t(100))
                + float_t(i * 100 / 5)
                + float_t( 0.5 )
                ;

            DOMElement * path_el = doc->createElement(XStr("path"));
            set_attribute(
                path_el,
                "d",
                BIO_MAKE_STRING(
                    "M " << int(position) - int(mid_position) << " " << y
                    << " L " << int(position) - int(mid_position) + pssm.size() << " " << y));
            set_attribute(path_el, "stroke", "orange");
            set_attribute(path_el, "stroke-width", "1");
            set_attribute(path_el, "opacity", "0.5");
            underline_el->appendChild(path_el);
        }
        hit_el->appendChild(underline_el);
    }

    //add the background
    {
        const string fill =
            PATHWAY_DATA != pathway_link.table_id
                ? "white"
                : get_pathway_colour_for(pathway_link);
        
        DOMElement * shape_el =
            create_shape(
                CIRCLE,
                coord_t(0, 0),
                7.5f,
                0.1f,
                fill,
                "black",
                "",
                "");
        set_attribute(shape_el, "opacity", "0.9");
        set_attribute(shape_el, "display", "none");
        set_attribute(shape_el, "id", BIO_MAKE_STRING("hit_" << hit_idx << "_bg"));
        hit_el->appendChild(shape_el);
    }

    //add the max chain background
    if( is_in_max_chain )
    {
        const string fill =
            PATHWAY_DATA != pathway_link.table_id
                ? "white"
                : get_pathway_colour_for(pathway_link);
        
        DOMElement * shape_el =
            create_shape(
                STAR,
                coord_t(0, 0),
                7.5f,
                0.1f,
                fill,
                "black",
                "",
                "");
        set_attribute(shape_el, "opacity", "0.5");
        set_attribute(shape_el, "pointer-events", "none");
        //set_attribute(shape_el, "display", "none");
        set_attribute(shape_el, "id", BIO_MAKE_STRING("hit_" << hit_idx << "_mc"));
        hit_el->appendChild(shape_el);
    }

    //add the shape
    {
        const string fill =
            PATHWAY_DATA != pathway_link.table_id
                ? "white"
                : get_pathway_colour_for(pathway_link);
        
        DOMElement * shape_el =
            create_shape(
                shape,
                coord_t(0, 0),
                1.5f,
                0.2f,
                fill,
                "black",
                pssm_name,
                pathway_description);
        hit_el->appendChild(shape_el);
    }

    //add the title
    if (show_titles)
    {
        /*
        <text style="font-size:7px" dx="5">
            V$GATA2_02
        </text>
        */
        DOMElement * text_el = doc->createElement(XStr("text"));
        set_attribute(text_el, "font-size", "6px");
        set_attribute(text_el, "dx", "5");
        text_el->setTextContent(XStr(pssm_name));
        hit_el->appendChild(text_el);
    }

    //add the info details
    {
        DOMElement * info_el = doc->createElement(XStr("g"));
        set_attribute(info_el, "id", BIO_MAKE_STRING("hit_" << hit_idx << "_info"));
        set_attribute(info_el, "display", "none");
        set_attribute(info_el, "opacity", "0.7");
        set_attribute(info_el, "transform", "translate(3 0) scale(1.5)");

        //add the background for the info
        {
            DOMElement * back_el = doc->createElement(XStr("rect"));
            set_attribute(back_el, "id", BIO_MAKE_STRING("hit_" << hit_idx << "_info_back"));
            set_attribute(back_el, "fill", "white");
            set_attribute(back_el, "rx", "1");
            info_el->appendChild(back_el);
        }

        //add the text to the info
        {
            DOMElement * text_el = doc->createElement(XStr("text"));
            set_attribute(text_el, "font-family", "Arial");
            set_attribute(text_el, "style", "font-size:10px");
            set_attribute(text_el, "transform", "scale(.2)");

            //the name
            {
                DOMElement * tspan_el = doc->createElement(XStr("tspan"));
                set_attribute(tspan_el, "font-weight", "bold");
                tspan_el->setTextContent(XStr(pssm_name));
                if( url != "" )
                {
                    set_attribute(tspan_el, "text-decoration", "underline");
                    DOMElement * hyperlink_el = doc->createElement(XStr("a"));
                    set_attribute(hyperlink_el, "xlink:href", url);
                    hyperlink_el->appendChild(tspan_el);
                    text_el->appendChild(hyperlink_el);
                } else {
                    text_el->appendChild(tspan_el);
                }
            }

            //the pathway
            {
                DOMElement * tspan_el = doc->createElement(XStr("tspan"));
                set_attribute(tspan_el, "dx", "1.2em");
                set_attribute(tspan_el, "fill", "blue");
                tspan_el->setTextContent(XStr(pathway_description));

                if (PATHWAY_DATA == pathway_link.table_id)
                {
                    DOMElement * hyperlink_el = doc->createElement(XStr("a"));
                    set_attribute(tspan_el, "text-decoration", "underline");
                    set_attribute(hyperlink_el, "xlink:href", BiobaseDb::singleton().get_entry<PATHWAY_DATA>(pathway_link)->get_link().get_url());
                    hyperlink_el->appendChild(tspan_el);
                    text_el->appendChild(hyperlink_el);
                }
                else
                {
                    text_el->appendChild(tspan_el);
                }
            }

            //consensus
            if( iupac != "" )
            {
                DOMElement * tspan_el = doc->createElement(XStr("tspan"));
                set_attribute(tspan_el, "dy", "1.2em");
                set_attribute(tspan_el, "x", "0");
                set_attribute(tspan_el, "fill", "blue");
                tspan_el->setTextContent(XStr(BIO_MAKE_STRING( "Consensus: " << iupac )));

                text_el->appendChild(tspan_el);
            }

            //the sequence
            string sequence;
            if (! complement)
            {
                sequence = seq.substr(position, pssm.size());
            }
            else
            {
                reverse_complement(seq.substr(position, pssm.size()), inserter(sequence, sequence.begin()));
            }
            const string sequence_text = "Sequence: " + sequence;
            {
                DOMElement * tspan_el = doc->createElement(XStr("tspan"));
                set_attribute(tspan_el, "dy", "1.2em");
                set_attribute(tspan_el, "x", "0");
                set_attribute(tspan_el, "fill", "blue");
                tspan_el->setTextContent(XStr(sequence_text));

                text_el->appendChild(tspan_el);
            }

            //the position
            {
                DOMElement * tspan_el = doc->createElement(XStr("tspan"));
                set_attribute(tspan_el, "dy", "1.2em");
                set_attribute(tspan_el, "x", "0");
                set_attribute(tspan_el, "fill", "blue");

                stringstream str_stream;
                if (! complement)
                {
                    str_stream
                        << "Position: ["
                        << position
                        << ","
                        << position + pssm.size()
                        << "]";
                }
                else
                {
                    str_stream
                        << "Position: ["
                        << position + pssm.size()
                        << ","
                        << position
                        << "]";
                }
                tspan_el->setTextContent(XStr(str_stream.str()));

                text_el->appendChild(tspan_el);
            }

            //factors
            {
                set<string> factor_names;
                for (FactorLinkList::const_iterator i = factors.begin();
                    factors.end() != i;
                    ++i)
                {
                    factor_names.insert(i->get()->name);
                }
                stringstream factors_stream;
                factors_stream << "Factors: ";
                copy(factor_names.begin(), factor_names.end(), ostream_iterator<string>(factors_stream, ", "));
                const string factors = factors_stream.str();

                DOMElement * tspan_el = doc->createElement(XStr("tspan"));
                set_attribute(tspan_el, "dy", "1.2em");
                set_attribute(tspan_el, "x", "0");
                set_attribute(tspan_el, "fill", "blue");
                tspan_el->setTextContent(XStr(factors));

                text_el->appendChild(tspan_el);
            }
            info_el->appendChild(text_el);

            //svg logo
            {
                DOMElement * svg_logo_el = create_svg_pssm_logo(
                    pssm,
                    doc,
                    sequence);
                set_attribute(svg_logo_el, "width", float_t( sequence.size() * 3 ));
                set_attribute(svg_logo_el, "height", float_t( 6 ));
                set_attribute(svg_logo_el, "preserveAspectRatio", "xMinYMin");
                set_attribute(svg_logo_el, "y", "-9");

                info_el->appendChild(svg_logo_el);
            }
        }

        hit_el->appendChild(info_el);
    }

    points_graphic->appendChild(hit_el);
}


DOMElement * SvgDomDocument::create_shape(
    Shape shape,
    coord_t position,
    float_t size,
    float_t stroke_width,
    const XStr & fill,
    const XStr & stroke,
    const string & title,
    const string & desc)
{
    DOMElement * element;
    switch(shape) {
        case CIRCLE:
            {
                element = doc->createElement(XStr("circle"));
                set_attribute(element, "cx", position.first);
                set_attribute(element, "cy", position.second);
                set_attribute(element, "r", size);
            }
            break;

        case SQUARE:
            {
                element = doc->createElement(XStr("rect"));
                set_attribute(element, "x", position.first - size);
                set_attribute(element, "y", position.second - size);
                set_attribute(element, "width", 2 * size);
                set_attribute(element, "height", 2 * size);
            }
            break;

        case VERT_ELLIPSE:
            {
                element = doc->createElement(XStr("ellipse"));
                set_attribute(element, "cx", position.first);
                set_attribute(element, "cy", position.second);
                set_attribute(element, "rx", size / 2);
                set_attribute(element, "ry", 3 * size / 2);
            }
            break;

        case HORIZ_ELLIPSE:
            {
                element = doc->createElement(XStr("ellipse"));
                set_attribute(element, "cx", position.first);
                set_attribute(element, "cy", position.second);
                set_attribute(element, "rx", 3 * size / 2);
                set_attribute(element, "ry", size / 2);
            }
            break;

        case UP_TRIANGLE:
            {
                const coord_t min (position.first - size, position.second - size);
                const coord_t max (position.first + size, position.second + size);
                vector<coord_t> points;
                points.push_back(coord_t(position.first, max.second));
                points.push_back(coord_t(min.first, min.second));
                points.push_back(coord_t(max.first, min.second));
                element = create_polygon(points.begin(), points.end());
            }
            break;

        case DOWN_TRIANGLE:
            {
                const coord_t min (position.first - size, position.second - size);
                const coord_t max (position.first + size, position.second + size);
                vector<coord_t> points;
                points.push_back(coord_t(position.first, min.second));
                points.push_back(coord_t(min.first, max.second));
                points.push_back(coord_t(max.first, max.second));
                element = create_polygon(points.begin(), points.end());
            }
            break;

        case RIGHT_TRIANGLE:
            {
                const coord_t min (position.first - size, position.second - size);
                const coord_t max (position.first + size, position.second + size);
                vector<coord_t> points;
                points.push_back(coord_t(max.first, position.second));
                points.push_back(coord_t(min.first, max.second));
                points.push_back(coord_t(min.first, min.second));
                element = create_polygon(points.begin(), points.end());
            }
            break;

        case LEFT_TRIANGLE:
            {
                const coord_t min (position.first - size, position.second - size);
                const coord_t max (position.first + size, position.second + size);
                vector<coord_t> points;
                points.push_back(coord_t(min.first, position.second));
                points.push_back(coord_t(max.first, max.second));
                points.push_back(coord_t(max.first, min.second));
                element = create_polygon(points.begin(), points.end());
            }
            break;

        case STAR:
            {
                const coord_t outer_min(float_t( position.first - size * 1.5 ), float_t( position.second - size * 1.5 ));
                const coord_t outer_max(float_t( position.first + size * 1.5 ), float_t( position.second + size * 1.5 ));
                const coord_t inner_min(float_t( position.first - size / 2 ), float_t( position.second - size / 2 ));
                const coord_t inner_max(float_t( position.first + size / 2 ), float_t( position.second + size / 2 ));
                vector<coord_t> points;
                points.push_back(coord_t(position.first, outer_max.second));
                points.push_back(coord_t(inner_max.first, inner_max.second));
                points.push_back(coord_t(outer_max.first, position.second));
                points.push_back(coord_t(inner_max.first, inner_min.second));
                points.push_back(coord_t(position.first, outer_min.second));
                points.push_back(coord_t(inner_min.first, inner_min.second));
                points.push_back(coord_t(outer_min.first, position.second));
                points.push_back(coord_t(inner_min.first, inner_max.second));
                element = create_polygon(points.begin(), points.end());
            }
            break;

        case ROTATED_STAR:
            {
                const coord_t outer_min(float_t( position.first - size * 1.5 ), float_t( position.second - size * 1.5 ));
                const coord_t outer_max(float_t( position.first + size * 1.5 ), float_t( position.second + size * 1.5 ));
                const coord_t inner_min(float_t( position.first - size / 2 ), float_t( position.second - size / 2 ));
                const coord_t inner_max(float_t( position.first + size / 2 ), float_t( position.second + size / 2 ));
                vector<coord_t> points;
                points.push_back(coord_t(outer_max.first, outer_max.second));
                points.push_back(coord_t(inner_max.first, position.second));
                points.push_back(coord_t(outer_max.first, outer_min.second));
                points.push_back(coord_t(position.first, inner_min.second));
                points.push_back(coord_t(outer_min.first, outer_min.second));
                points.push_back(coord_t(inner_min.first, position.second));
                points.push_back(coord_t(outer_min.first, outer_max.second));
                points.push_back(coord_t(position.first, inner_max.second));
                element = create_polygon(points.begin(), points.end());
            }
            break;

        default:
            throw std::logic_error( "Unsupported shape" );
    }

    set_attribute(element, "fill", fill);
    set_attribute(element, "stroke", stroke);
    set_attribute(element, "stroke-width", stroke_width);

    if (desc != "")
    {
        element->appendChild(create_text_element(doc, "desc", desc));
    }
    if (title != "")
    {
        element->appendChild(create_text_element(doc, "title", title));
    }

    return element;
}

DOMElement * SvgDomDocument::create_label(
    coord_t position,
    const XStr & text,
    const XStr & style)
{
    DOMElement * label = doc->createElement(XStr("tspan"));
    set_attribute(label, position);
    set_attribute(label, "style", style);
    label->setTextContent(text);
    return label;
}

DOMElement * SvgDomDocument::create_line(coord_t start, coord_t end)
{
    DOMElement * line = doc->createElement(XStr("line"));
    set_attribute(line, "x1", start.first);
    set_attribute(line, "y1", start.second);
    set_attribute(line, "x2", end.first);
    set_attribute(line, "y2", end.second);
    return line;
}

/** scale the significance from [threshold, 1] to [0,-1] */
float_t SvgDomDocument::scale_significance(float_t significance) const
{
    if (significance < min_threshold) {
        throw std::logic_error( "significance too small" );
    }
    if (significance > 1) {
        throw std::logic_error( "significance too big" );
    }
    return (significance - min_threshold) / get_range();
}

float_t SvgDomDocument::get_range() const
{
    return max_threshold - min_threshold;
}

coord_t SvgDomDocument::get_graph_coord(size_t idx, float_t significance) const
{
    if (idx > num_bases) {
        throw std::logic_error( "idx too high" );
    }

    const float_t scaled_significance = scale_significance(significance);

    return
        coord_t(
            float_t( idx * 100 ) / num_bases,
            float_t( 60 * (1.0 - scaled_significance)));
}

const string & SvgDomDocument::get_pathway_colour_for(const TableLink & pathway_link)
{
    assert(colours.size() > interesting_pathways.size()); //make sure we have enough (plus one for not identified pathways)

    //find the index of the pathway in interesting_pathways
    size_t index = 0;

    TableLinkVec::const_iterator interesting_pathway
        = find(interesting_pathways.begin(), interesting_pathways.end(), pathway_link);
    if (interesting_pathways.end() != interesting_pathway)
    {
        index = interesting_pathway - interesting_pathways.begin() + 1;
    }

    return colours[index];
}



string num2str(double num)
{
    stringstream ss;
    ss << setprecision(4) << fixed << num;
    return ss.str();
}

string num2str(float_t num)
{
    stringstream ss;
    ss << setprecision(4) << fixed << num;
    return ss.str();
}

string num2str(size_t num)
{
    stringstream ss;
    ss << setprecision(4) << fixed << num;
    return ss.str();
}

void set_attribute(DOMElement * el, const XStr & attr, const XStr & value) {
    el->setAttribute(attr, value);
}

void set_attribute(DOMElement * el, const XStr & attr, double value) {
    el->setAttribute(attr, XStr(num2str(value)));
}

void set_attribute(DOMElement * el, const XStr & attr, float_t value) {
    el->setAttribute(attr, XStr(num2str(value)));
}

void set_attribute(DOMElement * el, coord_t value) {
    el->setAttribute(XStr("x"), XStr(num2str(value.first)));
    el->setAttribute(XStr("y"), XStr(num2str(value.second)));
}


BIO_NS_END
