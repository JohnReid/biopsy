/* Copyright John Reid 2007-2010
*/

#ifndef BIO_SVG_MATCH_H_
#define BIO_SVG_MATCH_H_


#include "bio/defs.h"
#include "bio/x_str.h"
#include "bio/svg.h"
#include "bio/common.h"
#include "bio/match_hit.h"
#include "bio/sequence.h"
#include "bio/factor.h"
#include "bio/run_match.h"

#include <boost/filesystem/path.hpp>

#include <xercesc/dom/DOM.hpp>

#include <set>
#include <string>
#include <sstream>
#include <utility>



BIO_NS_START

void add_tooltip_support(
	XERCES_CPP_NAMESPACE::DOMDocument * doc,
	XERCES_CPP_NAMESPACE::DOMElement * svg_root);



class SvgDomDocument
{
public:
	SvgDomDocument(
		size_t num_bases,
		float_t min_threshold,
		float_t max_threshold,
		const std::string & title,
		const seq_t & seq,
		bool show_titles);

	coord_t get_graph_coord(size_t idx, float_t significance) const;
	float_t get_range() const;
	float_t scale_significance(float_t significance) const;

	enum Shape {
		CIRCLE,
		UP_TRIANGLE,
		DOWN_TRIANGLE,
		RIGHT_TRIANGLE,
		LEFT_TRIANGLE,
		STAR,
		ROTATED_STAR,
		SQUARE,
		VERT_ELLIPSE,
		HORIZ_ELLIPSE,
		NUM_SHAPES
	};

public:
	XERCES_CPP_NAMESPACE::DOMDocument * doc;
	XERCES_CPP_NAMESPACE::DOMElement * doc_root;
	XERCES_CPP_NAMESPACE::DOMElement * graph_graphic;
	XERCES_CPP_NAMESPACE::DOMElement * graph_area_graphic;
	XERCES_CPP_NAMESPACE::DOMElement * points_graphic;
	XERCES_CPP_NAMESPACE::DOMElement * defs_el;
	XERCES_CPP_NAMESPACE::DOMElement * factor_el;
	XERCES_CPP_NAMESPACE::DOMElement * notes_el;

	float_t height;
	float_t width;
	size_t num_bases;
	float_t min_threshold;
	float_t max_threshold;
	seq_t seq;
	bool show_titles;

	void add_result(
		const Hit & details,
		BiobaseTablePssmEntry & entry,
		size_t hit_idx,
		bool is_in_max_chain );

	void add_hit(
		const std::string & pssm_name,
		const Pssm & pssm,
		unsigned position,
		bool complement,
		float_t p_binding,
		size_t hit_idx,
		bool is_in_max_chain,
		Shape shape,
		const std::string & iupac = "",
		const FactorLinkList & factors = FactorLinkList(),
		TableLink pathway_link = TableLink(),
		const std::string & pathway_description = "No known pathway",
		const std::string & url = "" );

	void add_factor(Factor * factor, const std::set<size_t> & hits);

	void add_notes( const std::string & notes );

protected:

	void build_graph(const std::string & title);
	void build_pathways();

	XERCES_CPP_NAMESPACE::DOMElement * create_line(coord_t start, coord_t end);
	XERCES_CPP_NAMESPACE::DOMElement * create_label(
		coord_t position,
		const XStr & text,
		const XStr & style);

	XERCES_CPP_NAMESPACE::DOMElement * create_shape(
		Shape shape,
		coord_t position,
		float_t size,
		float_t stroke_width,
		const XStr & fill,
		const XStr & stroke,
		const std::string & desc,
		const std::string & title);

	template <class CoordIt>
	XERCES_CPP_NAMESPACE::DOMElement * create_polygon(CoordIt begin, CoordIt end)
	{
		XERCES_CPP_NAMESPACE::DOMElement * element = doc->createElement(XStr("polygon"));
		std::stringstream points;
		while (begin != end) {
			points << begin->first << "," << begin->second << " ";
			++begin;
		}
		set_attribute(element, "points", points.str().c_str());
		return element;
	}

	static const std::string & get_pathway_colour_for(const TableLink & link);
};


struct BuildSvgArgs
{
	std::string file;
	std::string title;
	float_t max_threshold;
	float_t min_threshold;
	unsigned max_num_factors;
	bool show_labels;
	bool open_file;
	std::string notes;

	BuildSvgArgs(
		const std::string & file = "matches.svg",
		const std::string & title = "BiFa Analysis",
		float_t max_threshold = 0.0,
		float_t min_threshold = 0.05,
		unsigned max_num_factors = 12,
		bool show_labels = false,
#ifdef WIN32
		bool open_file = true,
#else //WIN32
		bool open_file = false,
#endif //WIN32
		const std::string & notes = ""
	);

	void add_options(boost::program_options::options_description & options);

	/** Run the algorithm with the arguments. */
	void
	run_build_svg(
		const seq_t & seq,
		match_result_vec_t & results);
};

std::ostream &
operator<<(std::ostream & os, const BuildSvgArgs & args);

void
build_svg(
	const boost::filesystem::path & file,
	const std::string & title,
	const seq_t & seq,
	float_t min_threshold,
	match_result_vec_t & results,
	size_t max_num_factors,
	bool show_labels,
	bool open_file,
	match_result_vec_t * max_chain = 0,
	const std::string & notes = "",
	float max_threshold = 0.0 );



BIO_NS_END



#endif //BIO_SVG_MATCH_H_

