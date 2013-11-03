/* Copyright John Reid 2007
*/

#include "bio-pch.h"


#include "bio/defs.h"

#include "bio/svg_pssm_logo.h"
#include "bio/x_str.h"
#include "bio/svg.h"

using namespace boost;

XERCES_CPP_NAMESPACE_USE

#include <vector>
#include <limits>
using namespace std;


BIO_NS_START

using XERCES_CPP_NAMESPACE::DOMDocument;
using XERCES_CPP_NAMESPACE::DOMElement;

namespace impl
{
	static const char base_chars[] = "ACGT";
	static const char * base_strs[] = {
		"A",
		"C",
		"G",
		"T"
	};
	static const char * transforms[] = {
		"scale(1.5,1.4)",
		"scale(1.58,1.35) translate(-5,0)",
		"scale(1.5,1.35) translate(-5,0)",
		"scale(1.75,1.4) translate(-2,0)",
	};
	static const char * styles[] = {
		"text-anchor:start; font-family: arial; fill:chartreuse; font-size:100px",
		"text-anchor:start; font-family: arial; fill:orange; font-size:100px; baseline-shift:1%",
		"text-anchor:start; font-family: arial; fill:cyan; font-size:100px; baseline-shift:1%",
		"text-anchor:start; font-family: arial; fill:red; font-size:100px",
	};
}

using namespace impl;

DOMElement *
add_logo_defs(
	DOMDocument * doc)
{
	//if we can't find the right defs in the document then add them
	DOMElement * defs_el = doc->createElement(XStr("defs"));
	for (size_t i = 0; i != 4; ++i)
	{
		DOMElement * g_el = doc->createElement(XStr("g"));
		set_attribute(g_el, "id", base_strs[i]);
		set_attribute(g_el, "transform", transforms[i]);

		DOMElement * text_el = doc->createElement(XStr("text"));
		set_attribute(text_el, "style", styles[i]);
		text_el->setTextContent(XStr(base_strs[i]));
		g_el->appendChild(text_el);

		defs_el->appendChild(g_el);
	}
	doc->getDocumentElement()->appendChild(defs_el);
	return defs_el;
}

DOMElement *
create_svg_pssm_logo(
	const Pssm & pssm,
	DOMDocument * doc,
	const seq_t & seq)
{
	//calculate the heights of the bases
	typedef vector<float_t> base_heights_t;
	typedef vector<base_heights_t> base_heights_vec_t;
	base_heights_vec_t base_heights;
	//see http://www.lecb.ncifcrf.gov/~toms/paper/logopaper/paper/index.html for description of the following
	for (Pssm::const_iterator l = pssm.begin();
		pssm.end() != l;
		++l)
	{
		float_t H_l = 0.0;
		for (size_t i = 0; i != 4; ++i) //for each nucleotide
		{
			const float_t f_b_l = l->get_freq(base_chars[i]);
			if (f_b_l > 0.0)
			{
				H_l -= f_b_l * std::log(f_b_l);
			}
		}

		base_heights.push_back(vector<float>());
		const float_t R_sequence_l = float_t(2.0) - float_t(H_l);
		for (size_t i = 0; i != 4; ++i) //for each nucleotide
		{
			base_heights.rbegin()->push_back(l->get_freq(base_chars[i]) * R_sequence_l);
		}
	}

	size_t pos = 0;
	DOMElement * result = doc->createElement(XStr("svg"));
	set_attribute(
		result,
		"viewBox",
		BIO_MAKE_STRING("0 0 " << (pssm.size() * 100) << " 200"));
	for (base_heights_vec_t::iterator l = base_heights.begin();
		base_heights.end() != l;
		++l, ++pos)
	{
		DOMElement * base_el = doc->createElement(XStr("g"));
		set_attribute(
			base_el,
			"transform",
			BIO_MAKE_STRING("translate(" << pos * 100 << " 200)"));

		float_t vert_pos = 0.0;
		for (size_t i = 0; i != 4; ++i) //for each nucleotide
		{
			//find the smallest base
			base_heights_t::iterator b = min_element(l->begin(), l->end());
			const char base = base_chars[b - l->begin()];
			const char lower_base = base - ('A' - 'a'); //the base in lower case

			if (*b > 0.0)
			{
				DOMElement * nucleo_el = doc->createElement(XStr("g"));
				set_attribute(
					nucleo_el,
					"transform",
					BIO_MAKE_STRING(
						"translate(0 " << vert_pos << ") scale(1 " << *b << ")"));

				//if we have a sequence and it matches this base, highlight it with a background rect
				if (pos < seq.size() && (seq[pos] == base || seq[pos] == lower_base))
				{
					//<rect x="0" y="-100" width="100" height="100" fill="gray"/>
					DOMElement * rect_el = doc->createElement(XStr("rect"));
					set_attribute(rect_el, "rx", 10);
					set_attribute(rect_el, "y", -100);
					set_attribute(rect_el, "width", 100);
					set_attribute(rect_el, "height", 100);
					set_attribute(rect_el, "fill", "gray");
					nucleo_el->appendChild(rect_el);
				}

				DOMElement * use_el = doc->createElement(XStr("use"));
				set_attribute(
					use_el,
					"xlink:href",
					BIO_MAKE_STRING("#" << base));
				nucleo_el->appendChild(use_el);

				base_el->appendChild(nucleo_el);

				//adjust the vertical position
				vert_pos -= *b * 100;
			}

			//make sure this is not the min next time around
			*b = numeric_limits<float_t>::max();
		}

		result->appendChild(base_el);
	}

	return result;
}


BIO_NS_END

