/* Copyright John Reid 2007
*/

#include "bio-pch.h"


#include "bio/defs.h"
#include "bio/svg_match.h"

#include <xercesc/dom/DOMCDATASection.hpp>
XERCES_CPP_NAMESPACE_USE



BIO_NS_START

using XERCES_CPP_NAMESPACE::DOMDocument;

void
add_tooltip_support(
	DOMDocument * doc,
	DOMElement * svg_root)
{
	set_attribute(svg_root, "onload", "Init(evt)");
	set_attribute(svg_root, "onmousemove", "GetTrueCoords(evt); ShowTooltip(evt, true)");
	set_attribute(svg_root, "onmouseout", "ShowTooltip(evt, false)");
	set_attribute(svg_root, "onclick", "OnClick(evt)");
	set_attribute(svg_root, "onload", "Init(evt)");

	/** <script xlink:href="mouse_over_effects.js" type="text/ecmascript" /> */
	DOMElement * script = doc->createElement(XStr("script"));
	set_attribute(script, "type", "text/ecmascript");
	set_attribute(script, "xlink:href", "bifa.js");
	svg_root->appendChild(script);

/**
	   <g id="ToolTip" opacity="0.8" display="none" pointer-events="none" transform="translate(124,32)">
      <rect id="tipbox" x="0" y="5" width="85.5554" height="35.3797" rx="2" ry="2" fill="white" stroke="black" transform="scale(1,1)"/>
      <text id="tipText" x="5" y="20" font-size="12" startOffset="0" transform="scale(1,1)">
         <tspan id="tipTitle" x="5" font-weight="bold">Ellipse</tspan>
         <tspan id="tipDesc" x="5" dy="1.2em" fill="blue">A black ellipse</tspan>
      </text>
   </g>
*/
	DOMElement * tooltip_graphic = doc->createElement(XStr("g"));
	set_attribute(tooltip_graphic, "id", "ToolTip");
	set_attribute(tooltip_graphic, "opacity", 0.8f);
	set_attribute(tooltip_graphic, "display", "none");
	set_attribute(tooltip_graphic, "pointer-events", "none");
	set_attribute(tooltip_graphic, "transform", "translate(0,0)");
	svg_root->appendChild(tooltip_graphic);

	DOMElement * tipbox = doc->createElement(XStr("rect"));
	set_attribute(tipbox, "id", "tipbox");
	set_attribute(tipbox, "x", "5");
	set_attribute(tipbox, "y", "5");
	set_attribute(tipbox, "width", "80");
	set_attribute(tipbox, "height", "30");
	set_attribute(tipbox, "rx", "2");
	set_attribute(tipbox, "ry", "2");
	set_attribute(tipbox, "fill", "white");
	set_attribute(tipbox, "stroke", "black");
	set_attribute(tipbox, "transform", "scale(1,1)");
	tooltip_graphic->appendChild(tipbox);

	DOMElement * tip_text = doc->createElement(XStr("text"));
	set_attribute(tip_text, "id", "tipText");
	set_attribute(tip_text, "x", "10");
	set_attribute(tip_text, "y", "20");
	set_attribute(tip_text, "transform", "scale(1,1)");
	tooltip_graphic->appendChild(tip_text);

	{
		DOMElement * el = doc->createElement(XStr("tspan"));
		set_attribute(el, "id", "tipTitle");
		set_attribute(el, "x", "10");
		set_attribute(el, "font-weight", "bold");
		el->setTextContent(XStr("Title"));
		tip_text->appendChild(el);
	}

	{
		DOMElement * el = doc->createElement(XStr("tspan"));
		set_attribute(el, "id", "tipDesc");
		set_attribute(el, "dx", "1.2em");
		set_attribute(el, "fill", "blue");
		el->setTextContent(XStr("Description"));
		tip_text->appendChild(el);
	}

	{
		DOMElement * el = doc->createElement(XStr("tspan"));
		set_attribute(el, "id", "tipIupac");
		set_attribute(el, "dx", "1.2em");
		set_attribute(el, "fill", "blue");
		set_attribute(el, "font-style", "italic");
		el->setTextContent(XStr("Iupac"));
		tip_text->appendChild(el);
	}

	{
		DOMElement * el = doc->createElement(XStr("tspan"));
		set_attribute(el, "id", "tipSeq");
		set_attribute(el, "dx", "1.2em");
		set_attribute(el, "fill", "blue");
		el->setTextContent(XStr("Seq"));
		tip_text->appendChild(el);
	}

	{
		DOMElement * el = doc->createElement(XStr("tspan"));
		set_attribute(el, "id", "tipFactors");
		set_attribute(el, "x", "10");
		set_attribute(el, "dy", "1.2em");
		set_attribute(el, "fill", "blue");
		el->setTextContent(XStr("Factors"));
		tip_text->appendChild(el);
	}
}

BIO_NS_END

