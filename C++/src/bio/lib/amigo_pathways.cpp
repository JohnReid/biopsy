/* Copyright John Reid 2007
*/

#include "bio-pch.h"


#include "bio/defs.h"
#include "bio/amigo_pathways.h"
#include "bio/x_str.h"
USING_BIO_NS

#include <xercesc/sax2/SAX2XMLReader.hpp>
#include <xercesc/sax2/XMLReaderFactory.hpp>
#include <xercesc/sax2/DefaultHandler.hpp>
#include <xercesc/util/XMLString.hpp>
XERCES_CPP_NAMESPACE_USE

#include <iostream>
#include <stack>
#include <string>
#include <algorithm>
using namespace std;


#define TERM_NAME "term"
#define GENE_PRODUCT_NAME "gene_product"
#define NAME_NAME "name"
#define DATABASE_SYMBOL_NAME "database_symbol"
#define REFERENCE_NAME "reference"


BIO_NS_START

const char * pathway_xml_files [] = {
	"C:\\data\\amigo\\signaling\\cytokine-and-chemokine-mediated.xml",
	"C:\\data\\amigo\\signaling\\enzyme-linked-receptor-protein-signaling-pathway.xml",
	"C:\\data\\amigo\\signaling\\ER-nuclear.xml",
	"C:\\data\\amigo\\signaling\\glutamate.xml",
	"C:\\data\\amigo\\signaling\\notch.xml",
	"C:\\data\\amigo\\signaling\\receptor-guanylyl-cyclase.xml",
	"C:\\data\\amigo\\signaling\\transforming-growth-factor-beta-receptor.xml",
	"C:\\data\\amigo\\signaling\\Wnt-receptor.xml",
};

AmigoPathwaySet pathways;

bool
AmigoPathway::contains_database_ref(const DatabaseRef & database_ref) const
{
	return find(database_refs.begin(), database_refs.end(), database_ref) != database_refs.end();
}

void
AmigoPathwaySet::parse_default_xml_files() {
	parse_xml_files(pathway_xml_files, pathway_xml_files + sizeof(pathway_xml_files) / sizeof(const char *));
}

class MySAX2Handler : public DefaultHandler {
public:
	MySAX2Handler(AmigoPathwayPtr pathway)
		: pathway(pathway)
	{ }

	enum State {
		TERM_STATE,
		GENE_PRODUCT_STATE,
		NAME_STATE,
		DATABASE_SYMBOL_STATE,
		REFERENCE_STATE,
		NO_STATE,
	};

	static State get_state(const XMLCh * const localname) {
		const string name = (const char *) StrX(localname);
		if (name == TERM_NAME) {
			return TERM_STATE;
		} else if (name == GENE_PRODUCT_NAME) {
			return GENE_PRODUCT_STATE;
		} else if (name == NAME_NAME) {
			return NAME_STATE;
		} else if (name == DATABASE_SYMBOL_NAME) {
			return DATABASE_SYMBOL_STATE;
		} else if (name == REFERENCE_NAME) {
			return REFERENCE_STATE;
		} else {
			return NO_STATE;
		}
	}

	void startElement(
		const   XMLCh* const    uri,
		const   XMLCh* const    localname,
		const   XMLCh* const    qname,
		const   Attributes&     attrs)
	{
		const State state = get_state(localname);
		state_stack.push(state);
	}

	void endElement(
		const XMLCh * const uri,
		const XMLCh *const localname,
		const XMLCh *const qname)
	{
		assert(! state_stack.empty());
		const State state = get_state(localname);
		assert(state_stack.top() == state);

		switch(state) {
		case GENE_PRODUCT_STATE:
		{
			Database db = parse_database_string(database_symbol);
			if (UNKNOWN_DB != db) {
				pathway->database_refs.push_back( DatabaseRef( db, reference, 0 ) );
			} else {
				throw std::logic_error( "Unknown database" );
			}
			break;
		}
		default:
			break;
		}
		state_stack.pop();
	}

	void characters(const XMLCh * const chars, const unsigned int length)
	{
		switch (state_stack.top()) {
		case NAME_STATE:
			{
				State top = state_stack.top();
				state_stack.pop();
				if (TERM_STATE == state_stack.top()) {
					pathway->name = StrX(chars);
				} else {
					name = StrX(chars);
				}
				state_stack.push(top);
			}
			break;
		case DATABASE_SYMBOL_STATE:
			database_symbol = StrX(chars);
			break;
		case REFERENCE_STATE:
			reference = StrX(chars);
			break;
		default:
			break;
		}
	}

	void fatalError(const SAXParseException& exception)
	{
		char* message = XMLString::transcode(exception.getMessage());
		cout << "Fatal Error: " << message
			<< " at line: " << exception.getLineNumber()
			<< endl;
	}

	AmigoPathwayPtr pathway;
	stack<State> state_stack;
	string name;
	string database_symbol;
	string reference;
};


AmigoPathwayPtr
parse_pathway_xml(const std::string & xml_filename) {

	AmigoPathwayPtr pathway(new AmigoPathway);

	XMLPlatformUtils::Initialize();

	SAX2XMLReader* parser = XMLReaderFactory::createXMLReader();
	parser->setFeature(XMLUni::fgSAX2CoreValidation, true);   // optional
	parser->setFeature(XMLUni::fgSAX2CoreNameSpaces, true);   // optional

	MySAX2Handler* handler = new MySAX2Handler(pathway);
	parser->setContentHandler(handler);
	parser->setErrorHandler(handler);

	parser->parse(xml_filename.c_str());

	delete parser;
	delete handler;

	return pathway;
}


BIO_NS_END

