
#include <bio/svg_pssm_logo.h>
#include <bio/svg.h>
#include <bio/biobase_db.h>
#include <bio/biobase_match.h>
#include <bio/biobase_data_traits.h>
USING_BIO_NS

#include <boost/test/unit_test.hpp>
using namespace boost::unit_test;

XERCES_CPP_NAMESPACE_USE

#include <iostream>
#include <sstream>
using namespace std;


//#define VERBOSE_CHECKING

void
check_svg_pssm_logo()
{
	cout << "******* check_svg_pssm_logo()" << endl;

	XMLPlatformUtils::Initialize();

	DOMImplementation* impl
		= DOMImplementation::getImplementation();

	if (impl != NULL)
	{
		DOMDocumentType* doc_type
			= impl->createDocumentType(	XStr("svg"),
										NULL,
										XStr("svg-20000303-stylable.dtd") );
		XERCES_CPP_NAMESPACE::DOMDocument * doc = impl->createDocument(
					0,        // root element namespace URI.
					XStr("svg"),		// root element name
					doc_type);			// document type object (DTD).
		doc->setEncoding(XStr("UTF-8"));

		add_logo_defs(doc);

		TableLinkVec links;
		links.push_back(TableLink(MATRIX_DATA, 37));
		links.push_back(TableLink(MATRIX_DATA, 104));
		links.push_back(TableLink(MATRIX_DATA, 236));
		links.push_back(TableLink(MATRIX_DATA, 457));

		for (TableLinkVec::const_iterator i = links.begin();
			links.end() != i;
			++i)
		{
			const Matrix * matrix = BiobaseDb::singleton().get_matrices()[*i].get();
			const Pssm pssm = make_pssm(matrix);
			const seq_t sequence =
				(matrix->align_descs.begin() != matrix->align_descs.end())
					? matrix->align_descs.begin()->get()->sequence
					: "";
			DOMElement * pssm_logo = 
				create_svg_pssm_logo(
					pssm,
					doc,
					sequence);
			set_attribute(pssm_logo, "x", 0);
			set_attribute(pssm_logo, "y", BIO_NS::float_t(200 * (i - links.begin())));
			set_attribute(pssm_logo, "height", 200);
			set_attribute(pssm_logo, "width", BIO_NS::float_t(pssm.size() * 100));
			doc->getDocumentElement()->appendChild(pssm_logo);
		}

		dom_print(doc->getDocumentElement(), "logo_test.svg");
	}

	XMLPlatformUtils::Terminate();
}


void register_svg_tests(test_suite * test)
{
	test->add(BOOST_TEST_CASE(&check_svg_pssm_logo), 0);
}


