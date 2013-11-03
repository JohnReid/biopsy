#ifndef BIO_MATRIX_TEST_DATA
#define BIO_MATRIX_TEST_DATA

#include "bio/defs.h"
#include "bio/matrix.h"
#include "bio/test_data.h"
#include "bio/binding_model.h"


BIO_NS_START

struct MatrixTestData : public BiFaTestData
{
protected:
	/** The alignment used to generate the test sequence. */
	AlignDescPtr binding_site;

public:
	/** Generate some test data from a matrix by the following: Take a random sequence of the given length
	and replace a section of it with one of the sequences used to generate the matrix. */
	MatrixTestData(
		const Matrix * matrix,
		size_t seq_length);

	virtual ~MatrixTestData();

	/** The matrix that binds in the sequence. */
	const Matrix * matrix;

	Test::ptr_t get_test_for_algorithm( BiFaAlgorithm::ptr_t algorithm );
};


struct MatrixTest
	: Test
{
	MatrixTestData & data;
	BiFaAlgorithm::ptr_t algorithm;

	MatrixTest(
		MatrixTestData & data,
		BiFaAlgorithm::ptr_t algorithm );

	virtual ~MatrixTest();

	virtual BinaryTestResults get_results( double threshold );
};



BIO_NS_END

#endif //BIO_MATRIX_TEST_DATA
