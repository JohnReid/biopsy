

#include <bio/defs.h>
USING_BIO_NS

#include <boost/test/unit_test.hpp>
#include <boost/test/parameterized_test.hpp>
#include <boost/assign/list_of.hpp>
#undef max
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/matrix_sparse.hpp>
using namespace boost;
using namespace boost::assign;
using boost::unit_test::test_suite;

#include <iostream>
using namespace std;

#include <math.h>

#include "blas1c.h"
#include "lapackc.h"
#include "arlsmat.h"
//#include "arcomp.h"
#include "arlsmat.h"
#include "arlnsmat.h"
#include "arlssym.h"
#include "arlgsym.h"
#include "arlsnsym.h"
#include "arlgnsym.h"
#include "arlscomp.h"
#include "arlgcomp.h"



#define VERBOSE_CHECKING



template <class FLOAT>
int AREig(arcomplex<FLOAT> EigVal[], int n, int nnz, arcomplex<FLOAT> A[],
          int irow[], int pcol[], int nev, char* which = "LM", int ncv = 0,
          FLOAT tol = 0.0, int maxit = 0, arcomplex<FLOAT>* resid = 0,
          bool AutoShift = true)
{

  // Creating a matrix in ARPACK++ format.

  ARluNonSymMatrix<arcomplex<FLOAT> > matrix(n, nnz, A, irow, pcol);

  // Defining the eigenvalue problem.

  ARluCompStdEig<FLOAT> prob(nev, matrix, which, ncv, tol,
                             maxit, resid, AutoShift);

  // Finding eigenvalues.

  return prob.Eigenvalues(EigVal);

} // complex standard problem, only eigenvalues, regular mode.


template <class FLOAT>
int AREig(arcomplex<FLOAT> EigVal[], arcomplex<FLOAT> EigVec[], int n,
          int nnz, arcomplex<FLOAT> A[], int irow[], int pcol[],
          int nev, char* which = "LM", int ncv = 0, FLOAT tol = 0.0,
          int maxit = 0, arcomplex<FLOAT>* resid = 0, bool AutoShift = true)
{

  // Creating a matrix in ARPACK++ format.

  ARluNonSymMatrix<arcomplex<FLOAT> > matrix(n, nnz, A, irow, pcol);

  // Defining the eigenvalue problem.

  ARluCompStdEig<FLOAT> prob(nev, matrix, which, ncv, tol,
                             maxit, resid, AutoShift);

  // Finding eigenvalues and eigenvectors.

  return prob.EigenValVectors(EigVec, EigVal);

} // complex standard problem, values and vectors, regular mode.


template <class FLOAT>
int AREig(arcomplex<FLOAT> EigVal[], int n, int nnz, arcomplex<FLOAT> A[],
          int irow[], int pcol[], arcomplex<FLOAT> sigma, int nev,
          char* which = "LM", int ncv = 0, FLOAT tol = 0.0, int maxit = 0,
          arcomplex<FLOAT>* resid = 0, bool AutoShift = true)
{

  // Creating a matrix in ARPACK++ format.

  ARluNonSymMatrix<arcomplex<FLOAT> > matrix(n, nnz, A, irow, pcol);

  // Defining the eigenvalue problem.

  ARluCompStdEig<FLOAT> prob(nev, matrix, sigma, which, ncv,
                             tol, maxit, resid, AutoShift);

  // Finding eigenvalues.

  return prob.Eigenvalues(EigVal);

} // complex standard problem, only eigenvalues, shift-and-invert.


template <class FLOAT>
int AREig(arcomplex<FLOAT> EigVal[], arcomplex<FLOAT> EigVec[], int n,
          int nnz, arcomplex<FLOAT> A[], int irow[], int pcol[],
          arcomplex<FLOAT> sigma, int nev, char* which = "LM",
          int ncv = 0, FLOAT tol = 0.0, int maxit = 0,
          arcomplex<FLOAT>* resid = 0, bool AutoShift = true)
{

  // Creating a matrix in ARPACK++ format.

  ARluNonSymMatrix<arcomplex<FLOAT> > matrix(n, nnz, A, irow, pcol);

  // Defining the eigenvalue problem.

  ARluCompStdEig<FLOAT> prob(nev, matrix, sigma, which, ncv,
                             tol, maxit, resid, AutoShift);

  // Finding eigenvalues and eigenvectors.

  return prob.EigenValVectors(EigVec, EigVal);

} // complex standard problem, values and vectors, shift-and-invert.


template <class FLOAT>
int AREig(arcomplex<FLOAT> EigVal[], int n, int nnzA,
          arcomplex<FLOAT> A[], int irowA[], int pcolA[], int nnzB,
          arcomplex<FLOAT> B[], int irowB[], int pcolB[], int nev,
          char* which = "LM", int ncv = 0, FLOAT tol = 0.0, int maxit = 0,
          arcomplex<FLOAT>* resid = 0, bool AutoShift = true)
{

  // Creating two matrices in ARPACK++ format.

  ARluNonSymMatrix<arcomplex<FLOAT> > matrixA(n, nnzA, A, irowA, pcolA);
  ARluNonSymMatrix<arcomplex<FLOAT> > matrixB(n, nnzB, B, irowB, pcolB);

  // Defining the eigenvalue problem.

  ARluCompGenEig<FLOAT> prob(nev, matrixA, matrixB, which,
                             ncv, tol, maxit, resid, AutoShift);

  // Finding eigenvalues.

  return prob.Eigenvalues(EigVal);

} // complex generalized problem, only eigenvalues, regular mode.


template <class FLOAT>
int AREig(arcomplex<FLOAT> EigVal[], arcomplex<FLOAT> EigVec[], int n,
          int nnzA, arcomplex<FLOAT> A[], int irowA[], int pcolA[],
          int nnzB, arcomplex<FLOAT> B[], int irowB[], int pcolB[],
          int nev, char* which = "LM", int ncv = 0, FLOAT tol = 0.0,
          int maxit = 0, arcomplex<FLOAT>* resid = 0, bool AutoShift = true)
{

  // Creating two matrices in ARPACK++ format.

  ARluNonSymMatrix<arcomplex<FLOAT> > matrixA(n, nnzA, A, irowA, pcolA);
  ARluNonSymMatrix<arcomplex<FLOAT> > matrixB(n, nnzB, B, irowB, pcolB);

  // Defining the eigenvalue problem.

  ARluCompGenEig<FLOAT> prob(nev, matrixA, matrixB, which,
                             ncv, tol, maxit, resid, AutoShift);

  // Finding eigenvalues and eigenvectors.

  return prob.EigenValVectors(EigVec, EigVal);

} // complex generalized problem, values and vectors, regular mode.


template <class FLOAT>
int AREig(arcomplex<FLOAT> EigVal[], int n, int nnzA, arcomplex<FLOAT> A[],
          int irowA[], int pcolA[], int nnzB, arcomplex<FLOAT> B[],
          int irowB[], int pcolB[], arcomplex<FLOAT> sigma, int nev,
          char* which = "LM", int ncv = 0, FLOAT tol = 0.0, int maxit = 0,
          arcomplex<FLOAT>* resid = 0, bool AutoShift = true)
{

  // Creating two matrices in ARPACK++ format.

  ARluNonSymMatrix<arcomplex<FLOAT> > matrixA(n, nnzA, A, irowA, pcolA);
  ARluNonSymMatrix<arcomplex<FLOAT> > matrixB(n, nnzB, B, irowB, pcolB);

  // Defining the eigenvalue problem.

  ARluCompGenEig<FLOAT> prob(nev, matrixA, matrixB, sigma, which,
                             ncv, tol, maxit, resid, AutoShift);

  // Finding eigenvalues.

  return prob.Eigenvalues(EigVal);

} // complex generalized problem, only eigenvalues, shift-and-invert mode.


template <class FLOAT>
int AREig(arcomplex<FLOAT> EigVal[], arcomplex<FLOAT> EigVec[], int n,
          int nnzA, arcomplex<FLOAT> A[], int irowA[], int pcolA[],
          int nnzB, arcomplex<FLOAT> B[], int irowB[], int pcolB[],
          arcomplex<FLOAT> sigma, int nev, char* which = "LM",
          int ncv = 0, FLOAT tol = 0.0, int maxit = 0,
          arcomplex<FLOAT>* resid = 0, bool AutoShift = true)
{

  // Creating two matrices in ARPACK++ format.

  ARluNonSymMatrix<arcomplex<FLOAT> > matrixA(n, nnzA, A, irowA, pcolA);
  ARluNonSymMatrix<arcomplex<FLOAT> > matrixB(n, nnzB, B, irowB, pcolB);

  // Defining the eigenvalue problem.

  ARluCompGenEig<FLOAT> prob(nev, matrixA, matrixB, sigma, which,
                             ncv, tol, maxit, resid, AutoShift);

  // Finding eigenvalues and eigenvectors.

  return prob.EigenValVectors(EigVec, EigVal);

} // complex generalized problem, values and vectors, shift-and-invert mode.


template <class FLOAT>
int AREig(double EigValR[], FLOAT EigValI[], int n, int nnz,
          FLOAT A[], int irow[], int pcol[], int nev, char* which = "LM",
          int ncv = 0, FLOAT tol = 0.0, int maxit = 0, FLOAT* resid = 0,
          bool AutoShift = true)
{

  // Creating a matrix in ARPACK++ format.

  ARluNonSymMatrix<FLOAT> matrix(n, nnz, A, irow, pcol);

  // Defining the eigenvalue problem.

  ARluNonSymStdEig<FLOAT> prob(nev, matrix, which, ncv, tol,
                               maxit, resid, AutoShift);

  // Finding eigenvalues.

  return prob.Eigenvalues(EigValR, EigValI);

} // real nonsymmetric standard problem, only eigenvalues, regular mode.


template <class FLOAT>
int AREig(float EigValR[], FLOAT EigValI[], int n, int nnz,
          FLOAT A[], int irow[], int pcol[], int nev, char* which = "LM",
          int ncv = 0, FLOAT tol = 0.0, int maxit = 0, FLOAT* resid = 0,
          bool AutoShift = true)
{

  // Creating a matrix in ARPACK++ format.

  ARluNonSymMatrix<FLOAT> matrix(n, nnz, A, irow, pcol);

  // Defining the eigenvalue problem.

  ARluNonSymStdEig<FLOAT> prob(nev, matrix, which, ncv, tol,
                               maxit, resid, AutoShift);

  // Finding eigenvalues.

  return prob.Eigenvalues(EigValR, EigValI);

} // real nonsymmetric standard problem, only eigenvalues, regular mode.


template <class FLOAT>
int AREig(FLOAT EigValR[], FLOAT EigValI[], FLOAT EigVec[], int n, int nnz,
          FLOAT A[], int irow[], int pcol[], int nev, char* which = "LM",
          int ncv = 0, FLOAT tol = 0.0, int maxit = 0, FLOAT* resid = 0,
          bool AutoShift = true)
{

  // Creating a matrix in ARPACK++ format.

  ARluNonSymMatrix<FLOAT> matrix(n, nnz, A, irow, pcol);

  // Defining the eigenvalue problem.

  ARluNonSymStdEig<FLOAT> prob(nev, matrix, which, ncv, tol,
                               maxit, resid, AutoShift);

  // Finding eigenvalues and eigenvectors.

  return prob.EigenValVectors(EigVec, EigValR, EigValI);

} // real nonsymmetric standard problem, values and vectors, regular mode.


template <class FLOAT>
int AREig(double EigValR[], FLOAT EigValI[], int n, int nnz,
          FLOAT A[], int irow[], int pcol[], FLOAT sigma, int nev,
          char* which = "LM", int ncv = 0, FLOAT tol = 0.0,
          int maxit = 0, FLOAT* resid = 0, bool AutoShift = true)
{

  // Creating a matrix in ARPACK++ format.

  ARluNonSymMatrix<FLOAT> matrix(n, nnz, A, irow, pcol);

  // Defining the eigenvalue problem.

  ARluNonSymStdEig<FLOAT> prob(nev, matrix, sigma, which, ncv,
                               tol, maxit, resid, AutoShift);

  // Finding eigenvalues.

  return prob.Eigenvalues(EigValR, EigValI);

} // real nonsymmetric standard problem, only eigenvalues, shift-and-invert.


template <class FLOAT>
int AREig(float EigValR[], FLOAT EigValI[], int n, int nnz,
          FLOAT A[], int irow[], int pcol[], FLOAT sigma, int nev,
          char* which = "LM", int ncv = 0, FLOAT tol = 0.0,
          int maxit = 0, FLOAT* resid = 0, bool AutoShift = true)
{

  // Creating a matrix in ARPACK++ format.

  ARluNonSymMatrix<FLOAT> matrix(n, nnz, A, irow, pcol);

  // Defining the eigenvalue problem.

  ARluNonSymStdEig<FLOAT> prob(nev, matrix, sigma, which, ncv,
                               tol, maxit, resid, AutoShift);

  // Finding eigenvalues.

  return prob.Eigenvalues(EigValR, EigValI);

} // real nonsymmetric standard problem, only eigenvalues, shift-and-invert.


template <class FLOAT>
int AREig(FLOAT EigValR[], FLOAT EigValI[], FLOAT EigVec[], int n, int nnz,
          FLOAT A[], int irow[], int pcol[], FLOAT sigma, int nev,
          char* which = "LM", int ncv = 0, FLOAT tol = 0.0, int maxit = 0,
          FLOAT* resid = 0, bool AutoShift = true)
{

  // Creating a matrix in ARPACK++ format.

  ARluNonSymMatrix<FLOAT> matrix(n, nnz, A, irow, pcol);

  // Defining the eigenvalue problem.

  ARluNonSymStdEig<FLOAT> prob(nev, matrix, sigma, which, ncv,
                               tol, maxit, resid, AutoShift);

  // Finding eigenvalues and eigenvectors.

  return prob.EigenValVectors(EigVec, EigValR, EigValI);

} // real nonsymmetric standard problem, values and vectors, shift-and-invert.


template <class FLOAT>
int AREig(double EigValR[], FLOAT EigValI[], int n, int nnzA,
          FLOAT A[], int irowA[], int pcolA[], int nnzB,
          FLOAT B[], int irowB[], int pcolB[], int nev,
          char* which = "LM", int ncv = 0, FLOAT tol = 0.0,
          int maxit = 0, FLOAT* resid = 0, bool AutoShift = true)
{

  // Creating two matrices in ARPACK++ format.

  ARluNonSymMatrix<FLOAT> matrixA(n, nnzA, A, irowA, pcolA);
  ARluNonSymMatrix<FLOAT> matrixB(n, nnzB, B, irowB, pcolB);

  // Defining the eigenvalue problem.

  ARluNonSymGenEig<FLOAT> prob(nev, matrixA, matrixB, which,
                               ncv, tol, maxit, resid, AutoShift);

  // Finding eigenvalues.

  return prob.Eigenvalues(EigValR, EigValI);

} // real nonsymmetric generalized problem, only eigenvalues, regular mode.


template <class FLOAT>
int AREig(float EigValR[], FLOAT EigValI[], int n, int nnzA,
          FLOAT A[], int irowA[], int pcolA[], int nnzB,
          FLOAT B[], int irowB[], int pcolB[], int nev,
          char* which = "LM", int ncv = 0, FLOAT tol = 0.0,
          int maxit = 0, FLOAT* resid = 0, bool AutoShift = true)
{

  // Creating two matrices in ARPACK++ format.

  ARluNonSymMatrix<FLOAT> matrixA(n, nnzA, A, irowA, pcolA);
  ARluNonSymMatrix<FLOAT> matrixB(n, nnzB, B, irowB, pcolB);

  // Defining the eigenvalue problem.

  ARluNonSymGenEig<FLOAT> prob(nev, matrixA, matrixB, which,
                               ncv, tol, maxit, resid, AutoShift);

  // Finding eigenvalues.

  return prob.Eigenvalues(EigValR, EigValI);

} // real nonsymmetric generalized problem, only eigenvalues, regular mode.


template <class FLOAT>
int AREig(FLOAT EigValR[], FLOAT EigValI[], FLOAT EigVec[], int n,
          int nnzA, FLOAT A[], int irowA[], int pcolA[],
          int nnzB, FLOAT B[], int irowB[], int pcolB[],
          int nev, char* which = "LM", int ncv = 0, FLOAT tol = 0.0,
          int maxit = 0, FLOAT* resid = 0, bool AutoShift = true)
{

  // Creating two matrices in ARPACK++ format.

  ARluNonSymMatrix<FLOAT> matrixA(n, nnzA, A, irowA, pcolA);
  ARluNonSymMatrix<FLOAT> matrixB(n, nnzB, B, irowB, pcolB);

  // Defining the eigenvalue problem.

  ARluNonSymGenEig<FLOAT> prob(nev, matrixA, matrixB, which,
                               ncv, tol, maxit, resid, AutoShift);

  // Finding eigenvalues and eigenvectors.

  return prob.EigenValVectors(EigVec, EigValR, EigValI);

} // real nonsymmetric generalized problem, values and vectors, regular mode.


template <class FLOAT>
int AREig(double EigValR[], FLOAT EigValI[], int n, int nnzA,
          FLOAT A[], int irowA[], int pcolA[], int nnzB,
          FLOAT B[], int irowB[], int pcolB[], FLOAT sigma,
          int nev, char* which = "LM", int ncv = 0, FLOAT tol = 0.0,
          int maxit = 0, FLOAT* resid = 0, bool AutoShift = true)
{

  // Creating two matrices in ARPACK++ format.

  ARluNonSymMatrix<FLOAT> matrixA(n, nnzA, A, irowA, pcolA);
  ARluNonSymMatrix<FLOAT> matrixB(n, nnzB, B, irowB, pcolB);

  // Defining the eigenvalue problem.

  ARluNonSymGenEig<FLOAT> prob(nev, matrixA, matrixB, sigma, which,
                               ncv, tol, maxit, resid, AutoShift);

  // Finding eigenvalues.

  return prob.Eigenvalues(EigValR, EigValI);

} // real nonsymmetric generalized problem, only eigenvalues,
  // real shift-and-invert mode.


template <class FLOAT>
int AREig(float EigValR[], FLOAT EigValI[], int n, int nnzA,
          FLOAT A[], int irowA[], int pcolA[], int nnzB,
          FLOAT B[], int irowB[], int pcolB[], FLOAT sigma,
          int nev, char* which = "LM", int ncv = 0, FLOAT tol = 0.0,
          int maxit = 0, FLOAT* resid = 0, bool AutoShift = true)
{

  // Creating two matrices in ARPACK++ format.

  ARluNonSymMatrix<FLOAT> matrixA(n, nnzA, A, irowA, pcolA);
  ARluNonSymMatrix<FLOAT> matrixB(n, nnzB, B, irowB, pcolB);

  // Defining the eigenvalue problem.

  ARluNonSymGenEig<FLOAT> prob(nev, matrixA, matrixB, sigma, which,
                               ncv, tol, maxit, resid, AutoShift);

  // Finding eigenvalues.

  return prob.Eigenvalues(EigValR, EigValI);

} // real nonsymmetric generalized problem, only eigenvalues,
  // real shift-and-invert mode.


template <class FLOAT>
int AREig(FLOAT EigValR[], FLOAT EigValI[], FLOAT EigVec[], int n,
          int nnzA, FLOAT A[], int irowA[], int pcolA[], int nnzB,
          FLOAT B[], int irowB[], int pcolB[], FLOAT sigma, int nev,
          char* which = "LM", int ncv = 0, FLOAT tol = 0.0,
          int maxit = 0, FLOAT* resid = 0, bool AutoShift = true)
{

  // Creating two matrices in ARPACK++ format.

  ARluNonSymMatrix<FLOAT> matrixA(n, nnzA, A, irowA, pcolA);
  ARluNonSymMatrix<FLOAT> matrixB(n, nnzB, B, irowB, pcolB);

  // Defining the eigenvalue problem.

  ARluNonSymGenEig<FLOAT> prob(nev, matrixA, matrixB, sigma, which,
                               ncv, tol, maxit, resid, AutoShift);

  // Finding eigenvalues and eigenvectors.

  return prob.EigenValVectors(EigVec, EigValR, EigValI);

} // real nonsymmetric generalized problem, values and vectors,
  // real shift-and-invert mode.


template <class FLOAT>
int AREig(FLOAT EigValR[], FLOAT EigValI[], int n, int nnzA, FLOAT A[],
          int irowA[], int pcolA[], int nnzB, FLOAT B[], int irowB[],
          int pcolB[], char part, FLOAT sigmaR, FLOAT sigmaI,
          int nev, char* which = "LM", int ncv = 0, FLOAT tol = 0.0,
          int maxit = 0, FLOAT* resid = 0, bool AutoShift = true)
{

  // Creating two matrices in ARPACK++ format.

  ARluNonSymMatrix<FLOAT> matrixA(n, nnzA, A, irowA, pcolA);
  ARluNonSymMatrix<FLOAT> matrixB(n, nnzB, B, irowB, pcolB);

  // Defining the eigenvalue problem.

  ARluNonSymGenEig<FLOAT> prob(nev, matrixA, matrixB, part,
                               sigmaR, sigmaI, which, ncv,
                               tol, maxit, resid, AutoShift);

  // Finding eigenvalues.

  return prob.Eigenvalues(EigValR, EigValI);

} // real nonsymmetric generalized problem, only eigenvalues,
  // complex shift-and-invert mode.


template <class FLOAT>
int AREig(FLOAT EigValR[], FLOAT EigValI[], FLOAT EigVec[], int n, int nnzA,
          FLOAT A[], int irowA[], int pcolA[], int nnzB, FLOAT B[],
          int irowB[], int pcolB[], char part, FLOAT sigmaR, FLOAT sigmaI,
          int nev, char* which = "LM", int ncv = 0, FLOAT tol = 0.0,
          int maxit = 0, FLOAT* resid = 0, bool AutoShift = true)
{

  // Creating two matrices in ARPACK++ format.

  ARluNonSymMatrix<FLOAT> matrixA(n, nnzA, A, irowA, pcolA);
  ARluNonSymMatrix<FLOAT> matrixB(n, nnzB, B, irowB, pcolB);

  // Defining the eigenvalue problem.

  ARluNonSymGenEig<FLOAT> prob(nev, matrixA, matrixB, part,
                               sigmaR, sigmaI, which, ncv,
                               tol, maxit, resid, AutoShift);

  // Finding eigenvalues and eigenvectors.

  return prob.EigenValVectors(EigVec, EigValR, EigValI);

} // real nonsymmetric generalized problem, values and vectors,
  // complex shift-and-invert mode.


template <class FLOAT>
int AREig(FLOAT EigVal[], int n, int nnz, FLOAT A[], int irow[],
          int pcol[], char uplo, int nev, char* which = "LM",
          int ncv = 0, FLOAT tol = 0.0, int maxit = 0,
          FLOAT* resid = 0, bool AutoShift = true)
{

  // Creating a matrix in ARPACK++ format.

  ARluSymMatrix<FLOAT> matrix(n, nnz, A, irow, pcol, uplo);

  // Defining the eigenvalue problem.

  ARluSymStdEig<FLOAT> prob(nev, matrix, which, ncv, tol,
                            maxit, resid, AutoShift);

  // Finding eigenvalues.

  return prob.Eigenvalues(EigVal);

} // real symmetric standard problem, only eigenvalues, regular mode.


template <class FLOAT>
int AREig(FLOAT EigVal[], FLOAT EigVec[], int n, int nnz, FLOAT A[],
          int irow[], int pcol[], char uplo, int nev,
          char* which = "LM", int ncv = 0, FLOAT tol = 0.0,
          int maxit = 0, FLOAT* resid = 0, bool AutoShift = true)
{

  // Creating a matrix in ARPACK++ format.

  ARluSymMatrix<FLOAT> matrix(n, nnz, A, irow, pcol, uplo);

  // Defining the eigenvalue problem.

  ARluSymStdEig<FLOAT> prob(nev, matrix, which, ncv, tol,
                            maxit, resid, AutoShift);

  // Finding eigenvalues and eigenvectors.

  return prob.EigenValVectors(EigVec, EigVal);

} // real symmetric standard problem, values and vectors, regular mode.


template <class FLOAT>
int AREig(FLOAT EigVal[], int n, int nnz, FLOAT A[], int irow[],
          int pcol[], char uplo, FLOAT sigma, int nev,
          char* which = "LM", int ncv = 0, FLOAT tol = 0.0,
          int maxit = 0, FLOAT* resid = 0, bool AutoShift = true)
{

  // Creating a matrix in ARPACK++ format.

  ARluSymMatrix<FLOAT> matrix(n, nnz, A, irow, pcol, uplo);

  // Defining the eigenvalue problem.

  ARluSymStdEig<FLOAT> prob(nev, matrix, sigma, which, ncv,
                            tol, maxit, resid, AutoShift);

  // Finding eigenvalues.

  return prob.Eigenvalues(EigVal);

} // real symmetric standard problem, only eigenvalues, shift-and-invert.


template <class FLOAT>
int AREig(FLOAT EigVal[], FLOAT EigVec[], int n, int nnz, FLOAT A[],
          int irow[], int pcol[], char uplo, FLOAT sigma,
          int nev, char* which = "LM", int ncv = 0, FLOAT tol = 0.0,
          int maxit = 0, FLOAT* resid = 0, bool AutoShift = true)
{

  // Creating a matrix in ARPACK++ format.

  ARluSymMatrix<FLOAT> matrix(n, nnz, A, irow, pcol, uplo);

  // Defining the eigenvalue problem.

  ARluSymStdEig<FLOAT> prob(nev, matrix, sigma, which, ncv,
                            tol, maxit, resid, AutoShift);

  // Finding eigenvalues and eigenvectors.

  return prob.EigenValVectors(EigVec, EigVal);

} // real symmetric standard problem, values and vectors, shift-and-invert.


template <class FLOAT>
int AREig(FLOAT EigVal[], int n, int nnzA, FLOAT A[], int irowA[],
          int pcolA[], int nnzB, FLOAT B[], int irowB[], int pcolB[],
          char uplo, int nev, char* which = "LM", int ncv = 0,
          FLOAT tol = 0.0, int maxit = 0, FLOAT* resid = 0,
          bool AutoShift = true)
{

  // Creating two matrices in ARPACK++ format.

  ARluSymMatrix<FLOAT> matrixA(n, nnzA, A, irowA, pcolA, uplo);
  ARluSymMatrix<FLOAT> matrixB(n, nnzB, B, irowB, pcolB, uplo);

  // Defining the eigenvalue problem.

  ARluSymGenEig<FLOAT> prob(nev, matrixA, matrixB, which,
                            ncv, tol, maxit, resid, AutoShift);

  // Finding eigenvalues.

  return prob.Eigenvalues(EigVal);

} // real symmetric generalized problem, only eigenvalues, regular mode.


template <class FLOAT>
int AREig(FLOAT EigVal[], FLOAT EigVec[], int n, int nnzA, FLOAT A[],
          int irowA[], int pcolA[], int nnzB, FLOAT B[], int irowB[],
          int pcolB[], char uplo, int nev, char* which = "LM",
          int ncv = 0, FLOAT tol = 0.0, int maxit = 0,
          FLOAT* resid = 0, bool AutoShift = true)
{

  // Creating two matrices in ARPACK++ format.

  ARluSymMatrix<FLOAT> matrixA(n, nnzA, A, irowA, pcolA, uplo);
  ARluSymMatrix<FLOAT> matrixB(n, nnzB, B, irowB, pcolB, uplo);

  // Defining the eigenvalue problem.

  ARluSymGenEig<FLOAT> prob(nev, matrixA, matrixB, which,
                            ncv, tol, maxit, resid, AutoShift);

  // Finding eigenvalues and eigenvectors.

  return prob.EigenValVectors(EigVec, EigVal);

} // real symmetric generalized problem, values and vectors, regular mode.


template <class FLOAT>
int AREig(FLOAT EigVal[], int n, int nnzA, FLOAT A[], int irowA[],
          int pcolA[], int nnzB, FLOAT B[], int irowB[], int pcolB[],
          char uplo, char InvertMode, FLOAT sigma, int nev,
          char* which = "LM", int ncv = 0, FLOAT tol = 0.0,
          int maxit = 0, FLOAT* resid = 0, bool AutoShift = true)
{

  // Creating two matrices in ARPACK++ format.

  ARluSymMatrix<FLOAT> matrixA(n, nnzA, A, irowA, pcolA, uplo);
  ARluSymMatrix<FLOAT> matrixB(n, nnzB, B, irowB, pcolB, uplo);

  // Defining the eigenvalue problem.

  ARluSymGenEig<FLOAT> prob(InvertMode, nev, matrixA, matrixB, sigma,
                            which, ncv, tol, maxit, resid, AutoShift);

  // Finding eigenvalues.

  return prob.Eigenvalues(EigVal);

} // real symmetric generalized problem, only eigenvalues,
  // shift-and-invert, buckling and Cayley modes.


template <class FLOAT>
int AREig(FLOAT EigVal[], FLOAT EigVec[], int n, int nnzA, FLOAT A[],
          int irowA[], int pcolA[], int nnzB, FLOAT B[], int irowB[],
          int pcolB[], char uplo, char InvertMode, FLOAT sigma,
          int nev, char* which = "LM", int ncv = 0, FLOAT tol = 0.0,
          int maxit = 0, FLOAT* resid = 0, bool AutoShift = true)
{

  // Creating two matrices in ARPACK++ format.

  ARluSymMatrix<FLOAT> matrixA(n, nnzA, A, irowA, pcolA, uplo);
  ARluSymMatrix<FLOAT> matrixB(n, nnzB, B, irowB, pcolB, uplo);

  // Defining the eigenvalue problem.

  ARluSymGenEig<FLOAT> prob(InvertMode, nev, matrixA, matrixB, sigma,
                            which, ncv, tol, maxit, resid, AutoShift);

  // Finding eigenvalues and eigenvectors.

  return prob.EigenValVectors(EigVec, EigVal);

} // real symmetric generalized problem, values and vectors,
  // shift-and-invert, buckling and Cayley modes.



template<class FLOAT, class INT>
void SymmetricMatrixA(INT nx, INT& n, INT& nnz, FLOAT* &A, 
                      INT* &irow, INT* &pcol, char uplo = 'L')

{

  // Defining internal variables.

  INT    i, j;
  FLOAT  h2, df, dd;

  // Defining constants.

  h2  = 1.0/(FLOAT(nx+1)*FLOAT(nx+1));
  dd  = 4.0/h2;
  df  = -1.0/h2;

  // Defining the number of columns and nonzero elements of matrix.

  n   = nx*nx;
  nnz = 3*n-2*nx;

  // Creating output vectors.

  A    = new FLOAT[nnz];
  irow = new INT[nnz];
  pcol = new INT[n+1];

  // Defining  matrix A.

  pcol[0] = 0;
  i       = 0;

  if (uplo == 'U') {

    for (j = 0; j < n; j++) {
      if (j >= nx) {
        A[i] = df;   irow[i++] = j-nx;
      }
      if ((j%nx) != 0) {
        A[i] = df;   irow[i++] = j-1;
      }
      A[i] = dd;     irow[i++] = j;
      pcol[j+1] = i;
    }

  }
  else {

    for (j = 0; j < n; j++) {
      A[i] = dd;     irow[i++] = j;
      if (((j+1)%nx) != 0) {
        A[i] = df;   irow[i++] = j+1;
      }
      if (j < n-nx) {
        A[i] = df;   irow[i++] = j+nx;
      }
      pcol[j+1] = i;
    }

  }

} // SymmetricMatrixA.


template<class FLOAT, class INT>
void Solution(INT nconv, INT n, INT nnz, FLOAT A[], INT irow[], INT pcol[],
              char uplo, FLOAT EigVal[], FLOAT* EigVec = 0)
/*
  Prints eigenvalues and eigenvectors of symmetric eigen-problems
  on standard "cout" stream.
*/

{

  INT                  i;
  FLOAT*               Ax;
  FLOAT*               ResNorm;
  ARluSymMatrix<FLOAT> matrix(n, nnz, A, irow, pcol, uplo);

  cout << endl << endl << "Testing ARPACK++ function AREig" << endl;
  cout << "Real symmetric eigenvalue problem: A*x - lambda*x \n \n";

  cout << "Dimension of the system            : " << n     << endl;
  cout << "Number of 'converged' eigenvalues  : " << nconv << endl << endl;

  // Printing eigenvalues.

  cout << "Eigenvalues:" << endl;

  for (i=0; i<nconv; i++) {
    cout << "  lambda[" << (i+1) << "]: " << EigVal[i] << endl;
  }
  cout << endl;

  // Printing eigenvectors.

  if (EigVec != 0) {

    // Finding the residual norm || A*x - lambda*x ||
    // for the nconv accurately computed eigenvectors.

    Ax      = new FLOAT[n];
    ResNorm = new FLOAT[nconv+1];

    for (i=0; i<nconv; i++) {
      matrix.MultMv(&EigVec[i*n], Ax);
      axpy(n, -EigVal[i], &EigVec[i*n], 1, Ax, 1);
      ResNorm[i] = nrm2(n, Ax, 1)/fabs(EigVal[i]);
    }

    for (i=0; i<nconv; i++) {
      cout << "||A*x(" << (i+1) << ") - lambda(" << (i+1);
      cout << ")*x(" << (i+1) << ")||: " << ResNorm[i] << endl;
    }
    cout << endl;

    delete[] Ax;
    delete[] ResNorm;

  }

} // Solution.


template<class FLOAT, class INT>
void Solution(INT nconv, INT n, INT nnzA, FLOAT A[], INT irowA[],
              INT pcolA[], INT nnzB, FLOAT B[], INT irowB[], INT pcolB[],
              char uplo, FLOAT EigVal[], FLOAT* EigVec = 0)
/*
  Prints eigenvalues and eigenvectors of symmetric generalized
  eigen-problem on standard "cout" stream.
*/

{

  INT                  i;
  FLOAT                *Ax, *Bx;
  FLOAT                *ResNorm;
  ARluSymMatrix<FLOAT> matrixA(n, nnzA, A, irowA, pcolA, uplo);
  ARluSymMatrix<FLOAT> matrixB(n, nnzB, B, irowB, pcolB, uplo);

  cout << endl << endl << "Testing ARPACK++ function AREig" << endl;
  cout << "Real symmetric generalized eigenvalue problem: A*x - lambda*B*x";
  cout << endl << endl;

  cout << "Dimension of the system            : " << n     << endl;
  cout << "Number of 'converged' eigenvalues  : " << nconv << endl << endl;

  // Printing eigenvalues.

  cout << "Eigenvalues:" << endl;

  for (i=0; i<nconv; i++) {
    cout << "  lambda[" << (i+1) << "]: " << EigVal[i] << endl;
  }
  cout << endl;

  // Printing eigenvectors.

  if (EigVec != 0) {

    // Printing the residual norm || A*x - lambda*B*x ||
    // for the nconv accurately computed eigenvectors.

    Ax      = new FLOAT[n];
    Bx      = new FLOAT[n];
    ResNorm = new FLOAT[nconv+1];

    for (i=0; i<nconv; i++) {
      matrixA.MultMv(&EigVec[i*n], Ax);
      matrixB.MultMv(&EigVec[i*n], Bx);
      axpy(n, -EigVal[i], Bx, 1, Ax, 1);
      ResNorm[i] = nrm2(n, Ax, 1)/fabs(EigVal[i]);
    }

    for (i=0; i<nconv; i++) {
      cout << "||A*x(" << i << ") - lambda(" << i;
      cout << ")*B*x(" << i << ")||: " << ResNorm[i] << endl;
    }
    cout << endl;

    delete[] Ax;
    delete[] Bx;
    delete[] ResNorm;

  }

} // Solution.



void
check_eigen_solve()
{
	cout << "******* check_eigen_solve()" << endl;

	int     nx;
	int     n;           // Dimension of the problem.
	int     nconv;       // Number of "converged" eigenvalues.
	int     nnz;         // Number of nonzero elements in A.
	int*    irow;        // pointer to an array that stores the row
	// indices of the nonzeros in A.
	int*    pcol;        // pointer to an array of pointers to the
	// beginning of each column of A in vector A.
	double* A;           // pointer to an array that stores the
	// nonzero elements of A.
	double EigVal[101];  // Eigenvalues.
	double EigVec[1001]; // Eigenvectors stored sequentially.
	char    uplo;        // Variable that indicates whether the upper
	// (uplo='U') ot the lower (uplo='L') part of
	// A will be stored in A, irow and pcol.

	// Creating a 100x100 matrix.

	nx = 10;
	uplo = 'U';
	SymmetricMatrixA(nx, n, nnz, A, irow, pcol, uplo);

	// Finding the four eigenvalues with smallest magnitude and
	// the related eigenvectors.

	nconv = AREig(EigVal, EigVec, n, nnz, A, irow, pcol, uplo, 4, "SM");

	// Printing solution.

	Solution(nconv, n, nnz, A, irow, pcol, uplo, EigVal, EigVec);
}


void
register_eigen_solve_tests( boost::unit_test::test_suite * test )
{
	test->add( BOOST_TEST_CASE( &check_eigen_solve ), 0);
}
