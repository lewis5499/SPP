/*------------------------------------------------------------------------------
* matrix.h : SPP software matrix operations
*
*          Copyright (C) 2023 by H.Z. Liu, All rights reserved.
*
* options : none
*
* references :  [1]"RTK_Structs.h"
*
* version : $Revision: 1.1 $ $Date: 2023/10/2 11:53:00 $
*
* history : 2023/10/2 1.0 new
*
*-----------------------------------------------------------------------------*/
#ifndef _MATRIX_H_
#define _MATRIX_H_

#include <string>

using namespace std;

namespace matrix {
	/* find the biggest element */
	double maxElement(const double *a, const short row, const short col);
	/* find the smallest element */
	double minElement(const double *a, const short row, const short col);
	/* inner product of vectors */
	double dot(const double *a, const double *b, int n);
	/* euclid norm of vector */
	double norm(const double *a, int n);
	/* outer product of 3d vectors */
	void cross3(const double *a, const double *b, double *c);
	/* normalize 3d vector */
	int normv3(const double *a, double *b);
	/* matrix shallow copy */
	void matcpy(double *A, const double *B, int n, int m);
	/* deep copy a vec/matrix */
	void matdpcpy(double *dst, const double *ori, int length);
	/* perform scalar*mat */
	void scalar_matMul(const double scalar, const double *oriMat, double *dstMat, int n);
	/* new matrix */
	double *mat(int n, int m);
	/* new integer matrix */
	int *imat(int n, int m);
	/* zero matrix */
	double *zeros(int n, int m);
	/* identity matrix */
	double *eye(int n);
	void eye(double* mat, int height);

	/* Display of matrix: String */
	string matToStr(const double* iMatrix, const unsigned short row, const unsigned short column);		
	string imatToStr(const int* iMatrix, const unsigned short row, const unsigned short column);
	/* Transpose of matrix */
	bool matTran(const double* iMatrix, const unsigned short row, const unsigned short column, double* iNewMatrix);	
	/* Inverse of matrix */
	bool matInv(const double* iMatrix, const unsigned short size, double* iMatrixInv);
	/* inverse of matrix */
	int matinv(double *A, int n);
	/* Matrix addition */
	bool matAdd(const double* MatrixLin, const double* MatrixRin, double *MatrixResult,		
				const unsigned short row, const unsigned short column);
	/* Matrix subtraction */
	bool matSub(const double* MatrixLin, const double* MatrixRin, double* MatrixResult,		
				const unsigned short row, const unsigned short column);
	/*Matrix multiplication*/
	bool matMul(const double* MatrixLin, const unsigned short rowinL, const unsigned short colinL,		
				const double* MatrixRin, const unsigned short rowinR, const unsigned short colinR, double *MatrixResult);
	void matmul(const char *tr, int n, int k, int m, double alpha,
				const double *A, const double *B, double beta, double *C);
	/* LU: Determinant of a square matrix */
	double mat_nn_LU(const double* iMatrix, const unsigned short size);	
	/* Laplace: Determinant of a square matrix */
	double mat_nn_LA(const double* iMatrix, const unsigned short size);	
	/* Gaussian: Determinant of a square matrix */
	double mat_nn_Gauss(const double* iMatrix, const unsigned short size);	

	/* trim a str */
	static string trim(string s);
	/* LU decomposition */
	static int ludcmp(double *A, int n, int *indx, double *d);
	/* LU back-substitution */
	static void lubksb(const double *A, int n, const int *indx, double *b);
}
#endif