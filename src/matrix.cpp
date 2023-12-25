#include "matrix.h"
#include "consts.h"
#include <cmath>
#include <algorithm>

namespace matrix {

	/* max element ----------------------------------------------------------------
	* max element of vectors/matrix
	* args   : double *a        I   matrix a  
	*          short  row       I   num of rows of a
	*		   short  col		I   num of cols of a
	* return : max element of a
	*-----------------------------------------------------------------------------*/
	double maxElement(const double *a, const short row, const short col) {
		return *max_element(a, a + row * col);
	}
	
	/* min element ----------------------------------------------------------------
	* min element of vectors/matrix
	* args   : double *a        I   matrix a
	*          short  row       I   num of rows of matrix a
	*		   short  col		I   num of cols of matrix a
	* return : min element of matrix a
	*-----------------------------------------------------------------------------*/
	double minElement(const double *a, const short row, const short col) {
		return *min_element(a, a + row * col);
	}

	/* inner product ---------------------------------------------------------------
	* inner product of vectors
	* args   : double *a,*b     I   vector a,b (n x 1)
	*          int    n         I   size of vector a,b
	* return : a'*b
	*-----------------------------------------------------------------------------*/
	double dot(const double *a, const double *b, int n) {
		double c = 0.0;

		while (--n >= 0) c += a[n] * b[n];
		return c;
	}

	/* euclid norm -----------------------------------------------------------------
	* euclid norm of vector
	* args   : double *a        I   vector a (n x 1)
	*          int    n         I   size of vector a
	* return : || a ||
	*-----------------------------------------------------------------------------*/
	double norm(const double *a, int n)	{
		return sqrt(dot(a, a, n));
	}

	/* conmmon function for trimming a str */
	static string trim(string s) {
		if (s.empty()) {
			return s;
		}
		s.erase(0, s.find_first_not_of(" "));
		s.erase(s.find_last_not_of(" ") + 1);
		return s;
	}

	/* outer product of 3d vectors -------------------------------------------------
	* outer product of 3d vectors
	* args   : double *a,*b     I   vector a,b (3 x 1)
	*          double *c        O   outer product (a x b) (3 x 1)
	* return : none
	*-----------------------------------------------------------------------------*/
	void cross3(const double *a, const double *b, double *c) {
		c[0] = a[1] * b[2] - a[2] * b[1];
		c[1] = a[2] * b[0] - a[0] * b[2];
		c[2] = a[0] * b[1] - a[1] * b[0];
	}

	/* normalize 3d vector ---------------------------------------------------------
	* normalize 3d vector
	* args   : double *a        I   vector a (3 x 1)
	*          double *b        O   normlized vector (3 x 1) || b || = 1
	* return : status (1:ok,0:error)
	*-----------------------------------------------------------------------------*/
	int normv3(const double *a, double *b) {
		double r;
		if ((r = norm(a, 3)) <= 0.0) return 0;
		b[0] = a[0] / r;
		b[1] = a[1] / r;
		b[2] = a[2] / r;
		return 1;
	}

	/* copy matrix -----------------------------------------------------------------
	* copy matrix
	* args   : double *A        O   destination matrix A (n x m)
	*          double *B        I   source matrix B (n x m)
	*          int    n,m       I   number of rows and columns of matrix
	* return : none
	*-----------------------------------------------------------------------------*/
	void matcpy(double *A, const double *B, int n, int m) {
		memcpy(A, B, sizeof(double)*n*m);
	}

	/* deep copy nd vector ---------------------------------------------------------
	* deep copy nd vector
	* args   : int n			I	dimensions of vec
	*          double *ori      I   vector (n x 1)
	*          double *dst      O   vector (n x 1) 
	* return : none
	*-----------------------------------------------------------------------------*/
	void matdpcpy(double *dst, const double *ori, int n) {
		while (--n >= 0) {
			dst[n] = ori[n];
		}
	}

	/* deep copy nd vector ---------------------------------------------------------
	* deep copy nd vector
	* args   : int n			I	dimensions of vec
	*		   int scalar		I	the scalar to multiply
	*          double *ori      I   vector (n x 1)
	*          double *dst      O   vector (n x 1)
	* return : none
	*-----------------------------------------------------------------------------*/
	void scalar_matMul(const double scalar, const double *ori, double *dst, int n) {
		while (--n >= 0) {
			dst[n] = ori[n] * scalar;
		}
	}

	/* new matrix ------------------------------------------------------------------
	* allocate memory of matrix
	* args   : int    n,m       I   number of rows and columns of matrix
	* return : matrix pointer (if n<=0 or m<=0, return NULL)
	*-----------------------------------------------------------------------------*/
	extern double *mat(int n, int m) {
		double *p;

		if (n <= 0 || m <= 0) return NULL;
		if (!(p = (double *)malloc(sizeof(double)*n*m))) {
			printf("matrix memory allocation error: n=%d,m=%d\n", n, m);
		}
		return p;
	}
	/* new integer matrix ----------------------------------------------------------
	* allocate memory of integer matrix
	* args   : int    n,m       I   number of rows and columns of matrix
	* return : matrix pointer (if n<=0 or m<=0, return NULL)
	*-----------------------------------------------------------------------------*/
	extern int *imat(int n, int m) {
		int *p;

		if (n <= 0 || m <= 0) return NULL;
		if (!(p = (int *)malloc(sizeof(int)*n*m))) {
			printf("integer matrix memory allocation error: n=%d,m=%d\n", n, m);
		}
		return p;
	}

	/* zero matrix -----------------------------------------------------------------
	* generate new zero matrix
	* args   : int    n,m       I   number of rows and columns of matrix
	* return : matrix pointer (if n<=0 or m<=0, return NULL)
	*-----------------------------------------------------------------------------*/
	double *zeros(int n, int m)	{
		double *p;

	#if NOCALLOC 
		if ((p = mat(n, m))) for (n = n * m - 1; n >= 0; n--) p[n] = 0.0;
	#else
		if (n <= 0 || m <= 0) return NULL;
		if (!(p = (double *)calloc(sizeof(double), n*m))) {
			printf("matrix memory allocation error: n=%d,m=%d\n", n, m);
		}
	#endif
		return p;
	}

	/* identity matrix -------------------------------------------------------------
	* generate new identity matrix
	* args   : int    n         I   number of rows and columns of matrix
	* return : matrix pointer (if n<=0, return NULL)
	*-----------------------------------------------------------------------------*/
	double *eye(int n) {
		double *p;
		int i;

		if ((p = zeros(n, n))) for (i = 0; i < n; i++) p[i + i * n] = 1.0;
		return p;
	}

	void eye(double *mat, int height)
	{
		int i, n;
		n = height * height;
		for (i = 0; i < n; i++)	{ mat[i] = 0.0;	}
		for (i = 0; i < height; i++) { mat[i + i * height] = 1.0; }
	}

	/* Display of matrix: String */
	string matToStr(const double* iMatrix, const unsigned short row, const unsigned short column) {
		string strmat = "";
		for (unsigned short i = 0; i < row; ++i) {
			for (unsigned short j = 0; j < column; ++j) {
				strmat += to_string(iMatrix[column*i + j]);
				strmat += ' ';
			}
			strmat = trim(strmat);
			strmat += "\r\n";
		}
		return strmat;
	}
	string imatToStr(const int* iMatrix, const unsigned short row, const unsigned short column) {
		string strmat = "";
		for (unsigned short i = 0; i < row; ++i) {
			for (unsigned short j = 0; j < column; ++j) {
				strmat += to_string(iMatrix[column*i + j]);
				strmat += ' ';
			}
			strmat = trim(strmat);
			strmat += "\r\n";
		}
		return strmat;
	}

	/* Transpose of matrix */
	bool matTran(const double* iMatrix, const unsigned short row, const unsigned short column, double* iNewMatrix) {
		for (unsigned short i = 0; i < column; ++i) {
			for (unsigned short j = 0; j < row; ++j) {
				iNewMatrix[i * row + j] = iMatrix[j * column + i];
			}
		}
		return true;
	}

	/* Inverse of matrix: Gaussian-Jordan elimination with full pivoting */
	bool matInv(const double* iMatrix, const unsigned short size, double* iMatrixInv) {
		// Copy the input matrix and get unitized inverse matrix
		double *imatrixTemp = new double[size * size];
		for (unsigned short i = 0; i < size * size; ++i) {
			imatrixTemp[i] = iMatrix[i];
			iMatrixInv[i] = (i % (size + 1) == 0) ? 1.0 : 0.0;
		}
		// save the consequences of row/column swaps
		unsigned short *indices = new unsigned short[size*2]; 
		// algorithm begins
		for (unsigned short i = 0; i < size; ++i) {
			// select pivot
			unsigned short pivotRow = i;
			unsigned short pivotCol = i;
			double maxVal = abs(imatrixTemp[i * size + i]);
			for (unsigned short row = i; row < size; ++row) {
				for (unsigned short col = i; col < size; ++col) {
					double absVal = abs(imatrixTemp[row * size + col]);
					if (absVal > maxVal) {
						maxVal = absVal;
						pivotRow = row;
						pivotCol = col;
					}
				}
			}
			// determine if the matrix is singular
			if (maxVal < EPSILON) {
				delete[] imatrixTemp;
				delete[] indices;
				printf("matInv:'iMatrix' RANK Deficiency!\n");
				return false; // a singular matrix cannot be inverted
			}
			// recording the order of row swaps
			indices[i] = pivotRow;
			if (pivotRow != i) {
				// swap rows
				for (unsigned short col = 0; col < size; ++col) {
					swap(imatrixTemp[i * size + col], imatrixTemp[pivotRow * size + col]);
					swap(iMatrixInv[i * size + col], iMatrixInv[pivotRow * size + col]);
				}
			}
			indices[size + i] = pivotCol;
			if (pivotCol != i) {
				// swap columns
				for (unsigned short row = 0; row < size; ++row) {
					swap(imatrixTemp[row * size + i], imatrixTemp[row * size + pivotCol]);
					swap(iMatrixInv[row * size + i], iMatrixInv[row * size + pivotCol]);
				}
			}
			// unitize pivot 
			double pivot = imatrixTemp[i * size + i];
			for (unsigned short col = 0; col < size; ++col) {
				imatrixTemp[i * size + col] /= pivot;
				iMatrixInv[i * size + col] /= pivot;
			}
			// eliminate non-pivot elements
			for (unsigned short row = 0; row < size; ++row) {
				if (row != i) {
					double factor = imatrixTemp[row * size + i];
					for (unsigned short col = 0; col < size; ++col) {
						imatrixTemp[row * size + col] -= factor * imatrixTemp[i * size + col];
						iMatrixInv[row * size + col] -= factor * iMatrixInv[i * size + col];
					}
				}
			}
		}
		// recover row/column swaps
		for (int i = size - 1; i >= 0; --i) {
			if (indices[i] != i) {
				// row
				for (unsigned short col = 0; col < size; ++col) {
					swap(iMatrixInv[i * size + col], iMatrixInv[indices[i] * size + col]);
				}
			}
			if (indices[size + i] != i) {
				// column
				for (unsigned short row = 0; row < size; ++row) {
					swap(iMatrixInv[row * size + i], iMatrixInv[row * size + indices[size + i]]);
				}
			}
		}
		delete[] imatrixTemp;
		delete[] indices;
		return true;
	}

	/* inverse of matrix ---------------------------------------------------------
	*  note: A is I&O, after this operation the value of matrix 'A' will get changed
	*---------------------------------------------------------------------------*/
	int matinv(double *A, int n) {
		double d, *B;
		int i, j, *indx;

		indx = imat(n, 1); B = mat(n, n); matcpy(B, A, n, n);
		if (ludcmp(B, n, indx, &d)) { free(indx); free(B); return -1; }
		for (j = 0; j < n; j++) {
			for (i = 0; i < n; i++) A[i + j * n] = 0.0;
			A[j + j * n] = 1.0;
			lubksb(B, n, indx, A + j * n);
		}
		free(indx); free(B);
		return 0;
	}
	
	/* Matrix addition */
	bool matAdd(const double* MatrixLin, const double* MatrixRin, double *MatrixResult,
		const unsigned short row, const unsigned short column) {
		for (unsigned short i = 0; i < row; ++i) {
			for (unsigned short j = 0; j < column; ++j) {
				MatrixResult[column*i + j] = MatrixLin[column*i + j] + MatrixRin[column*i + j];
			}
		}
		return true;
	}
	
	/* Matrix subtraction */
	bool matSub(const double* MatrixLin, const double* MatrixRin, double* MatrixResult,
		const unsigned short row, const unsigned short column) {
		for (unsigned short i = 0; i < row; ++i) {
			for (unsigned short j = 0; j < column; ++j) {
				MatrixResult[column*i + j] = MatrixLin[column*i + j] - MatrixRin[column*i + j];
			}
		}
		return true;
	}

	/*Matrix multiplication*/
	bool matMul(const double* MatrixLin, const unsigned short rowinL, const unsigned short colinL,
		const double* MatrixRin, const unsigned short rowinR, const unsigned short colinR, double *MatrixResult) {
		if (colinL != rowinR) {
			printf("matMul:The dimensions of 'MatrixLin' and 'MatrixRin' mismatch!\n");
			return false;
		}
		for (unsigned short i = 0; i < rowinL; ++i) {
			for (unsigned short j = 0; j < colinR; ++j) {
				MatrixResult[colinR*i + j] = 0.0;
				for (unsigned short k = 0; k < colinL; ++k) {	
					MatrixResult[colinR*i + j] += MatrixLin[colinL*i + k] * MatrixRin[colinR*k + j];
					double a = MatrixResult[colinR*i + j];
				}
			}
		}
		return true;
	}
	/* multiply matrix -----------------------------------------------------------*/
	void matmul(const char *tr, int n, int k, int m, double alpha,
		const double *A, const double *B, double beta, double *C) {
		double d;
		int i, j, x, f = tr[0] == 'N' ? (tr[1] == 'N' ? 1 : 2) : (tr[1] == 'N' ? 3 : 4);

		for (i = 0; i < n; i++) for (j = 0; j < k; j++) {
			d = 0.0;
			switch (f) {
			case 1: for (x = 0; x < m; x++) d += A[i + x * n] * B[x + j * m]; break;
			case 2: for (x = 0; x < m; x++) d += A[i + x * n] * B[j + x * k]; break;
			case 3: for (x = 0; x < m; x++) d += A[x + i * m] * B[x + j * m]; break;
			case 4: for (x = 0; x < m; x++) d += A[x + i * m] * B[j + x * k]; break;
			}
			if (beta == 0.0) C[i + j * n] = alpha * d; else C[i + j * n] = alpha * d + beta * C[i + j * n];
		}
	}

	/* LU: Determinant of a square matrix */
	double mat_nn_LU(const double* iMatrix, const unsigned short size) {
		double* LU = new double[size * size];
		for (int i = 0; i < size * size; ++i) {
			LU[i] = iMatrix[i];
		}

		double det = 1.0;
		for (int i = 0; i < size; ++i) {
			if (abs(LU[i * size + i]) < EPSILON) {
				delete[] LU;
				return 0.0;
			}
			det *= LU[i * size + i];
			for (int j = i + 1; j < size; ++j) {
				LU[j * size + i] /= LU[i * size + i];
				for (int k = i + 1; k < size; ++k) {
					LU[j * size + k] -= LU[j * size + i] * LU[i * size + k];
				}
			}
		}
		delete[] LU;
		return det;
	}
	
	
	/* Laplace: Determinant of a square matrix */
	double mat_nn_LA(const double* iMatrix, const unsigned short size) {
		if (size == 1) {
			return iMatrix[0];
		}
		else if (size == 2) {
			return iMatrix[0] * iMatrix[3] - iMatrix[1] * iMatrix[2];
		}
		else {
			double det = 0;
			for (unsigned short i = 0; i < size; ++i) {
				double* submatrix = new double[(size - 1) * (size - 1)];
				unsigned short sub_i = 0;
				for (unsigned short row = 1; row < size; ++row) {
					unsigned short sub_j = 0;
					for (unsigned short col = 0; col < size; ++col) {
						if (col != i) {
							submatrix[sub_i * (size - 1) + sub_j] = iMatrix[row * size + col];
							++sub_j;
						}
					}
					++sub_i;
				}
				double sub_det = mat_nn_LA(submatrix, size - 1);
				int sign = (i % 2 == 0) ? 1 : -1;
				det += sign * iMatrix[i] * sub_det;
				delete[] submatrix;
			}
			return det;
		}
	}

	/* Gaussian: Determinant of a square matrix */
	double mat_nn_Gauss(const double* iMatrix, const unsigned short size) {
		double* matrix = new double[size * size];
		for (int i = 0; i < size * size; i++) {
			matrix[i] = iMatrix[i];
		}

		double det = 1.0;
		for (int i = 0; i < size; i++) {
			int max = i;
			for (int j = i + 1; j < size; j++) {
				if (abs(matrix[j * size + i]) > abs(matrix[max * size + i])) {
					max = j;
				}
			}

			if (i != max) {
				for (int j = 0; j < size; j++) {
					swap(matrix[i * size + j], matrix[max * size + j]);
				}
				det *= -1;
			}

			det *= matrix[i * size + i];
			if (abs(det) < EPSILON) {
				delete[] matrix;
				return 0;
			}

			for (int j = i + 1; j < size; j++) {
				matrix[j * size + i] /= matrix[i * size + i];
				for (int k = i + 1; k < size; k++) {
					matrix[j * size + k] -= matrix[j * size + i] * matrix[i * size + k];
				}
			}
		}

		delete[] matrix;
		return det;
	}

	/* LU decomposition ----------------------------------------------------------*/
	static int ludcmp(double *A, int n, int *indx, double *d) {
		double big, s, tmp, *vv = mat(n, 1);
		int i, imax = 0, j, k;

		*d = 1.0;
		for (i = 0; i < n; i++) {
			big = 0.0; for (j = 0; j < n; j++) if ((tmp = fabs(A[i + j * n])) > big) big = tmp;
			if (big > 0.0) vv[i] = 1.0 / big; else { free(vv); return -1; }
		}
		for (j = 0; j < n; j++) {
			for (i = 0; i < j; i++) {
				s = A[i + j * n]; for (k = 0; k < i; k++) s -= A[i + k * n] * A[k + j * n]; A[i + j * n] = s;
			}
			big = 0.0;
			for (i = j; i < n; i++) {
				s = A[i + j * n]; for (k = 0; k < j; k++) s -= A[i + k * n] * A[k + j * n]; A[i + j * n] = s;
				if ((tmp = vv[i] * fabs(s)) >= big) { big = tmp; imax = i; }
			}
			if (j != imax) {
				for (k = 0; k < n; k++) {
					tmp = A[imax + k * n]; A[imax + k * n] = A[j + k * n]; A[j + k * n] = tmp;
				}
				*d = -(*d); vv[imax] = vv[j];
			}
			indx[j] = imax;
			if (A[j + j * n] == 0.0) { free(vv); return -1; }
			if (j != n - 1) {
				tmp = 1.0 / A[j + j * n]; for (i = j + 1; i < n; i++) A[i + j * n] *= tmp;
			}
		}
		free(vv);
		return 0;
	}

	/* LU back-substitution ------------------------------------------------------*/
	static void lubksb(const double *A, int n, const int *indx, double *b) {
		double s;
		int i, ii = -1, ip, j;

		for (i = 0; i < n; i++) {
			ip = indx[i]; s = b[ip]; b[ip] = b[i];
			if (ii >= 0) for (j = ii; j < i; j++) s -= A[i + j * n] * b[j]; else if (s) ii = i;
			b[i] = s;
		}
		for (i = n - 1; i >= 0; i--) {
			s = b[i]; for (j = i + 1; j < n; j++) s -= A[i + j * n] * b[j]; b[i] = s / A[i + i * n];
		}
	}

}