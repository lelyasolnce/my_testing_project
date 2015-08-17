#ifndef __LINBOX_matrix_SlicedPolynomialMatrix_SlicedPolynomialMatrix_H
#define __LINBOX_matrix_SlicedPolynomialMatrix_SlicedPolynomialMatrix_H

#include <vector>
#include <istream>
#include <ostream>
/*!#include "linbox/matrix/DenseMatrix/blas-matrix.h"*/
/**/
#include "/home/slb/fakeit/myBM.h"/*!*/
/*!#include <givaro/modular.h>*/
/**/
#include "/home/slb/fakeit/my_modular.h"
//!#include <givaro/givpoly1dense.h>
//!#include <givaro/givpoly1denseops.inl>

namespace LinBox
{
	template <class _Field, class _MatrixElement = double>
	class SlicedPolynomialMatrix
	/*Requirement: variable class _Field must be a representation of a field of order p^n
	 *and have functions that return p (characteristic()) and p^n(cardinality()).
	 *For example, Givaro::GfqDom.
	 */

        /*Description: this class implements representation of a matrix over GF
	 *as a vector of length n of matrices over F (BlasMatrices),
	 *which corresponds to representation of this matrix as a polynomial of degree e-1 with matrix-coefficients.
	 */
	{
	public:
		typedef _Field Field;
		typedef typename Field::Element Element;
		typedef _MatrixElement MatrixElement;
		typedef std::vector<MatrixElement> Rep;
		typedef Givaro::Modular<MatrixElement> IntField; //???
		typedef SlicedPolynomialMatrix<Field, MatrixElement> Self_t;
		//!typedef Givaro::Poly1Dom<IntField, Dense>::Rep polynomial;
	private:
		_Field *GF;
		IntField F;
	private:
		int e;//GF.cardinality() == p^e
		std::vector<BlasMatrix<IntField, Rep>> V;
	//!public:
		//!polynomial irreducible;
	private:
		//polynomial& modulo(polynomial& g, polynomial&h, polynomial& f);
		//bool rabin(int n, int q, int a, int b);
		//void setIrreduciblePlynomial(int max_steps = 1000000);

						////////////////
		        			//Constructors//
						////////////////

	public:
		/*! Allocates a vector of new \f$ 0 \times 0\f$ matrices (shaped and ready). Irreducible polynomial is chosen randomly.
		 */
		SlicedPolynomialMatrix (_Field &BF);

		/*Allocates a vector of $ m1 \times m2\f$ zero matrices (shaped and ready). Irreducible polynomial is chosen randomly.
		 */
		SlicedPolynomialMatrix (_Field &BF, const size_t & m1, const size_t &m2);

						///////////////
						// Destructor//
						///////////////

	public:
		~SlicedPolynomialMatrix() {}

						////////////////////////
	          				//dimensions of vector//
						////////////////////////

	public:
		/*Get length of V.
		 * @returns length of V
		 */
		size_t length() const;
		
		/*Get the number of rows in a matrix.
		 * @returns Number of rows in a matrix
		 */
		size_t rowdim() const;

		/* Get the number of columns in a matrix.
		 * @returns Number of columns in a matrix
		 */
		size_t coldim() const;
		
	                    			/////////////////
	                    			//return fields//
	                    			/////////////////
	
	public:
		Field& fieldGF();

		IntField& fieldF();

						/////////////////////////
		        			//functions for entries//
						/////////////////////////

	public:
		/* Set the entry of the m-th matrix-coefficient at the (i, j) position to a_mij.
		 * @param m matrix-coefficient number, 0...length() - 1
		 * @param i Row number, 0...rowdim () - 1
		 * @param j Column number 0...coldim () - 1
		 * @param a_mij Element to set
		 */
		void setEntry (size_t m, size_t i, size_t j, _MatrixElement &a_mij);
		
	private:
		/* Get a writeable reference to the m-th matrix-coefficient at the (i, j) position.
		 * @param m matrix-coefficient number, 0...length() - 1
		 * @param i Row index of entry
		 * @param j Column index of entry
		 * @returns Reference to matrix entry
		 */
		_MatrixElement &refEntry (size_t m, size_t i, size_t j);

	public:
		/* Get a read-only reference to the m-th matrix-coefficient at the (i, j) position.
		 * @param m matrix-coefficient number, 0...length() - 1
		 * @param i Row index
		 * @param j Column index
		 * @returns Const reference to matrix entry
		 */
		_MatrixElement &getEntry (size_t m, size_t i, size_t j);
		
						/////////////////////////////////////
		                		//functions for matrix-coefficients//
						/////////////////////////////////////
						
	public:
		/* Set the m-th matrix-coefficient to V_m.
		 * @param m matrix-coefficient number, 0...length() - 1
		 * @param V_m matrix to set
		 */
		void setMatrixCoefficient (size_t m, const BlasMatrix<IntField, Rep> &V_m) ;

	private:
		/* Get a writeable reference to the m-th matrix-coefficient.
		 * @param m matrix-coefficient number, 0...length() - 1
		 * @returns Reference to matrix-coefficent
		 */
		BlasMatrix<IntField, Rep> &refMatrixCoefficient (size_t m) ;

	public:
		/** Get a read-only reference to the m-th matrix-coefficient
		 * @param m matrix-coefficient number, 0...length() - 1
		 * @returns Const reference to matrix-coefficent
		 */
		const BlasMatrix<IntField, Rep> &getMatrixCoefficient (size_t m) const ;

						/////////
		                		//swaps//
						/////////

	public:
		/* Swap i1-th and i2-th rows of matrices.
		 * This is done inplace.
		 */
		void swapRows(size_t i1, size_t i2);

		/* Swap j1-th and j2-th columns of matrices.
		 * This is done inplace.
		 */
		void swapCols(size_t j1, size_t j2);
		
								/////////////
		                				//transpose//
								/////////////

	public:
		/* Creates a transposed polynomial matrix of this.
		 * @param[in] tV
		 * @return the transposed polynomial matrix of this.
		 */
		Self_t transpose(Self_t & tV) const;
		
								//////////////////
		                				//input / output//
								//////////////////
	
	public:
		std::istream &read (std::istream &file);

		std::ostream &write (std::ostream &os)
		{
			int K = this->length();
			int I = this->rowdim();
			int J = this->coldim();
			for (int k = 0; k < K; k++)
			{
				for (int i = 0; i < I; i++)
				{
					for (int j = 0; j < J; j++)
					{
						os << this->getEntry(k, i, j) << " ";
					}
					os << std::endl;
				}
				os << std::endl;
			}
			return os;
		}
	};
}

/*#include "SlicedPolynomialMatrix.inl"*/
/**/
#include "/home/slb/fakeit/mySlicedPolynomialMatrix.inl" /**/

#endif
