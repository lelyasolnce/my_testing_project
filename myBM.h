#ifndef __myBM_H
#define __myBM_H

#include <iostream>

namespace LinBox
{ /*  Blas Matrix */
	

	/*! Dense matrix representation.
	 * @ingroup matrix
	 * A \p BlasMatrix is a matrix of \p _Field::Element, with the structure of BLAS matrices.
	 * It is basically a vector of \p _Field::Element.
	 * In the Mother model, a \p BlasMatrix is allocated by the user.
	 *@bug why not BlasMatrixDomain ?
	 */
	template <class _Field, class _Storage>
	class BlasMatrix {
		// private :

	public:
		typedef _Field                                  Field;
		typedef typename Field::Element               Element;    //!< Element type
		typedef _Storage                                  Rep;    //!< Actually a <code>std::vector<Element></code> (or alike.)
		typedef BlasMatrix<Field,Rep>                  Self_t;    //!< Self typeype
		typedef const BlasMatrix<Field,Rep>       constSelf_t;    //!< Self typeype

                //typedef BlasSubmatrix<Self_t>           subMatrixType;    //!< Submatrix type
		//typedef BlasSubmatrix<constSelf_t> constSubMatrixType;    //!< Submatrix type
                typedef Self_t                             matrixType;    //!< matrix type
                typedef constSelf_t                   constMatrixType;    //!< matrix type
                typedef Self_t                               blasType;    //!< blas matrix type

	protected:
		size_t _row;
		size_t _col;
		Rep	_rep;
	public:
	       	const Field		    * _field; //! @bug why public ?
		//MatrixDomain<Field>    _MD; //! @bug why public ?
		//VectorDomain<Field>    _VD;


	public:

		//////////////////
		// CONSTRUCTORS //
		//////////////////


		/*! Allocates a new \f$ 0 \times 0\f$ matrix (shaped and ready).*/
		BlasMatrix (const _Field &F)
		{
			_row = 0; _col = 0; _field = &F;
		}
		

		/*! Allocates a new \f$ m \times n\f$ zero matrix (shaped and ready).		 */
		BlasMatrix (const _Field &F, const size_t & m, const size_t &n)
		{
			//_row = m; _col = n; _rep(_row*_col, F.zero); _filed = &F;
			_row = m; _col = n;
			for (size_t j = 0; j < _row * _col; j++)
				_rep.push_back(F.zero);
			_field = &F;
		}

		/*
		// Allocates a new bare \f$ 0 \times 0\f$ matrix (unshaped, unready).//
		// BlasMatrix () ;

		/// (Re)allocates a new \f$ m \times n\f$ zero matrix (shaped and ready).
		void init(const _Field & F, const size_t & r = 0, const size_t & c = 0);

		//! Generic copy constructor from either a blackbox or a matrix container.
		template <class Matrix>
		BlasMatrix (const Matrix &A) ;

		//! Create a BlasMatrix from a vector of elements
		BlasMatrix (const _Field &F, const std::vector<Element>& v,
			    const size_t &m , const size_t &n) ;
		*/

		/// Destructor.
		~BlasMatrix () {}

		//////////////////
		//  DIMENSIONS  //
		//////////////////

		size_t rowdim() const {return _row;}
		size_t coldim() const {return _col;}
		size_t getStride() const {return _col;}
		size_t stride() const { return getStride() ;}

		//////////////////
		//   ELEMENTS   //
		//////////////////

	
		void setEntry (size_t i, size_t j, const Element &a_ij) {_rep[i * _col + j] = a_ij;}
		Element &refEntry (size_t i, size_t j) {return _rep[i * _col + j];}
		Element &getEntry (size_t i, size_t j) {return _rep[i * _col + j];}
		Element &getEntry (Element &x, size_t i, size_t j) {return _rep[i * _col + j];}

		///////
		//i/o//
		///////
		void write()
		{
			for (size_t i = 0; i < _row; i++)
			{
				for (size_t j = 0; j < _col; j++) std::cout<<_rep[i*_col+j]<<' ';
				std::cout<<'\n';
			}
		}

		///////////////////
		// TRANSPOSE &AL //
		///////////////////

		/*! Creates a transposed matrix of \c *this.
		 * @param[in] tM
		 * @return the transposed matrix of this.
		 */
		//Self_t transpose(Self_t & tM) const ;
		//void transpose() ; //inplace

	/*
		void random()
		{
			subMatrixType B(*this, 0, 0, rowdim(), coldim());
			B.random();
		}

		template<class Rand>
		void random(const Rand&)
		{
			return random();
		}
		*/
	}; // end of class BlasMatrix
} // end of namespace LinBox

#endif
