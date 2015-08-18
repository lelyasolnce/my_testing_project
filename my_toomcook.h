#ifndef __LINBOX_matrix_SlicedPolynomialMatrix_SlicedPolynomialMatrixMulToomCook_H
#define __LINBOX_matrix_SlicedPolynomialMatrix_SlicedPolynomialMatrixMulToomCook_H

//#include "SlicedPolynomialMatrix.h"
#include "/home/slb/fakeit/mySlicedPolynomialMatrix.h"
#include "/home/slb/fakeit/myBM.h"

namespace LinBox
{ 
	template< class GField, class Operand1, class Operand2, class Operand3>
	class SlicedPolynomialMatrixMulToomCook
	{
	private:
		typedef typename Operand1::IntField IntField;
		typedef Givaro::Poly1Dom<IntField> PolyDom;
		typedef typename Operand1::Rep Rep;
		typedef BlasMatrix<IntField, Rep> Matrix;
		typedef typename Operand1::polynomial polynomial;
		Matrix& EvaluationInterpolationMatrices (Matrix& TC, Matrix& iTC);
		Matrix& mul (IntField& F, Matrix& CMatBloc,  Matrix& AMatBloc,  Matrix& BMatBloc,
								    size_t m,  size_t k,  size_t n,  size_t e, polynomial irreducible);
	public:
		Operand1 &operator() ( GField &GF, Operand1 &C,  Operand2 &A,  Operand3 &B) ;
	}; 
} /* end of namespace LinBox */

//#include "SlicedPolynomialMatrixMulToomCook.inl"
#include "/home/slb/fakeit/my_toomcook.inl"

#endif

