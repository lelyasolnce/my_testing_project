#ifndef __LINBOX_matrix_SlicedPolynomialMatrix_SlicedPolynomialMatrixMulKaratsuba_H
#define __LINBOX_matrix_SlicedPolynomialMatrix_SlicedPolynomialMatrixMulKaratsuba_H

//#include "SlicedPolynomialMatrix.h"
#include "/home/slb/fakeit/mySlicedPolynomialMatrix.h"
#include "/home/slb/fakeit/myBM.h"
#include <vector>

namespace LinBox
{ 
	template< class GField, class Operand1, class Operand2, class Operand3>
	class SlicedPolynomialMatrixMulKaratsuba
	{
	private:
		typedef typename Operand1::IntField IF;
		typedef typename Operand1::MatrixElement ME;
		typedef BlasMatrix<IF, std::vector<ME>> BM;
		typedef std::vector<BM> vec;
		typedef typename Operand1::polynomial polynomial;
		typedef Givaro::Poly1Dom<IF> PolyDom;
		vec& karatsuba(IF& F, vec& A, vec& B, vec& C);
	public:
		Operand1 &operator() (GField &GF, Operand1 &C, Operand2 &A, Operand3 &B);
	}; 
} /* end of namespace LinBox */

//#include "SlicedPolynomialMatrixMulKaratsuba.inl"
#include "/home/slb/fakeit/my_karatsuba.inl"

#endif
