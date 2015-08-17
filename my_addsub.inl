#ifndef __LINBOX_matrix_SlicedPolynomialMatrix_SlicedPolynomialMatrixAddSub_INL
#define __LINBOX_matrix_SlicedPolynomialMatrix_SlicedPolynomialMatrixAddSub_INL

#include <vector>
//#include "linbox/matrix/MatrixDomain/blas-matrix-domain.h"
#include "/home/slb/fakeit/my_bm_addsub.h"

namespace LinBox
{
	template<class GField, class Vector1, class Vector2, class Vector3>
	Vector1& SlicedPolynomialMatrixAdd<GField, Vector1, Vector2, Vector3 >::operator()(GField& F, Vector1& C, Vector2& A, Vector3& B)
	{
		BM1 M(C.fieldF(), C.rowdim(), C.coldim());
		BM2 M2(A.fieldF(), A.rowdim(), A.coldim());
		BM3 M3(B.fieldF(), B.rowdim(), B.coldim());
		for (int m = 0; m < C.length(); m++)
		{
			M2 =  A.getMatrixCoefficient(m);
			M3 =  B.getMatrixCoefficient(m);
			BlasMatrixDomainAdd<IF1, BM1, BM2, BM3>()(C.fieldF(), M, M2, M3);
			C.setMatrixCoefficient(m, M);
		}
		return C;
	}

	template<class GField, class Vector1, class Vector2, class Vector3>
	Vector1& SlicedPolynomialMatrixSub<GField, Vector1, Vector2, Vector3 >::operator()(GField& F, Vector1& C, Vector2& A, Vector3& B)
	{
		BM1 M(C.fieldF(), C.rowdim(), C.coldim());
		BM2 M2(A.fieldF(), A.rowdim(), A.coldim());
		BM3 M3(B.fieldF(), B.rowdim(), B.coldim());
		for (int m = 0; m < C.length(); m++)
		{
			M2 =  A.getMatrixCoefficient(m);
			M3 =  B.getMatrixCoefficient(m);
			BlasMatrixDomainSub<IF1, BM1, BM2, BM3>()(C.fieldF(), M, M2, M3);
			C.setMatrixCoefficient(m, M);
		}
		return C;
	}

	template<class GField, class Vector1, class Vector2>
	Vector1& SlicedPolynomialMatrixAddin<GField, Vector1, Vector2>::operator()(GField& F, Vector1& C, Vector2& A)
	{
		BM1 M1(C.fieldF(), C.rowdim(), C.coldim());
		BM2 M2(A.fieldF(), A.rowdim(), A.coldim());
		for (int m = 0; m < C.length(); m++)
		{
			M1 =  C.getMatrixCoefficient(m);
			M2 =  A.getMatrixCoefficient(m);
			BlasMatrixDomainAddin<IF1, BM1, BM2>()(C.fieldF(), M1, M2);
			C.setMatrixCoefficient(m, M1);
		}
		return C;
	}

	template<class GField, class Vector1, class Vector2>
	Vector1& SlicedPolynomialMatrixSubin<GField, Vector1, Vector2>::operator()(GField& F, Vector1& C, Vector2& A)
	{
		BM1 M1(C.fieldF(), C.rowdim(), C.coldim());
		BM2 M2(A.fieldF(), A.rowdim(), A.coldim());
		for (int m = 0; m < C.length(); m++)
		{
			M1 =  C.getMatrixCoefficient(m);
			M2 =  A.getMatrixCoefficient(m);
			BlasMatrixDomainSubin<IF1, BM1, BM2>()(C.fieldF(), M1, M2);
			C.setMatrixCoefficient(m, M1);
		}
		return C;
	}
} // LinBox

#endif
