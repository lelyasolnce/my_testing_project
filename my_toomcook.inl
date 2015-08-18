#ifndef __LINBOX_matrix_SlicedPolynomialMatrix_SlicedPolynomialMatrixMulToomCook_INL
#define __LINBOX_matrix_SlicedPolynomialMatrix_SlicedPolynomialMatrixMulToomCook_INL

//#include <fflas-ffpack/fflas/fflas.h>
#include "/home/slb/fakeit/my_fflas.h"
//#include <fflas-ffpack/ffpack/ffpack.h>
#include "/home/slb/fakeit/my_ffpack.h"
//#include "linbox/algorithms/matrix-hom.h"
//#include <givaro/givpoly1denseops.inl>
#include "/home/slb/fakeit/GivPol.h"
//#include "linbox/matrix/SlicedPolynomialMatrix/conversion.h"
#include <math.h>

namespace LinBox
{
	template<class GField, class Operand1, class Operand2, class Operand3>
	BlasMatrix<typename Operand1::IntField, typename Operand1::Rep>&
	SlicedPolynomialMatrixMulToomCook<GField, Operand1, Operand2, Operand3 >::EvaluationInterpolationMatrices
				(Matrix& TC, Matrix& iTC)
	{
		size_t E = TC.rowdim();
		for (size_t i = 0 ; i < E ; ++i)
		{
			for (size_t j = 0 ; j < E ; ++j)
			{
				TC.setEntry(i, j, pow(i, j));
			}
		}
		int null;
		FFPACK::Invert<IntField>(TC.field(),E,TC.getPointer(),E,iTC.getWritePointer(),E,null);
		return TC;
	}

	template<class GField, class Operand1, class Operand2, class Operand3>
		BlasMatrix<typename Operand1::IntField, typename Operand1::Rep>&
		SlicedPolynomialMatrixMulToomCook<GField, Operand1, Operand2, Operand3 >::mul
			(IntField& F, Matrix& CMatBloc, Matrix& AMatBloc, Matrix& BMatBloc,
								    size_t m,
								    size_t k,
								    size_t n,  size_t e,
								   polynomial irreducible)
	{
		//#if (__LINBOX_FFLAS_FFPACK_VERSION < 10501)F
				//#warning "Invert is buggy in your fflas-ffpack version. please consider upgrading to >=1.5.1."
		//#endif
		size_t E = 2*e - 1 ;

		Matrix TC    (F, E, E);
		Matrix iTC   (F, E, E);
		Matrix iEval (F, E, E);
		EvaluationInterpolationMatrices(TC, iTC);

		// each row is a result matrix
		Matrix TMatBloc( F, E, m*n);
		Matrix AEval( F , E, m*k);
		Matrix BEval( F , E, k*n);

		FFLAS::fgemm(F,
					 FFLAS::FflasNoTrans, FFLAS::FflasNoTrans,
					 E, m*k, e,
					 F.one,
					 TC.getPointer(),E,
					 AMatBloc.getPointer(), m*k,
					 F.zero,
					 AEval.getWritePointer(), m*k);

		FFLAS::fgemm(F,
					 FFLAS::FflasNoTrans, FFLAS::FflasNoTrans,
					 E ,n*k, e,
					 F.one,
					 TC.getPointer(),E,
					 BMatBloc.getPointer(), n*k,
					 F.zero,
					 BEval.getWritePointer(), n*k);

		for (size_t i = 0 ; i < E ; ++i)
		{
			FFLAS::fgemm(F,
						 FFLAS::FflasNoTrans, FFLAS::FflasNoTrans,
						 m, n, k,
						 F.one,
						 AEval.getPointer()+i*m*k, k,
						 BEval.getPointer()+i*n*k, n,
						 F.zero,
						 TMatBloc.getWritePointer()+i*m*n, n);


		FFLAS::fgemm(F,
					 FFLAS::FflasNoTrans, FFLAS::FflasNoTrans,
				     E, m * n, E,
					 F.one,
					 iTC.getPointer(),E, //lda
					 TMatBloc.getWritePointer()+i*m*n, m*n, //ldb
					 F.zero,
					 CMatBloc.getWritePointer()+i*m*n, m*n);
		}

		return CMatBloc;
	}

	// all matrix classes should be SlicedPolynomialMatrices
	template<class GField, class Operand1, class Operand2, class Operand3>
	Operand1& SlicedPolynomialMatrixMulToomCook<GField, Operand1, Operand2, Operand3 >::operator()
									   ( GField& GF,
									   Operand1& C,
									    Operand2& A,
									    Operand3& B)
	{
		PolyDom polydom(C.fieldF());
		size_t e = C.length();
		size_t m = C.rowdim();
		size_t k = B.rowdim();
		size_t n = C.coldim();

		IntField F = A.fieldF();

		if (e == 1) {
			Matrix Am(F, A.rowdim(), A.coldim());
			Am = A.getMatrixCoefficient(0);//conversionAtoAm()(Am, A);
			Matrix Bm(F, B.rowdim(), B.coldim());
			Bm = B.getMatrixCoefficient(0); //conversionBtoBm()(Bm, B);
			Matrix Cm(F, C.rowdim(), C.coldim());
			FFLAS::fgemm(F,
						 FFLAS::FflasNoTrans, FFLAS::FflasNoTrans,
						 m,n,k,
						 F.one,
						 Am.getPointer(), Am.getStride(), //lda
						 Bm.getPointer(), Bm.getStride(), //ldb
						 F.zero,
						 Cm.getWritePointer(), Cm.getStride());
			C.setMatrixCoefficient(0, Cm);
			return C;
		}

		Matrix Cbloc(F,e,m*n);
		Matrix Abloc(F,e,m*k);
		Matrix Bbloc(F,e,k*n);


		for (size_t l = 0 ; l < e ; ++l)
		{
			for (size_t i = 0 ; i < m ; ++i)
			{
				for (size_t j = 0 ; j < k ; ++j)
				{
					Abloc.setEntry(l, i*k+j, A.getEntry(l, i, j));
				}
			}
		}

		for (size_t l = 0 ; l < e ; ++l)
		{
			for (size_t i = 0 ; i < k ; ++i)
			{
				for (size_t j = 0 ; j < n ; ++j)
				{
					Bbloc.setEntry(l, i*n+j, B.getEntry(l, i, j));
				}
			}
		}

		mul(C.fieldF(), Cbloc, Abloc, Bbloc, m, k, n, e, C.irreducible());
		
		size_t E = 2 * E - 1;
		polynomial x; for (int l = 0; l < E; l++) x.push_back(F.zero);
		for (size_t i = 0 ; i < m ; ++i)
		{
			for (size_t j = 0 ; j < n ; ++j)
			{
				for (size_t l = 0 ; l < E ; ++l)
				{
					x[l] = Cbloc.getEntry(l, i*n+j);
				}
				polydom.modin(x, C.irreducible());
				for (size_t l = 0 ; l < e ; ++l)
				{
					C.setEntry(l, i, j, x[l]);
				}
			}
		}

		return C ;
	}
} // LinBox

#endif

