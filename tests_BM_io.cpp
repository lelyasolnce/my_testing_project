#include <iostream>
#include <fstream>
#include "/home/slb/fakeit/myGF.h"
#include "/home/slb/fakeit/myBM.h"
#include "/home/slb/fakeit/mySlicedPolynomialMatrix.h"
#include "/home/slb/fakeit/my_modular.h"
#include "/home/slb/fakeit/my_bm_addsub.h"

typedef int64_t TT;
typedef Givaro::Modular<TT> Field;
typedef Givaro::GFqDom<TT> GField;
typedef LinBox::BlasMatrix<Field, std::vector<Field::Element>> BM;
int main()
{
	std::filebuf fb3; fb3.open("file3.txt", std::ios::out); std::ostream file3(&fb3);

	TT p = 7; Field F(p); int m = 5; int n = 3;
	BM bm1(F, m, n);
	for (int i = 0; i < m; i++) for (int j = 0; j < n; j++) bm1.setEntry(i, j, i+j);
	BM bm2(F, m, n);
	for (int i = 0; i < m; i++) for (int j = 0; j < n; j++) bm2.setEntry(i, j, 5);
	BM bm3(F, m, n);
	LinBox::BlasMatrixDomainAdd<Field, BM, BM, BM>()(F, bm3, bm1, bm2);
	file3<<"\nC = A + B\n";
	for (int i = 0; i < m; i++)
	{
		for (int j = 0; j < n; j++) file3<<bm3.getEntry(i, j)<<" ";
		file3<<"\n";
	}
	LinBox::BlasMatrixDomainSub<Field, BM, BM, BM>()(F, bm3, bm1, bm2);
	file3<<"\nC = A - B\n";
	for (int i = 0; i < m; i++)
	{
		for (int j = 0; j < n; j++) file3<<bm3.getEntry(i, j)<<" ";
		file3<<"\n";
	}
	LinBox::BlasMatrixDomainAddin<Field, BM, BM>()(F, bm3, bm1);
	file3<<"\nC += A\n";
	for (int i = 0; i < m; i++)
	{
		for (int j = 0; j < n; j++) file3<<bm3.getEntry(i, j)<<" ";
		file3<<"\n";
	}
	LinBox::BlasMatrixDomainSubin<Field, BM, BM>()(F, bm3, bm1);
	file3<<"\nC -= A\n";
	for (int i = 0; i < m; i++)
	{
		for (int j = 0; j < n; j++) file3<<bm3.getEntry(i, j)<<" ";
		file3<<"\n";
	}
}
