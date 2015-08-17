#include <iostream>
#include <fstream>
#include "/home/slb/fakeit/myGF.h"
#include "/home/slb/fakeit/myBM.h"
#include "/home/slb/fakeit/mySlicedPolynomialMatrix.h"
#include "/home/slb/fakeit/my_modular.h"

typedef int64_t TT;
typedef Givaro::Modular<TT> Field;
typedef Givaro::GFqDom<TT> GField;
typedef LinBox::BlasMatrix<Field, std::vector<Field::Element>> BM;
typedef LinBox::SlicedPolynomialMatrix<GField, TT> SPM;

int main()
{
	GField _gf(7, 3);

	//constructor without dimensions
	SPM spm1(_gf);
	
	//constructor with dimensions
	int m = 2, n = 5;
	SPM spm2(_gf, m, n);

	//input
	std::filebuf fb1; fb1.open("file1.txt", std::ios::in); std::istream file1(&fb1);
	spm2.read(file1);

	//output
	std::filebuf fb2; fb2.open("file2.txt", std::ios::out); std::ostream file2(&fb2);
	spm2.write(file2);

	//dimensions
	file2<<"rowdim = "<<spm2.rowdim()<<"; coldim = "<<spm2.coldim()<<"; length = "<<spm2.length();

	//access to the entries
	TT el; el = spm2.getEntry(2, 0, 0); spm2.setEntry(2, 1, 4, el);
	file2<<"\nafter setEntry\n"; spm2.write(file2);

	//swaps
	int i1 = 0, i2 = 1; spm2.swapRows(i1, i2);
	file2<<"\nswapped rows "<<i1<<" and "<<i2<<"\n"; spm2.write(file2);
	int j1 = 2, j2 = 4; spm2.swapCols(j1, j2);
	file2<<"\nswapped columns "<<j1<<" and "<<j2<<"\n"; spm2.write(file2);

	//access to the matrix coefficients
	TT p = 7; Field F(p);
	BM bm2(F, m, n);
	for (int i = 0; i < m; i++) for (int j = 0; j < n; j++) bm2.setEntry(i, j, i+j);
	int m0 = 2; spm2.setMatrixCoefficient(m0, bm2); bm2 = spm2.getMatrixCoefficient(m0);
	file2<<"\nchanged leading matrix coefficient\n"; spm2.write(file2);
}
