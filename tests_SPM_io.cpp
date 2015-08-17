#include <iostream>
#include <fstream>
#include "/home/slb/fakeit/myGF.h"
#include "/home/slb/fakeit/myBM.h"
#include "/home/slb/fakeit/mySlicedPolynomialMatrix.h"
#include "/home/slb/fakeit/my_modular.h"
#include "/home/slb/fakeit/my_addsub.h"

typedef int64_t TT;
typedef Givaro::Modular<TT> Field;
typedef Givaro::GFqDom<TT> GField;
typedef LinBox::BlasMatrix<Field, std::vector<Field::Element>> BM;
typedef LinBox::SlicedPolynomialMatrix<GField, TT> SPM;
int main()
{
	GField _gf(7, 3); int m = 2, n = 5;
	std::filebuf fb1; fb1.open("file1.txt", std::ios::in); std::istream file1(&fb1);
	std::filebuf fb2; fb2.open("file2.txt", std::ios::in); std::istream file2(&fb2);
	std::filebuf fb3; fb3.open("file3.txt", std::ios::out); std::ostream file3(&fb3);
	SPM spm1(_gf, m, n); SPM spm2(_gf, m, n); spm1.read(file1); spm2.read(file2);
	SPM spm3(_gf, m, n); 
	
	
	LinBox::SlicedPolynomialMatrixAdd<GField, SPM, SPM, SPM>()(_gf, spm3, spm1, spm2);
	file3<<"\nC = A + B\n"; spm3.write(file3);
	LinBox::SlicedPolynomialMatrixSub<GField, SPM, SPM, SPM>()(_gf, spm3, spm1, spm2);
	file3<<"\nC = A - B\n"; spm3.write(file3);
	LinBox::SlicedPolynomialMatrixAddin<GField, SPM, SPM>()(_gf, spm3, spm1);
	file3<<"\nC += A\n"; spm3.write(file3);
	LinBox::SlicedPolynomialMatrixSubin<GField, SPM, SPM>()(_gf, spm3, spm1);
	file3<<"\nC -= A\n"; spm3.write(file3);
}
