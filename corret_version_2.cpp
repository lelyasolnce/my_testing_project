#include <givaro/modular-double.h>
#include <givaro/modular-int64.h>
#include <iostream>
#include "/home/slb/fakeit/myGF.h"
#include "/home/slb/fakeit/myBM.h"

typedef int64_t TT;
typedef Givaro::Modular<TT> Field;
typedef GF<TT> GField; //instead of Givaro::GFqDom<TT>
int main()
{
	TT p = 7;
	Field F(p);

	GField _gf(5, 3, 0); //GField _gf(5, 3);
	int m = 7, n = 9;
	LinBox::BlasMatrix<Field, std::vector<Field::Element>> _bm(F, m, n);
std::cout<<_bm.coldim();
int el(6);
_bm.setEntry(3, 7, el);
}
