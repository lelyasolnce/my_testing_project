#include <iostream>
#include <fstream>
#include "/home/slb/fakeit/myGF.h"
#include "/home/slb/fakeit/my_modular.h"
#include "/home/slb/fakeit/GivPol.h"

typedef int64_t TT;
typedef Givaro::Modular<TT> Field;
typedef Givaro::GFqDom<TT> GField;
typedef Givaro::Poly1Dom<Field> PolyDom;
typedef typename PolyDom::Rep polynomial;

int main()
{
	TT p = 7; Field F(p);
	PolyDom polydom(F);
	polynomial _irred;
	for (int i = 0; i < 4; i++) _irred.push_back(0); _irred.push_back(1);
	polynomial entry;
	for (int i = 0; i < 10; i++) entry.push_back((i*2 + 5) % p);
	polynomial result;
	polydom.mod(result, entry, _irred);
for (int i = 0; i < result.size(); i++) std::cout<<result[i]<<" ";
std::cout<"done";
}
