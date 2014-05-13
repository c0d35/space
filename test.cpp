#include "le_space.h"
#include <vector>
#include <iostream>

int main()
{


    SimpleEuklidianSpace< 6 >::Vector v;
	AbstractSimplex< 3, LinkType::Single, AccessScheme::Index, std::vector, std::allocator, Set< 3 >  > atetra;
	Simplex< 6 > simplex6;
	std::cout << "atetra::d = " << decltype(atetra)::d << std::endl;
	std::cout << "simplex6::d = " << decltype(simplex6)::d << std::endl;
    std::cout << " SimpleEuklidianMetricTrait< 6, float >::Vector::d  = " << decltype(v)::d << std::endl;
    std::cout << "ibinom< 6, 3 >::eval = " << ibinom< 6, 3 >::eval << std::endl;
    std::cout << "ipow< 2, 10 >::eval = " << ipow< 2, 10 >::eval << std::endl;
    std::cout << "ifact< 6 >::eval = " << ifact< 6 >::eval << std::endl;
	return 0;
}
