#include "le_space.h"
#include <vector>
#include <iostream>

int main()
{


    SimpleEuklidianSpace< 6 >::Vector< double > v;
	AbstractSimplex< 3, LinkType::Single, AccessScheme::Index, std::vector, std::allocator, Set< 3 >  > atetra;
	Simplex< 6 > simplex6;
	std::cout << "atetra::d = " << decltype(atetra)::d << std::endl;
	std::cout << "simplex6::d = " << decltype(simplex6)::d << std::endl;
    std::cout << " SimpleEuklidianMetricTrait< 6, float >::Vector::d  = " << decltype(v)::d << std::endl;
	return 0;
}
