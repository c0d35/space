#include "le_space.h"
#include <vector>
#include <iostream>

int main()
{
	AbstractSimplex< 3, LinkType::Single, AccessScheme::Index, std::vector, std::allocator, Set< 3 >  > atetra;
	AbstractSimplex< 6, LinkType::Single, AccessScheme::Index, std::vector, std::allocator, TopologicalSpace< 6 > > simplex6;
	std::cout << "atetra::d = " << decltype(atetra)::d << std::endl;
	std::cout << "simplex6::d = " << decltype(simplex6)::d << std::endl;
	return 0;
}
