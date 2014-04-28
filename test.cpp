#include "le_space.h"
#include <vector>
#include <iostream>

int main()
{
    AbstractSimplex< 3, LinkType::Single, AccessScheme::Index, std::vector, std::allocator > tetra;
    std::cout << "atetra::d = " << decltype(tetra)::d << std::endl;
    return 0;
}
