#include "le_space.h"
#include <vector>
#include <iostream>
typedef HyperCubeTree< 3, 0, SimpleEuklidianSpaceFuint16 , std::vector, std::allocator > HCTree3;
//typedef HyperCubeTree< 3, 0, unsigned short, SimplePoint, DefaultEuklidianMetricTrait, std::vector, std::allocator >::PointT HCTree3Point;
int main()
{
    
    HCTree3 tree;
    SimpleEuklidianSpaceFdouble< 6 >::Vector v;
    EuklidianSpaceCompressedFuint16< 4 > space;
    EuklidianSpaceCompressedFuint16< 4 >::PointT point = {1, 2, 3, 4};
    //typedef EuklidianSpaceCompressedFuint16< 4 > Space3;
    space.insert(point);

	AbstractHalfSimplex< 3, LinkType::Single, AccessScheme::Index, std::vector, std::allocator, Set< 3 >  > atetra;
	HalfSimplex< 6 > simplex6;
	std::cout << "atetra::d = " << decltype(atetra)::d << std::endl;
	std::cout << "simplex6::d = " << decltype(simplex6)::d << std::endl;
    std::cout << " SimpleEuklidianMetricTrait< 6, float >::Vector::d  = " << decltype(v)::d << std::endl;
    std::cout << "nchoosek< 6, 3 >::eval = " << nchoosek< 6, 3 >::eval << std::endl;
    std::cout << "ipow< 2, 10 >::eval = " << ipow< 2, 10 >::eval << std::endl;
    std::cout << "ifact< 6 >::eval = " << ifact< 6 >::eval << std::endl;
    //std::cout << "space.simplicial_complex[1] = " << space.simplicial_complex[
    //    (SimplicialComplex<3>::Simplex<1>::d)] << std::endl;
	return 0;
}
