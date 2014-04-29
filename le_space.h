/**
 *
 * coder@computer.org
 * coder@0xc0d3.org
 *
 * kinda note-pad foo
 */

#ifndef LE_SPACE_H
#define LE_SPACE_H

#include <cstddef>


enum class LinkType: bool
{
	Single,
	Double,
};

enum class AccessScheme: bool
{
	Pointer,
	Index,
};


template< int _Dim, typename _Topo >
struct TopologicalSpace {};
template< int _Dim, typename _Topo  >
struct UniformSpace: TopologicalSpace< _Dim, _Topo > {};
template< int _Dim, typename _Metric >
struct MetricSpace: UniformSpace < _Dim, _Metric > {};

/**
 * kinda abstract simplex(?), one example for an element
 */
template< int _Dim, typename _Trait, LinkType _LType, AccessScheme _AScheme, template< class U , class V > class Containment , template< class U > class Allocator> class AbstractSimplicialComplex{};
template< typename _ASD > struct AbstractSimplicialComplexTopologyTrait{};
template< typename _Iterator > struct AbstractSimplicialComplexIteratorFunctionTrait{};
template< typename _ASD > struct AbstractSimplicialComplexIterator{};
template< int _Dim, LinkType _LType, AccessScheme _AScheme,  template< class U, class V > class _Containment, template< class U > class _Allocator > struct AbstractSimplex {};


///why d + 1 you may ask. 'cause \f$ {d+1 \choose d} = d+1\f$ , u kno'. thus we're terminating the recursion at dimension = -1
template< int _Dim,  template< class U, class V > class _Containment, template< class U > class _Allocator >
struct AbstractSimplex< _Dim, LinkType::Single, AccessScheme::Index, _Containment, _Allocator >: TopologicalSpace< _Dim, AbstractSimplex< _Dim, LinkType::Single, AccessScheme::Index, _Containment, _Allocator > >
{
        enum {d = _Dim};
        ptrdiff_t upper, opponent, next;
        ptrdiff_t lower[_Dim + 1]; 
};

///terminate the recursion in the empty simplex (simplicial set) with dimension -1
template< template< class U, class V > class _Containment, template< class U > class _Allocator >
struct AbstractSimplex< -1, LinkType::Single, AccessScheme::Index, _Containment, _Allocator >: TopologicalSpace< -1, AbstractSimplex< -1, LinkType::Single, AccessScheme::Index, _Containment, _Allocator > >
{
        enum { d = -1, };
};

template< int _D, class _ASD > struct ContainerFiller
{ 
	static inline void fill(_ASD& s)
	{
		typedef AbstractSimplex< _D - 1, LinkType::Single, AccessScheme::Index, _ASD::template _Containment, _ASD::template _Allocator > AbstractSimplexT;
		typedef typename _ASD::_Trait::template SimplexT< _D, AbstractSimplexT > Simplex;
		typedef typename _ASD::template _Containment< Simplex, _ASD::template _Allocator< Simplex > > Simplices;
		s.simplex_containers[_D - 1] = new Simplices;
		ContainerFiller< _D - 1, _ASD >::fill(s);

	};	

};



template< class _ASD > struct ContainerFiller< 0, _ASD >{ static inline void fill(_ASD&){}};


template< int _D, class _IT, class  _ASD > struct IterFiller
{
	static inline void fill(_IT& i)
	{
		typedef  AbstractSimplex< _D - 1, LinkType::Single, AccessScheme::Index, _ASD::template _Containment, _ASD::template _Allocator > AbstractSimplexT;
		i.iterdata[_D] = new AbstractSimplexT;
		IterFiller< _D - 1, _IT, _ASD >::fill(i);
	}
};



//

#endif
