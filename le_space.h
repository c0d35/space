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

template< int _Dim, LinkType _LType, AccessScheme _AScheme,  template< class U, class V > class _Containment, template< class U > class _Allocator > struct AbstractSimplex {};


///why Dim + 1 you may ask? 'cause \f$ {d+1 \choose d} = d + 1\f$ , u kno?
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

//

#endif
