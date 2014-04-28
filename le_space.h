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
	PointerAccessScheme,
	IndexAccessScheme,
};


template< int _Dim, typename _Topo, typename _ElementT >
class TopologicalSpace {};
template< int _Dim, typename _Topo, typename _ElementT >
class UniformSpace: public TopologicalSpace< _Dim, _Topo, _ElementT > {};
template< int _Dim, typename _Topo, typename _Metric >
class MetricSpace: public UniformSpace < _Dim, _Topo, _Metric > {};

/**
 * kinda abstract simplex(?), one example for an element
 */

template< int _Dim, LinkType _LType, AccessScheme _AScheme,  template< class U, class V > class _Containment, template< class U > class _Allocator > struct AbstractSimplex {};


///why Dim + 1 you may ask? 'cause \f$ {d+1 \choose d} = d + 1\f$ , u kno?
template< int _Dim,  template< class U, class V > class _Containment, template< class U > class _Allocator >
struct AbstractSimplex< _Dim, LinkType::Single, AccessScheme::IndexAccessScheme, _Containment, _Allocator >
{
    public:
        enum {d = _Dim};
        ptrdiff_t upper, opponent, next;
        ptrdiff_t lower[_Dim + 1]; 
};

//

#endif
