/**
 *
 * coder@computer.org
 * coder@0xc0d3.org
 *
 * kinda note-pad foo
 */

#ifndef LE_SPACE_H
#define LE_SPACE_H

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


template< int _Dim, class _Topo, class _ElementT >
class TopologicalSpace {};
template< int _Dim, class _Topo, class _ElementT >
class UniformSpace: public TopologicalSpace< _Dim, _Topo, _ElementT > {};
template< int _Dim, class _Topo, class _Metric >
class MetricSPace: public UniformSpace < _Dim, _Topo, _Metric > {};

//kinda abstract simplex(?), one example for an element
template< int _Dim, LinkType _LType, AccessScheme _AScheme,  template< class U, class V > class _Containment, template< class U > class _Allocator > class AbstractSimplex {};

template< int _Dim,  template< class U, class V > class _Containment, template< class U > class _Allocator >
class AbstractSimplex< _Dim, LinkType::Single, AccessScheme::IndexAccessScheme, _Containment, _Allocator >
{
    public:
        ptrdiff_t upper, opponent, next;
};
//

#endif
