/**
 *
 * coder@computer.org
 * coder@0xc0d3.org
 *
 * kinda note-pad foo
 */

#ifndef LE_SPACE_H
#define LE_SPACE_H

template< int _Dim, class _Topo, class _ElementT >
class TopologicalSpace {};
template< int _Dim, class _Topo, class _ElementT >
class UniformSpace: public TopologicalSpace< _Dim, _Topo, _ElementT > {};
template< int _Dim, class _Topo, class _Metric >: public UniformSpace < _Dim, _Topo, _Metric >
class MetricSpace {};

//

#endif
