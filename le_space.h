/**
 *
 * coder@computer.org
 * coder@0xc0d3.org
 *
 * this is only a sketch for some foo i'm thinking about
 */

//abstract simplex with topological space as default space
//metrical traits for every space (quasi-, pseudo-, semi- metrics included ... divergences & stuff)

#ifndef LE_SPACE_H
#define LE_SPACE_H

#include <cstddef>
//some helpers, rewrite them with constexpr

template< bool _cond, class _then, class _else >
struct IF
{
	typedef _then RET;
};

template< class _then, class _else >
struct IF< false, _then, _else >
{
	typedef _else RET;
};

template< int N, int M >
struct ipow
{
	enum { eval = N * ipow< N , M - 1 >::eval};
};

template < int N >
struct ipow< N, 0 >
{
	enum { eval = 1 };
};

template< typename _T, int _I >
struct simd_sum
{
	static inline _T eval(_T &a)
	{
		return a[_I] + simd_sum< _T, _I - 1 >::eval(a);
	}
};

template< typename _T >
struct simd_sum< _T, 0 >
{
	static inline _T eval(_T &a)
	{
		return a[0];
	}
};

template< typename _T, int _I >
struct simd_cp
{
	static inline void eval(_T &d, _T &s)
	{
		d[_I] = s[_I]; simd_cp< _T, _I - 1 >::eval(d, s);
	}
};

template< typename _T >
struct simd_cp< _T, 0 >
{
	static inline void eval(_T &d, _T &s)
	{
		return d[0] = s[0];
	}
};


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
//
//some structures

template< int N > struct Algebra {};
template< int M, int N > struct GradedAlgebra {};
template< int _Dim >
struct Set {};
template< int _Dim >
struct TopologicalSpace {};
template< int _Dim >
struct UniformSpace: TopologicalSpace< _Dim > {};
template< int _Dim >
struct MetricSpace: UniformSpace < _Dim > {};

/**
 * kinda abstract simplex(?), one example for an element
 */
template< int _Dim, typename _Trait, LinkType _LType, AccessScheme _AScheme, template< class U , class V > class Containment , template< class U > class Allocator, class _Space > class AbstractSimplicialComplex{};
template< typename _ASD > struct AbstractSimplicialComplexTopologyTrait{};
template< typename _Iterator > struct AbstractSimplicialComplexIteratorFunctionTrait{};
template< typename _ASD > struct AbstractSimplicialComplexIterator{};
template< int _Dim, LinkType _LType, AccessScheme _AScheme,  template< class U, class V > class _Containment, template< class U > class _Allocator, class _Space > struct AbstractSimplex: _Space {};

///why d + 1 you may ask. 'cause \f$ {d+1 \choose d} = d+1\f$ , u kno'. thus we're terminating the recursion at dimension = -1
template< int _Dim,  template< class U, class V > class _Containment, template< class U > class _Allocator, template < int D > class _Space >
struct AbstractSimplex< _Dim, LinkType::Single, AccessScheme::Index, _Containment, _Allocator, _Space< _Dim > >: _Space < _Dim >
{
    enum {d = _Dim};
    ptrdiff_t upper, opponent, next;
    ptrdiff_t lower[_Dim + 1]; 
};

///terminate the recursion in the empty simplex (simplicial set) with dimension -1
template< template< class U, class V > class _Containment, template< class U > class _Allocator, template < int D > class _Space >
struct AbstractSimplex< -1, LinkType::Single, AccessScheme::Index, _Containment, _Allocator, _Space< -1 > >: _Space< -1 >
{
    enum { d = -1, };
};


///why d + 1 you may ask. 'cause \f$ {d+1 \choose d} = d+1\f$ , u kno'. thus we're terminating the recursion at dimension = -1
template< int _Dim,  template< class U, class V > class _Containment, template< class U > class _Allocator >
struct AbstractSimplex< _Dim, LinkType::Single, AccessScheme::Index, _Containment, _Allocator, Set< _Dim > >: Set < _Dim >
{
    enum {d = _Dim};
    ptrdiff_t upper, opponent, next;
    ptrdiff_t lower[_Dim + 1]; 
};

///terminate the recursion in the empty simplex (simplicial set) with dimension -1
template< template< class U, class V > class _Containment, template< class U > class _Allocator >
struct AbstractSimplex< -1, LinkType::Single, AccessScheme::Index, _Containment, _Allocator, Set< -1 > >: Set< -1 >
{
    enum { d = -1, };
};
/*
template< int _Dim,  template< class U, class V > class _Containment, template< class U > class _Allocator >
struct AbstractSimplex< _Dim, LinkType::Single, AccessScheme::Index, _Containment, _Allocator, TopologicalSpace< _Dim > >: TopologicalSpace < _Dim >
{
    enum {d = _Dim};
    ptrdiff_t upper, opponent, next;
    ptrdiff_t lower[_Dim + 1]; 
};

///terminate the recursion in the empty simplex (simplicial set) with dimension -1
template< template< class U, class V > class _Containment, template< class U > class _Allocator >
struct AbstractSimplex< -1, LinkType::Single, AccessScheme::Index, _Containment, _Allocator, TopologicalSpace< -1 > >: TopologicalSpace< -1 >
{
    enum { d = -1, };
};
*/

template< int _D, class _SC > struct ContainerFiller
{ 
    static inline void fill(_SC& s)
    {
        typedef AbstractSimplex< _D - 1, LinkType::Single, AccessScheme::Index, _SC::template _Containment, _SC::template _Allocator, typename _SC::template Space< _D - 1 > > AbstractSimplexT;
        typedef typename _SC::_Trait::template SimplexT< _D, AbstractSimplexT > Simplex;
        typedef typename _SC::template _Containment< Simplex, _SC::template _Allocator< Simplex > > Simplices;
        s.simplex_containers[_D - 1] = new Simplices;
        ContainerFiller< _D - 1, _SC >::fill(s);

    };	

};

template< class _SC > struct ContainerFiller< 0, _SC >{ static inline void fill(_SC&){}};

template< int _D, class _IT, class  _SC > struct IterFiller
{
    static inline void fill(_IT& i)
    {
        typedef  AbstractSimplex< _D - 1, LinkType::Single, AccessScheme::Index, _SC::template _Containment, _SC::template _Allocator, typename _SC::template Space < _D - 1 > > AbstractSimplexT;
        i.iterdata[_D] = new AbstractSimplexT;
        IterFiller< _D - 1, _IT, _SC >::fill(i);
    }
};

template< class _IT, class _SC > struct IterFiller< 0, _IT, _SC >{ static inline void fill(_IT&){}};

template< int _Dim, class _Trait,
    template< class U , class V > class _Containment , template< class U > class _Allocator, class _Space> 
    class AbstractSimplicialComplex< _Dim, _Trait, LinkType::Single, AccessScheme::Index, _Containment, _Allocator, _Space >
{
    public:
        enum { d = _Dim,};

        AbstractSimplicialComplex()
        {
            ContainerFiller< _Dim, AbstractSimplicialComplex >::fill(this);
        }
        ~AbstractSimplicialComplex()
        {
            for(int i = 0; i < _Dim; i++) delete simplex_containers[i];

        }
        void* simplex_containers[_Dim];
};
template< int _Dim, class _Trait, LinkType _LType, AccessScheme _AScheme, template< class U , class V > class _Containment , template< class U > class _Allocator, class _Space> 
class AbstractSimplicialComplexIterator< AbstractSimplicialComplex< _Dim, _Trait, _LType, _AScheme, _Containment, _Allocator, _Space > >
{
    public:
        typedef AbstractSimplicialComplex< _Dim, _Trait, _LType, _AScheme, _Containment, _Allocator, _Space > ASCT;
        typedef AbstractSimplicialComplexTopologyTrait< ASCT > ASCTopoT;
        typedef AbstractSimplicialComplexIteratorFunctionTrait< AbstractSimplicialComplexIterator > IterFuncT;

        AbstractSimplicialComplexIterator()
        {
            IterFiller< _Dim, AbstractSimplicialComplexIterator, ASCT >::fill(this);
        }

        ~AbstractSimplicialComplexIterator()
        {
            for(int i = 0; i < _Dim; i++) delete iterdata[i];
        }

        inline bool isvalid()
        {
            return IterFuncT::isvalid(*this);
        }
        template < int _D >
            inline AbstractSimplicialComplexIterator& simplexCCW()
            {
                ASCTopoT::template simplexCCV<_D, AbstractSimplicialComplexIterator >::doit(this);	
                return *this;
            }

        inline ptrdiff_t& operator [] (ptrdiff_t i)
        {
            return simplicesindices[i];
        }

        void*	    iterdata[_Dim];
        ptrdiff_t	simplicesindices[_Dim];
        ASCT	*m_sd;
};

template< int _Dim, class _Trait, template< class U , class V > class _Containment , template< class U > class _Allocator, class _Space> 
class AbstractSimplicialComplexTopologyTrait<  AbstractSimplicialComplex< _Dim, _Trait, LinkType::Single, AccessScheme::Index, _Containment, _Allocator, _Space > >
{
    public:

        typedef AbstractSimplicialComplex< _Dim, _Trait, LinkType::Single, AccessScheme::Index, _Containment, _Allocator, _Space > ASCT;
        typedef AbstractSimplicialComplexIterator< ASCT > IterT;
        //use const_expr for operations
        template< int _D, class _It >
            struct simplexFlip
            {
                inline bool doit(IterT &iter)
                {
                    bool succ;
                    ptrdiff_t opponent, parent, upper, max_dim;

                    opponent = iter.m_sd->simplex_containers[_D][iter.simplicesindices[_D]].opponent;
                    if(opponent == -1) return false;

                    iter.simplicesindices[_D] = opponent;

                    upper = iter.m_sd->simplex_containers[_D][iter.simplicesindices[_D]].upper;

                    for(int i = 0; i < _Dim - _D; i++)
                    {
                        iter.simplicesindices[_D + i + 1] = iter.m_sd->simplex_containers[_D + i][iter.simplicesindices[_D + i]].upper; 

                    }

                    ptrdiff_t n[_D];
                    simd_cp< int, _D >::eval(n, iter.m_sh->simplex_containers[_D][iter.simplicesindices[_D - 1]].vertices);

                    //... searching for odd permutation of vertices in the (_D ) (_D - 1) simplices
                    // it's enough to search for any permutation of the vertices, since two half simplices
                    // sharing the same vertices are the maximum - only two orientations (even and odd permutation)
                    // but the algebraic structure so the simple product addition is non-ambiguous
                    // (need a proof based on the eilenberg-zilber theorem and the kuenneth theorem)
                    //

                    simplexAlign< _D, _It>::doit(iter);
                    simplexFlip< _D - 2, _It >::doit(iter);
                    simplexCCW< _D - 1, _It >::doit(iter);
                    return succ;
                }
            };

        template< int _D, class _It >
            struct simplexAlign
            {
                inline bool doit(IterT &iter)
                {
                    bool succ;
                    int inv = simd_sum< int, _D >::eval(iter.m_sh->simplex_containers[_D - 1][iter.simplicesindices[_D - 1]].vertices);
                    IterT start, run;
                    start = iter;
                    iter.simplicesindices[ _D - 1] = iter.m_sd->simplex_containers[_D][iter.simplicesindices[_D]].lower[_D - 1];

                    do{

                        simplexCCW< _D - 1, _It >::doit(iter);

                    }while( simd_sum< int, _D >::eval(run.spimplicesindices[_D - 1].vertices) != inv);

                    simplexAlign< _D - 1, _It >::doit(iter);

                    return succ;

                }
            };

        template< class _It >
            struct simplexAlign< 0, _It >
            {
                inline bool doit(IterT &iter)
                {
                    return true;
                }
            };

        template< int _D, class _It >
            struct simplexCCW
            {
                inline bool doit(IterT &iter)
                {
                    bool succ;
                    simplexFlip< _D - 2, _It >::doit(iter);
                    simplexCCW< _D - 1, _It >::doit(iter);
                    return succ;
                }
            };

        template <class _It>
            struct simplexFlip< 1, _It>
            {
                inline bool doit(IterT &iter)
                {	
                    ptrdiff_t oppo, high;

                    return false;
                }

            };

        template <class _It>
            struct simplexCCW< 1, _It>
            {
                inline bool doit(IterT &iter)
                {
                    return false;
                }

            };

        template <class _It>
            struct simplexFlip< 0, _It>
            {
                inline bool doit(IterT &iter)
                {
                    return false;
                }
            };

        template <class _It>
            struct simplexCCW< 0, _It>
            {
                inline bool doit(IterT &iter)
                {
                    return false;
                }

            };

};

//template< template< class U, class V > class _Containment, template< class U > class _Allocator >

//test concretisations
#include <vector>

typedef AbstractSimplex< 0, LinkType::Single, AccessScheme::Index, std::vector, std::allocator, TopologicalSpace< 0 > > Simplex0;
typedef AbstractSimplex< 1, LinkType::Single, AccessScheme::Index, std::vector, std::allocator, TopologicalSpace< 1 > > Simplex1;
typedef AbstractSimplex< 2, LinkType::Single, AccessScheme::Index, std::vector, std::allocator, TopologicalSpace< 2 > > Simplex2;
typedef AbstractSimplex< 3, LinkType::Single, AccessScheme::Index, std::vector, std::allocator, TopologicalSpace< 3 > > Simplex3;
typedef AbstractSimplex< 4, LinkType::Single, AccessScheme::Index, std::vector, std::allocator, TopologicalSpace< 4 > > Simplex4;
typedef AbstractSimplex< 5, LinkType::Single, AccessScheme::Index, std::vector, std::allocator, TopologicalSpace< 5 > > Simplex5;
typedef AbstractSimplex< 6, LinkType::Single, AccessScheme::Index, std::vector, std::allocator, TopologicalSpace< 6 > > Simplex6;


//simple implementation of a point
//kinda n-vector stuff



//todo AS <-> multivectors/pseudoscalar (wedge product & other clifford algebra stuff) half simplices -> SO(n)
//wedge product, wedge sum ... bouquet of circles

#endif
