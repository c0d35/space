/**
 *
 * coder@computer.org
 * coder@0xc0d3.org
 *
 * this is only a sketch for some foo i'm thinking about
 */

//abstract simplex with topological space as default space
//metrical traits for every space (quasi-, pseudo-, semi- metrics 
//included divergences & stuff) bregman divergences as an example
//empty sets are void :3

#ifndef LE_SPACE_H
#define LE_SPACE_H

#include <cstddef>
#include <cstdint>
#include <iostream>
#include <array>
#include <cmath>
#include <memory>
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

template< int N > struct ifact{ enum { eval = N * ifact< N - 1>::eval }; };
template <> struct ifact< 0 > { enum { eval = 1 }; };
template< int N, int K >
struct nchoosek{ 
    enum { 
        eval = ifact< N >::eval / ( ifact< K >::eval *  ifact< N - K >::eval)
    };
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
/*
template< typename T >
struct fifo< T >
{
    T
}
*/
//quick hack integer chooser stuff
//
template< int S > struct SimpleLargerInteger{
	unsigned long long values[S/64];
};

template< int S > struct IntegerChooser
{
	typedef typename 
        IF< S <= sizeof(uint8_t), uint8_t,
		IF< S <= sizeof(uint16_t), uint16_t,
		IF< S <= sizeof(uint32_t), uint32_t,
		IF< S <= sizeof(uint64_t), uint64_t,
        SimpleLargerInteger< S > > > > >::RET IntegerType;
};

template< > struct IntegerChooser< 3 >
{
    enum { size = 3 };
    typedef uint32_t IntegerType;
};

template< > struct IntegerChooser< 4 >
{
	enum { size = 4};
	typedef uint32_t IntegerType;
};


template< > struct IntegerChooser< 6 >
{
	enum { size = 6};
    typedef __attribute__ ((aligned (8))) uint64_t IntegerType;
};

template< > struct IntegerChooser< 8 >
{
	enum { size = 8};
	typedef __attribute__ ((aligned (8))) uint64_t IntegerType;
};



enum class ArchType: short
{
    GENERIC,
    AMD64,
    AARCH64,
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

template < int K, int D, class S>
struct ExteriorPower;
template < int D, class S >
struct ExteriorAlgebra;
template< int N > struct Algebra {};
template< int N > struct GradedAlgebra {};
template< int K, int N > struct QuotientSpace {};
template< int _Dim >
struct Set { enum { d = _Dim }; };
template< int _Dim >
struct TopologicalSpace: Set< _Dim > {};
template< int _Dim >
struct UniformSpace: TopologicalSpace< _Dim > {};
template< int _Dim, template < int __D > class _M  >
struct MetricSpace: UniformSpace < _Dim > {};

/**
 * kinda abstract simplex(?), one example for an element
 */
template< int _Dim, LinkType _LType, AccessScheme _AScheme,
    template< class U , class V > class Containment , template< class U >
    class Allocator, class _Space > class AbstractHalfSimplicialComplex{};
template< typename _ASD > struct AbstractHalfSimplicialComplexTopologyTrait{};
template< typename _ASD > struct
AbstractHalfSimplicialComplexTopologicalOperations{};
template< typename _Iterator > 
struct AbstractHalfSimplicialComplexIteratorFunctionTrait{};
template< typename _ASD > struct AbstractHalfSimplicialComplexIterator{};
template< int _Dim, LinkType _LType, AccessScheme _AScheme,  
    template< class U, class V > class _Containment,
    template< class U > class _Allocator, class _Space >
    struct AbstractHalfSimplex: _Space {};


///why d + 1 you may ask. 'cause \f$ {d+1 \choose d} = d+1\f$ , u kno'. 
//thus we're terminating the recursion at dimension = -1
template< int _Dim,  template< class U, class V > class _Containment,
    template< class U > class _Allocator,  class _Space >
    struct AbstractHalfSimplex< _Dim, LinkType::Single, AccessScheme::Index,
    _Containment, _Allocator, _Space >: _Space
{
    enum {d = _Dim};
    ptrdiff_t upper, opponent,
              //opponext,
              next;
    //ptrdiff_t lower[_Dim + 1]; 
    ptrdiff_t lower;
    AbstractHalfSimplex()
    {
        upper = -1;
        opponent = -1;
        next = -1;
        lower = -1;
    };
    AbstractHalfSimplex(const AbstractHalfSimplex &s): _Space(s)
    {
        upper = s.upper;
        opponent = s.opponent;
        //opponext = s.opponext;
        next = s.next;
        lower = s.lower;
    }
};

///terminate the recursion in the empty simplex (simplicial set) with
//dimension -1
template< template< class U, class V > class _Containment, 
    template< class U > class _Allocator, 
    class _Space >
    struct AbstractHalfSimplex< -1, LinkType::Single,
    AccessScheme::Index, _Containment,
    _Allocator, _Space >: _Space
{
    enum { d = -1, };
    AbstractHalfSimplex(){};
    AbstractHalfSimplex(AbstractHalfSimplex &s){};
};

///why d + 1 you may ask. 'cause \f$ {d+1 \choose d} = d+1\f$ , u kno'. thus
//we're terminating the recursion at dimension = -1
template< int _Dim,  template< class U, class V > class _Containment,
    template< class U > class _Allocator >
    struct AbstractHalfSimplex< _Dim, LinkType::Single,
    AccessScheme::Index, _Containment, _Allocator, Set< _Dim > >: Set < _Dim >
{
    enum {d = _Dim};
    ptrdiff_t upper, opponent, 
              //opponext,
              next;
    //ptrdiff_t lower[_Dim + 1];
    ptrdiff_t lower; 
    AbstractHalfSimplex(){};
    AbstractHalfSimplex(AbstractHalfSimplex &s){};
};

///terminate the recursion in the empty simplex (simplicial set) with dimension -1
template< template< class U, class V > class _Containment,
    template< class U > class _Allocator >
struct AbstractHalfSimplex< -1, LinkType::Single, AccessScheme::Index, _Containment,
    _Allocator, Set< -1 > >: Set< -1 >
{
    enum { d = -1, };
    AbstractHalfSimplex(){};
    AbstractHalfSimplex(AbstractHalfSimplex &s){};
};

template< int D, class SC > struct makeHalfSimplex
{
    static inline void make(SC &s)
    {
        
    }

};
//default filling with complete simplex of dimension D
template< int D, class _SC > struct ContainerFiller
{ 
    static inline void fill(_SC& s)
    {
        /*
        typedef AbstractHalfSimplex< _D - 1, LinkType::Single, 
                AccessScheme::Index, _SC::template Containment,
                _SC::template Allocator,
                typename _SC::template Space< _D - 1 > > AbstractHalfSimplexT;
        typedef typename _SC::_Trait::template 
            HalfSimplexT< _D, AbstractHalfSimplexT > HalfSimplex;
        typedef typename _SC::template 
            _Containment< HalfSimplex, 
            _SC::template _Allocator< HalfSimplex > > Simplices;
        s.simplex_containers[_D - 1] = new Simplices;
        */
        /*typedef AbstractHalfSimplex< _D, LinkType::Single, 
                AccessScheme::Index, _SC::template Containment,
                _SC::template Allocator, typename _SC::template Space< _D > > HalfSimplex;*/
        s.simplex_containers[D] = new typename _SC::template Containment<
            typename _SC::template HalfSimplex< D >, typename _SC::template Allocator< 
            typename _SC::template HalfSimplex< D > > >;
        ContainerFiller< D - 1, _SC >::fill(s);
    };	
};

template< class _SC > 
struct ContainerFiller< -1, _SC >{ static inline void fill(_SC&){}};
//template< class _SC > 
//struct ContainerFiller< 0, _SC >{ static inline void fill(_SC&){}};


template< int _D, class _IT, class  _SC > struct IterFiller
{
    static inline void fill(_IT* i)
    {
        typedef  AbstractHalfSimplex< _D, LinkType::Single,
                 AccessScheme::Index, _SC::template Containment,
                 _SC::template Allocator, 
                 typename _SC::Space > HalfSimplex;
        i->iterdata[_D ] = new HalfSimplex;
        //void* foo = new HalfSimplex;
        IterFiller< _D - 1, _IT, _SC >::fill(i);
    }
};

template< class _IT, class _SC > struct IterFiller< -1, _IT, _SC >{
    static inline void fill(_IT* i){}
};

template< int _Dim,
    template< class U , class V > class _Containment,
    template< class U > class _Allocator,  class _Space> 
    class AbstractHalfSimplicialComplex< _Dim,
    LinkType::Single, AccessScheme::Index,
    _Containment, _Allocator, _Space >
{
    public:
        enum { d = _Dim,};
        template< class U, class V >
            using Containment = _Containment< U, V >;
        template< class U >
            using Allocator = _Allocator< U >;
        //template< int D >
        using Space = _Space;
        template< int D >
          using HalfSimplex = AbstractHalfSimplex< D, LinkType::Single,
                 AccessScheme::Index, _Containment,
                 _Allocator, Space >;
        template< int D >
            using Container = Containment< HalfSimplex< D >,
                  Allocator< HalfSimplex< D > > >;

  
        AbstractHalfSimplicialComplex()
        {
            ContainerFiller< _Dim, AbstractHalfSimplicialComplex >::fill(*this);
        }
        ~AbstractHalfSimplicialComplex()
        {
            //for(int i = 0; i < _Dim; i++) delete simplex_containers[i];

        }
        
        /*
        template < int D >
        static  Container< D >& getContainer(int i = D)
            {
                return *((Container< D > *)simplex_containers[D]);
            }
         */
        inline void* foo(ptrdiff_t i)
        {
            return simplex_containers[i];
        }


        template < int D >
        inline AbstractHalfSimplicialComplexIterator < AbstractHalfSimplicialComplex >
        insert( HalfSimplex< D > &s )
        {

            AbstractHalfSimplicialComplexIterator<AbstractHalfSimplicialComplex>
                iter(this);
            iter.insert(s);
            return iter;

        }
        /*
        template< int i >
        constexpr  Containment< HalfSimplex< i >,
                  Allocator< HalfSimplex < i > > > *
                      fick( int j = i )
                      {
                          return static_cast< 
                              Containment< HalfSimplex< i >,
                  Allocator< HalfSimplex < i > > > *
>(simplex_containers[i]);
                      }
                      */

        inline void* operator [] (const int i)
        {
            //Container< i > C;
            return simplex_containers[i];
            //return *((Container< i > *) simplex_containers[D]);
            //return new F;
        }
  
        void* simplex_containers[_Dim + 1];
        int bla;
};

//is an HalfSimplex< 1 > homoeomorph to a proper minimal iterator valid
//on every simplicial complex?
//
template< int _Dim, LinkType _LType, AccessScheme _AScheme,
    template< class U , class V > class _Containment,
    template< class U > class _Allocator, class _Space> 
class AbstractHalfSimplicialComplexIterator<
AbstractHalfSimplicialComplex< _Dim, _LType,
    _AScheme, _Containment, _Allocator, _Space > >
{
    public:
        typedef AbstractHalfSimplicialComplex< _Dim, _LType,
                _AScheme, _Containment, _Allocator, _Space > ASCT;
        typedef AbstractHalfSimplicialComplexTopologyTrait< ASCT > ASCTopoT;
        typedef AbstractHalfSimplicialComplexTopologicalOperations < ASCT >
            ASCTopoOp;
        typedef AbstractHalfSimplicialComplexIteratorFunctionTrait<
            AbstractHalfSimplicialComplexIterator > IterFuncT;
        template< int D >
          using HalfSimplex = AbstractHalfSimplex< D, LinkType::Single,
                 AccessScheme::Index, _Containment,
                 _Allocator, _Space >;
        template< int D >
        using Container = _Containment< HalfSimplex< D >,
              _Allocator< HalfSimplex < D > > >;


        AbstractHalfSimplicialComplexIterator(ASCT* sc)
        {
            m_sd = sc;
            IterFiller< _Dim, AbstractHalfSimplicialComplexIterator,
                ASCT >::fill(this);
        }
        AbstractHalfSimplicialComplexIterator()
        {
            IterFiller< _Dim, AbstractHalfSimplicialComplexIterator,
                ASCT >::fill(this);
        }
        AbstractHalfSimplicialComplexIterator(const 
                AbstractHalfSimplicialComplexIterator &iter)
        {
            for(int i = 0; i < _Dim + 1; i++)
            {
                iterdata[i] = iter.iterdata[i];
                simplicesindices[i] = iter.simplicesindices[i];
                m_sd = iter.m_sd;
            }

        }
   
        ~AbstractHalfSimplicialComplexIterator()
        {
            //for(int i = 0; i < _Dim; i++) delete iterdata[i];
        }

        inline bool isvalid()
        {
            return IterFuncT::isvalid(*this);
        }
        inline AbstractHalfSimplicialComplexIterator&
            make(HalfSimplex< 0 > s[ _Dim + 1])
            {
                return ASCTopoOp::template makeAbstractSimplex<
                       AbstractHalfSimplicialComplexIterator, _Dim >
                       ::doit(*this, s);
            }
        template < int D >
            inline AbstractHalfSimplicialComplexIterator& 
            insert(HalfSimplex< D > &s)
            {
                ASCTopoT::template insert<D,
                AbstractHalfSimplicialComplexIterator >::doit(*this, s);	
                return *this;
            }
        template < int _D >
            inline AbstractHalfSimplicialComplexIterator& simplexCCW()
            {
                ASCTopoT::template simplexCCV<_D,
                    AbstractHalfSimplicialComplexIterator >::doit(*this);	
                return *this;
            }
        template < int D >
            inline Container< D > &
            getSimplices(HalfSimplex< D > &hs)
            {
                return *(static_cast< Container< D > * >
                        (m_sd->simplex_containers[D]));

            }
        template < int D >
            inline HalfSimplex< D > &
            get(HalfSimplex< D > &hs)
            {

                return 
                (*(static_cast< Container< D > * >
                        (m_sd->simplex_containers[D])))
                    [simplicesindices[ D ]];


            }
        inline ptrdiff_t getID(ptrdiff_t i)
        {
            return simplicesindices[i];
        }

        inline ptrdiff_t& operator [] (ptrdiff_t i)
        {
            return simplicesindices[i];
        }

        template < int D >
            inline HalfSimplex < D > &
            makeHalfSimplex()
            {
                /*
                      for(int i = 0; i < nchoosek< ASCT::d, D >::eval; i++)
          (static_cast< Container< _D > * >(
                                            s.simplex_containers[D]
                                           ).push_back(
                                               new typename 
                                               _SC::template HalfSimplex< D >
                                               );
*/
            }
        void*	    iterdata[_Dim + 1];
        ptrdiff_t	simplicesindices[_Dim + 1];
        ASCT	*m_sd;
};

template< int _Dim,
    template< class U , class V > class _Containment,
    template< class U > class _Allocator, class _Space> 
class AbstractHalfSimplicialComplexTopologyTrait<  
AbstractHalfSimplicialComplex< _Dim, LinkType::Single,
    AccessScheme::Index, _Containment, _Allocator, _Space > >
{
    public:

        typedef AbstractHalfSimplicialComplex< _Dim,
                LinkType::Single, AccessScheme::Index, _Containment,
                _Allocator, _Space > ASCT;
        typedef AbstractHalfSimplicialComplexIterator< ASCT > IterT;

        template< int D >
          using HalfSimplex = AbstractHalfSimplex< D, LinkType::Single,
                AccessScheme::Index, _Containment,
                _Allocator, _Space >;

        template< int D >
        using Container = _Containment< HalfSimplex< D >,
              _Allocator< HalfSimplex < D > > >;

        //use const_expr for operations
/*
        template< int D, class _It >
            struct makeSimplex
            {
                static inline bool doit(_It &iter)
                {
                    bool succ = true;
                    return succ;
                }
            };
        template< class _It >
            struct makeSimplex< 0, _It >
            {
                static inline bool doit(_It &iter)
                {
                    return true;
                }
            };
            */
        template< int _D, class _It >
            struct simplexFlip
            {
                static inline bool doit(IterT &iter)
                {
                    bool succ;
                    ptrdiff_t opponent, parent, upper, max_dim;

                    opponent = iter.m_sd->simplex_containers[_D]
                        [iter.simplicesindices[_D]].opponent;
                    if(opponent == -1) return false;

                    iter.simplicesindices[_D] = opponent;

                    upper = iter.m_sd->simplex_containers[_D]
                        [iter.simplicesindices[_D]].upper;

                    for(int i = 0; i < _Dim - _D; i++)
                    {
                        iter.simplicesindices[_D + i + 1] =
                            iter.m_sd->simplex_containers[_D + i]
                            [iter.simplicesindices[_D + i]].upper; 

                    }

                    ptrdiff_t n[_D];
                    simd_cp< int, _D >::eval(n,
                            iter.m_sh->simplex_containers[_D]
                            [iter.simplicesindices[_D - 1]].vertices);

                    //... searching for odd permutation of vertices in the
                    //(_D ) (_D - 1) simplices
                    // it's enough to search for any permutation of the
                    // vertices, since two half simplices
                    // sharing the same vertices are the maximum - only two
                    // orientations (even and odd permutation)
                    // but the algebraic structure so the simple product
                    // addition is non-ambiguous
                    // (need a proof based on the eilenberg-zilber theorem and
                    // the kuenneth theorem)
                    //

                    simplexAlign< _D, _It>::doit(iter);
                    simplexFlip< _D - 2, _It >::doit(iter);
                    simplexCCW< _D - 1, _It >::doit(iter);
                    return succ;
                }
            };
        template< int _D, class _It >
            struct insert
            {
                static inline bool doit(IterT &iter, 
                        typename ASCT::template HalfSimplex< _D > s)
                {
                    (static_cast< Container< _D > * >(
                                                      iter.m_sd
                                                      ->
                                                      simplex_containers[_D]
                                                     ))
                        ->push_back(s);
                    iter.simplicesindices[ _D ] = 
                        (static_cast< Container< _D > * >(
                                                          iter.m_sd
                                                          ->
                                                          simplex_containers[_D]
                                                         ))
                        ->size() - 1;


                  /* 
                    bool succ;
                    int inv = simd_sum< int,
                        _D >::eval((*(static_cast< Container< _D - 1 > * >(
                                    iter.m_sd->simplex_containers[_D - 1]))
                                [iter.simplicesindices[_D - 1]]).vertices);
                    IterT start, run;
                    start = iter;
                    iter.simplicesindices[ _D - 1] = 
                        iter.m_sd->simplex_containers[_D]
                        [iter.simplicesindices[_D]].lower[_D - 1];

                    do{

                        simplexCCW< _D - 1, _It >::doit(iter);

                    }while( simd_sum< int, _D >::eval(
                                run.spimplicesindices[_D - 1].vertices) != inv);
                    simplexAlign< _D - 1, _It >::doit(iter);
                    //
                    //(static_cast< Container< _D > * >(
                    //                                  iter.m_sd
                    //                                  ->
                    //                                  simplex_containers[_D]
                    //                                  ))
                    //   ->push_back(s);
                        
                    */
                    return true;

                }
            };


        template< int _D, class _It >
            struct simplexAlign
            {
                static inline bool doit(IterT &iter)
                {
                    bool succ;
                    int inv = simd_sum< int,
                        _D >::eval(iter.m_sh->simplex_containers[_D - 1]
                                [iter.simplicesindices[_D - 1]].vertices);
                    IterT start, run;
                    start = iter;
                    iter.simplicesindices[ _D - 1] = 
                        iter.m_sd->simplex_containers[_D]
                        [iter.simplicesindices[_D]].lower[_D - 1];

                    do{

                        simplexCCW< _D - 1, _It >::doit(iter);

                    }while( simd_sum< int, _D >::eval(
                                run.spimplicesindices[_D - 1].vertices) != inv);
                    simplexAlign< _D - 1, _It >::doit(iter);
                    return succ;

                }
            };

        template< class _It >
            struct insert< 0, _It >
            {
                static inline bool doit(IterT &iter, 
                        typename ASCT::template HalfSimplex< 0 > s)
                {
                    (static_cast< Container< 0 > * >(
                                                      iter.m_sd
                                                      ->
                                                      simplex_containers[0]
                                                      ))
                        ->push_back(s);

                    iter.simplicesindices[ 0 ] = 
                        (static_cast< Container< 0 > * >(
                                                          iter.m_sd
                                                          ->
                                                          simplex_containers[0]
                                                         ))
                        ->size() - 1;


                    return true;

                }

            };

        template< class _It >
            struct simplexAlign< 0, _It >
            {
                static inline bool doit(IterT &iter)
                {
                    return true;
                }
            };

        template< int _D, class _It >
            struct simplexCCW
            {
                static inline bool doit(IterT &iter)
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
                static inline bool doit(IterT &iter)
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

//generators should be independend
//make full simplex (?)
template< int _Dim,
    template< class U , class V > class _Containment,
    template< class U > class _Allocator, class _Space >
    struct      AbstractHalfSimplicialComplexTopologicalOperations<  
    AbstractHalfSimplicialComplex< _Dim, LinkType::Single,
    AccessScheme::Index, _Containment, _Allocator, _Space > >
{
    public:

        typedef AbstractHalfSimplicialComplex< _Dim,
                LinkType::Single, AccessScheme::Index, _Containment,
                _Allocator, _Space > ASCT;
        typedef AbstractHalfSimplicialComplexIterator< ASCT > Iter;

        template< int D >
            using HalfSimplex = AbstractHalfSimplex< D, LinkType::Single,
                  AccessScheme::Index, _Containment,
                  _Allocator, _Space >;

        template< int D >
            using Container = _Containment< HalfSimplex< D >,
                  _Allocator< HalfSimplex < D > > >;
//assuming that the HalfSimplex<0> are already inserted and part
//of the simplicial complex referenced by the iterator
        template< int D, class _It >
            struct appendAbstractSimplex
            {
                static inline _It& doit(_It &iter, HalfSimplex< 0 > s[_It::d])
                {
                    bool succ = true;
                    return succ;
                }
            };
        template< class _It >
            struct appendAbstractSimplex< 0, _It >
            {
                static inline _It doit(_It &iter, HalfSimplex< 0 > s[1])
                {
                    return true;
                }
            };
        template< class _It, int D >
            struct setOpponent
            {
                static inline void doit(_It iterators[D])
                {
                    HalfSimplex< D - 2 > append;

                    for(int i = 0; i < (D + 1); i++)
                    {
                        iterators[i].get(append).opponent = 
                            iterators[(i + 1) % (D + 1)][D - 2];
                        iterators[(i + 1) % (D + 1)].get(append). opponent =
                            iterators[i][D - 2];
                    }
                }
            };
        template< class _It >
            struct setOpponent< _It, 0 >
            {
                static inline void doit(_It* iterators)
                {

                }
            };
        template< class _It >
            struct setOpponent< _It, 1 >
            {
                static inline void doit(_It* iterators)
                {

                }
            };

        template< class _It, int D >
            struct setNext
            {
                static inline void doit(_It iterators[D + 1])
                {
                    HalfSimplex< D - 1 > append;

                    for(int i = 0; i < (D + 1); i++)
                    {
                        iterators[i].get(append).next = 
                            iterators[(i + 1) % (D + 1)][D - 1];
                    }
                }
            };
        template< class _It >
            struct setNext< _It, 0 >
            {
                static inline void doit(_It* iterators)
                {

                }
            };
        template< class _It, int D >
            struct setUpperLower
            {
                static inline void doit(_It iterators[D + 1])
                {
                    HalfSimplex< D > upper;
                    HalfSimplex< D - 1 > lower;
                    /*
                                 (*(static_cast< Container< D > * >
                        (m_sd->simplex_containers[D])))
                    [simplicesindices[ D ]];*/



                    for(int i = 0; i < (D); i++)
                    {
                        Container< D > *cp = static_cast< Container< D > *>(
                                iterators[i].m_sd->simplex_containers[D]);
                        iterators[i].get(upper).lower = 
                            iterators[i][D - 1];

                        iterators[i].get(lower).upper =
                            iterators[i][D]; 
 
                    }
                }
            };

        template< class _It >
            struct setUpperLower< _It, -1 >
            {
                static inline void doit(_It iterators[0])
                {
                    /*
                    HalfSimplex< D > upper;
              //      HalfSimplex< D - 1 > lower;

                    for(int i = 0; i < (1 + 1); i++)
                    {
                        iterators[i].get(upper).lower = 
                            iterators[i][1 - 1];
                        iterators[i].get(lower).upper =
                            iterators[i][1];

                    }
                    */
                    //nothin to do, cause upper is set by 
                    //setUpperLower<1> 
                }
            };


        template< class _It , int D >
            struct makeAbstractSimplex
            {
                static inline Iter& doit(Iter &it, HalfSimplex< 0 > s[D + 1])
                {
                    //directed complete partial order
                    bool succ = true;
                    HalfSimplex< D > hs;
                    it.insert(hs);
                    //it.make(s);
                    Iter iterators[D + 1];
                    for(int i = 0; i < (D + 1); i++) iterators[i] = it;
                    HalfSimplex< 0 > sub[D];
                    for(int i = 0; i < (D + 1); i++)
                    {
                        for(int j = 0; j < D; j++)
                        {
                            sub[j] = s[(j + i) % D];
                        }
                        if(i & 1)std::swap(sub[0], sub[D - 1]);//loop separation
                        //it.make(
                        iterators[i] = makeAbstractSimplex< _It, D - 1 >
                        ::doit(it, sub); //iterator now contains handles
                    }
                      setOpponent< _It, D >::doit(iterators);
                    setNext< _It, D >::doit(iterators);
                    setUpperLower< _It, D >::doit(iterators);
                    /*
                    if((D - 2) > 0)
                    {
                        HalfSimplex< D - 2 > append;

                        for(int i = 0; i < (D + 1); i++)
                        {
                            iterators[i].get(append).opponent = 
                                iterators[(i + 1) % (D + 1)][D - 2];
                            iterators[(i + 1) % (D + 1)].get(append). opponent =
                                iterators[i][D - 2];
                        }
                    }
                    */
                    return it;
                }
            };
       template< class _It >
            struct makeAbstractSimplex< _It, 0 >
            {
                static inline _It doit(_It &it, HalfSimplex< 0 > s[1])
                {
                    HalfSimplex< 0 > hs;
                    it.insert(hs);
                    return it;
                }
            };
       template< class It, int D >
           struct getHamiltonianPath
           {
               static inline Iter& doit(It &it, std::unique_ptr< HalfSimplex< D > > circle)
               { 
                   return it;

               }

           };
       template< class It >
           struct getHamiltonPathVertices
           {
               static inline Iter& doit(It &it, std::unique_ptr< HalfSimplex< 0 > > circle)
               {
                   return it;
               }
           };

};

#include <vector>
#include <map>

//specialise the simplices for proper spaces
template< int D > using HalfSimplex = AbstractHalfSimplex< D, LinkType::Single,
    AccessScheme::Index, std::vector, std::allocator, TopologicalSpace< D + 1 > >;
template< int D > using HalfSimplicialComplex = AbstractHalfSimplicialComplex< D,
    LinkType::Single, AccessScheme::Index, std::vector, std::allocator,
    TopologicalSpace< D + 1 > >;
template< int D > using HalfSimplicialComplexTopologyTrait =
AbstractHalfSimplicialComplexTopologyTrait < HalfSimplicialComplex< D + 1 > >;
template< int D > using HalfSimplicialComplexIterator = AbstractHalfSimplicialComplexIterator< HalfSimplicialComplex< D > >;




//simple implementation of a point
//kinda n-vector stuff if used in a vector-space
#include <iostream>
template< int D, typename _Type, typename ... _Args >
struct SimplePoint
{
	public:
		enum { d = D };
		typedef _Type Type;
        typedef Type ValueType;
		_Type values[D];
		inline _Type& operator [](ptrdiff_t n) { return values[n];}

		SimplePoint()
		{
			for(ptrdiff_t i = 0; i < D; i++) values[i] = 0;
		}
		SimplePoint( std::initializer_list< Type > val)
		{
			ptrdiff_t i = 0;
			for(auto v : val){ values[i++] = v;}
		}
		SimplePoint(const SimplePoint& p)
		{
			for(ptrdiff_t i = 0; i < D; i++)
			{
				values[i] = p.values[i];
			}
		}
		inline SimplePoint& operator = (const SimplePoint& p)
		{
			for(ptrdiff_t i = 0; i < D; i++)
			{
				values[i] = p.values[i];
			}
			return *this;
		}
        inline SimplePoint& operator += (const SimplePoint& p)
        {
            for(ptrdiff_t i = 0; i < D; i++)
            {
                values[i] += p.values[i];
            }
            return *this;
        }
        inline SimplePoint& operator -= (const SimplePoint& p)
        {
            for(ptrdiff_t i = 0; i < D; i++)
            {
                values[i] -= p.values[i];
            }
            return *this;
        }
        inline const SimplePoint operator + (const SimplePoint& p) const
        {
            return SimplePoint(*this) += p;
        }
        inline const SimplePoint operator - (const SimplePoint& p) const
        {
            return SimplePoint(*this) -= p;
        }
};

template< int D, typename CompType, typename MachineType >
struct ManhattanPoint
{
	public:
		enum { d = D };
		typedef CompType Type;
        typedef Type ValueType; //machine type
        //typedef Type MachineType;
		Type values[D];
		inline Type& operator [](ptrdiff_t n) { return values[n];}
		ManhattanPoint()
		{
			for(ptrdiff_t i = 0; i < D; i++) values[i] = 0;
		}
		ManhattanPoint( std::initializer_list< Type > val)
		{
			ptrdiff_t i = 0;
			for(auto v : val) values[i++] = v;
		}
		ManhattanPoint(const ManhattanPoint& p)
		{
			for(ptrdiff_t i = 0; i < D; i++)
			{
				values[i] = p.values[i];
			}
		}
		inline ManhattanPoint& operator = (const ManhattanPoint& p)
		{
			for(ptrdiff_t i = 0; i < D; i++)
			{
				values[i] = p.values[i];
			}
			return *this;
		}
};


template < int D, template< int _D > class _M >
struct LinearSpace: public MetricSpace< D, _M >
{
    enum{ d = D };

    typedef _M< D > Metric;
    typedef typename Metric::KeyType KeyType;
    typedef typename Metric::ValueType ValueType;
    //typedef double ValueType;
    //should i pull the Tensor Definitions from the exterior algebra?
    //ok, let's do it :|
    typedef typename ExteriorPower< 0, 0, LinearSpace>::Vector Scalar;
    /*
    struct Scalar: HalfSimplex< 0 >
    {
        enum { k = 0, d = D };

    };*/
    struct Vertex: HalfSimplex< 0 >//derive from simplex, cause every
                   //linear space is hausdorffian (kolomogorov T_2),
                   //so i hope that's okay
    {
        enum { k = 1, d = D};
        Scalar v[D];
        inline Scalar& operator [](ptrdiff_t n) { return v[n];}
        Vertex()
        {
            //for(ptrdiff_t i = 0; i < D; i++) v[i] = 0;
        }
        Vertex( std::initializer_list< Scalar > val)
        {
            ptrdiff_t i = 0;
          //  for(auto v : val) v[i++] = v;
        }
        Vertex(const Vertex& p)
        {
            for(ptrdiff_t i = 0; i < D; i++)
            {
                v[i] = p.v[i];
            }
        }
        inline Vertex& operator = (const Vertex& p)
        {
            for(ptrdiff_t i = 0; i < D; i++)
            {
                v[i] = p.v[i];
            }
            return *this;
        }
    };

    typedef typename ExteriorPower< 1, d, LinearSpace >::Vector Vector;
    //typedef Vertex Vector; //i guess, every vector is a 1-cell, so ...

    HalfSimplicialComplex< D > simplicial_decomposition;
    Vector e[D]; //basis

};

//metrics and norms
//
template< int D, typename _Type, template < int _D > class _Point >
struct EuklidianMetric
{
    enum{ d = D, };
    typedef _Type ValueType;
    typedef _Point< d > Point;
    typedef typename _Point< d >::ValueType PValueType;
    template < int N > using ElementT = _Point < N >;
    typedef typename IntegerChooser< D * sizeof(ValueType) >::IntegerType 
        KeyType;
    typedef Point Vector;
    //struct Vector: Point {};
    static inline KeyType dist(Point p0, Point p1)
    {
        KeyType dist = 0;
        for(ptrdiff_t i = 0; i < D; i++)dist += 
            p1[i] - p0[i] *
            p1[i] - p0[i];
        dist = std::sqrt(dist);
        return dist;
    }
    static inline Vector distance(Point p0, Point p1)
    {
        Vector d;
        for(ptrdiff_t i = 0; i < D; i++)d[i] = p1[i] - p0[i];
        return d;
    }
    /*
    inline Vector operator - (Vector v0)
    {
        Vector v;
        for(ptrdiff_t i = 0; i < D; i++) v[i] = v0[i]
    }
    */
    //draftly putting Morton Code into the metric
    //only makes sense with integer

    static inline KeyType morton_encode(Point p)
	{
		KeyType pre_key = 0;
		KeyType mask = 1;
        for(int i = 1; i <= (sizeof(ValueType) * 8); i++)
		{

            mask = 1;
            for(int j = 0; j < d; j++)
            {
                pre_key <<= 1;
                pre_key |= (p[j] >> ((sizeof(ValueType) * 8)- i)) & (mask);
            }
        }
		return pre_key;
	}
    
    static inline Point morton_decode(KeyType k)
    {
        KeyType pre_key = k;
        KeyType mask = 1;
        Point p;
        for(int i = 0; i < (sizeof(KeyType) * 8); i++)
        {
                p[(i) % d] <<= 1;
                p[(i) % d] |= (pre_key >>
                        ((sizeof(KeyType) * 8) - 1))  & mask;
                //p[(i + d - 1) % d] <<= 1;
                //p[(i + d - 1) % d] |= (pre_key >>
                //        ((sizeof(KeyType) * 8) - 1))  & mask;
                pre_key <<= 1;
        }
        return p;
    }

};



template< int D, typename _Type, template < int _D > class _E >
struct ManhattanMetric
{
    enum{ d = D, };
    typedef _Type ValueType;
    typedef _E< d > Point;
    template < int N > using ElementT = _E< N >;
    typedef typename IntegerChooser< D * sizeof(ValueType) >::IntegerType
        KeyType;

};

//hm'kay currently no how to correctly classify
template< int D, typename MachineType >
struct MortonMetric
{
    enum{ d = D, };
    typedef MachineType Element;
   // template < int N > using ElementT = _E< N >;
  //  typedef typename IntegerChooser< D * sizeof(ValueType) >::IntegerType
   //     KeyType;
};

//to be used in the constructors of spaces
template< class X, class Y >
struct MorphismTrait
{
    static inline void morph(X x, Y y)
    {
    }
};


template < ArchType a, int D, typename Type > struct EuklidianMetricTrait { };

template< int D > using SimplePointFdouble = SimplePoint< D, double >;
template< int D > using SimplePointFuint16 = SimplePoint< D, uint16_t >;
template< int D > using SimplePointFuint32 = SimplePoint< D, uint32_t >;
template< int D > using SimplePointFuint64 = SimplePoint< D, uint64_t >;
template< int D > using SimpleEuklidianMetricFdouble = EuklidianMetric< D,
    double,
    SimplePointFdouble >;
template< int D > using SimpleEuklidianMetricFuint16 = EuklidianMetric< D,
    int16_t, 
    SimplePointFuint16 >;
template< int D > using SimpleEuklidianMetricFuint32 = EuklidianMetric< D,
    int32_t, 
    SimplePointFuint16 >;
template< int D > using SimpleEuklidianMetricFuint64 = EuklidianMetric< D,
    int64_t, 
    SimplePointFuint64 >;

//template< int D >
//using EuklidianSpace = LinearSpace<D, EuklidianMetric>;
template< int D >
using SimpleEuklidianSpaceFuint16 = LinearSpace< D,
      SimpleEuklidianMetricFuint16 >;
template< int D >
using SimpleEuklidianSpaceFuint32 = LinearSpace< D,
      SimpleEuklidianMetricFuint32 >;
template< int D >
using SimpleEuklidianSpaceFuint64 = LinearSpace< D,
      SimpleEuklidianMetricFuint64 >;
template< int D >
using SimpleEuklidianSpaceFdouble = LinearSpace< D,
      SimpleEuklidianMetricFdouble >;

//specialise the simplices for euklidian spaces
/*
template< int D, class M > using HalfSimplex = AbstractHalfSimplex< D, LinkType::Single,
    AccessScheme::Index, std::vector, std::allocator,
    EuklidianSpace< D, M > >;
template< int D, class M > using HalfSimplicialComplex =
AbstractHalfSimplicialComplex< D, LinkType::Single, AccessScheme::Index,
    std::vector, std::allocator, EuklidianSpace< D, M > >;
template< int D, class M > using HalfSimplicialComplexTopologyTrait =
AbstractHalfSimplicialComplexTopologyTrait < HalfSimplicialComplex< D, M > >;
template< int D, class M > using HalfSimplicialComplexIterator = AbstractHalfSimplicialComplexIterator< HalfSimplicialComplex< D, M > >;
*/


template< int K, int D, class S > struct ExteriorPower: QuotientSpace< K, D >
{
    enum{ d = nchoosek< D, K >::eval };
    typedef typename S::Metric Metric;
    struct Vector;

    typedef typename ExteriorPower< 0, 0, S >::Vector Scalar;
    //typedef typename ExteriorPower< 0, 0, S >::Vector Scalar;

    struct Vector: HalfSimplex< K >
    {
        typename Metric::template ElementT< d > v;
        //or should i use an vector of scalars (would be more canonical)
        //todo:
        //- check performance penalty
        Vector(){};
        inline typename Metric::ValueType operator [](ptrdiff_t n)
        {
            return v[n];
        }
    };
    struct Null: Vector
    {
        Null(){};//nullify
    };

    //static constexpr Vector e[d] = {}; // basis // need to be calced by S::e
    inline Vector operator ^ (const Vector& v) { return new Vector; }//wedge
    inline typename ExteriorPower< D - K, D, S>::Vector operator * ()
    { return  ExteriorPower< D - K, D, S>::Vector; } //Hodge-*
    //\f$ \star : \bigwedge^{k} V \to \bigwedge^{n-k} V \choosekVector e[d];
    //base

};

template< int D, class S >
struct ExteriorAlgebra: GradedAlgebra< D >
{
    enum{ d = ipow< 2, D >::eval };
};

/**
 * begin LP-Tree Stuff
 *
 */
/**
 * Hyper Cube foo
 */

template< int _Dim, LinkType _LType, AccessScheme _AScheme,  
    template< class U, class V > class _Containment,
    template< class U > class _Allocator, class _Space >
    struct AbstractHyperCube: _Space {};

template< int _Dim,  template< class U, class V > class _Containment,
    template< class U > class _Allocator, template < int D > class _Space >
    struct AbstractHyperCube< _Dim, LinkType::Single, AccessScheme::Index,
    _Containment, _Allocator, _Space< _Dim > >: _Space < _Dim >
{
    enum {d = _Dim};
    ptrdiff_t upper, opponent, next;
    ptrdiff_t lower[_Dim * 2]; 
};
///terminate the recursion in the empty simplex (simplicial set) with
//dimension -1
template< template< class U, class V > class _Containment, 
    template< class U > class _Allocator, 
    template < int D > class _Space >
    struct AbstractHyperCube< -1, LinkType::Single,
    AccessScheme::Index, _Containment,
    _Allocator, _Space< -1 > >: _Space< -1 >
{
    enum { d = -1, };
};

template< int D, int M,
    template< int __D > class Space, template< class U, class V >
    class Containment, template< class U > class Allocator  >
    class	HyperCubeTree
{
    public:
        typedef Space< D > MetricTraitT;
        typedef Space< D > SpaceT;
        typedef typename MetricTraitT::ValueType ValueType;
        typedef typename MetricTraitT::KeyType KeyType;
        typedef typename MetricTraitT::Vector PointT;
        enum {
            d = D, alex_dimension = M,
            numofsubkeys = (sizeof(ValueType) * 8) * D / (D + M), 
            numoflevels = (sizeof(KeyType) * 8 / D), //(sizeof(ValueType) * 8),
            numofchilds = ipow< 2, d + alex_dimension >::eval,
            numofneighbours = ipow< 3, D >::eval - 1 };
        HyperCubeTree()
        {
            KeyType k = 0;
            PointT l_p; //l_p[0] = 0; l_p[1] = 0;
            reserve();
            HyperCube cube(k, l_p, 0);
            hypercubes[0].push_back(cube);
            counter[0] = 1;
        }
        struct HyperCube: SpaceT::Vector
        {
            public:
                typename IntegerChooser< 
                    (ipow< 3, D >::eval - 1) / 8 
                    >::IntegerType neighbourpattern;
                //using lossy state/phase space
                //communication for the intrawebz
                int childs[numofchilds];
                bool isleaf;
                KeyType key;
                //typename SpaceT::Vector key;
                int level;
                int weight;
                int age;
                int feature_index;
                PointT p; //quick hack reference/index
                //Containment::difference_type point_id;
                ptrdiff_t vertex_id; //reference to 0-simplex the cube is
                //containing; assuming one simplicial
                //decomposition/graph per tree layer

                HyperCube(const HyperCube& hc)
                {
                    for(int i = 0; i < numofchilds; i++)
                    {
                        childs[i] = hc.childs[i];
                    }
                    level = hc.level;
                    key = hc.key;
                    isleaf = hc.isleaf;
                    weight = hc.weight;
                    vertex_id = hc.vertex_id;
                    age = hc.age;

                }
                HyperCube()
                {
                    for(int i = 0; i < numofchilds; i++)
                    {
                        childs[i] = -1;
                    }
                    level = -1;
                    isleaf = false;
                    key = 0;
                    weight = 0;
                    age = 0;

                }

                HyperCube(KeyType &k, PointT _p, int _level)
                {
                    for(int i = 0; i < numofchilds; i++)
                    {
                        childs[i] = -1;
                    }
                    level = _level;
                    key = k;
                    p = _p;
                    weight = 0;

                }
        };

        inline void resize()
        {
            for(int i = 0; i < numoflevels; i++)
            {

                long long int l_size = 0x1; l_size <<= (i * 3);
                l_size = (l_size <= 1000000) ? l_size : 1000000;
                hypercubes[i].resize(l_size);
                counter[i] = 0;

            }
        }
        inline void reserve()
        {
            for(int i = 0; i < numoflevels; i++)
            {

                long long int l_size = 0x1; l_size <<= (i * 3);
                l_size = (l_size <= 100000) ? l_size : 100000;
                hypercubes[i].reserve(l_size);
                counter[i] = 0;

            }
        }

        inline void clear()
        {
            for(int i = 0; i < numoflevels; i++)
            {
                hypercubes[i].clear();
                counter[i] = 0;

            }
            KeyType k = 0;
            PointT l_p; //l_p[0] = 0; l_p[1] = 0;
            HyperCube cube(k, l_p, 0);
            hypercubes[0].push_back(cube);
            counter[0] = 1;
        }
        typedef typename Containment< HyperCube,
                Allocator< HyperCube > >::difference_type diff_type;
        typedef Containment< diff_type, Allocator< diff_type > > Indices;
        typedef Containment< KeyType, Allocator< KeyType > > Keys;
        typedef Containment< std::pair< KeyType, ptrdiff_t >, Allocator<
            std::pair< KeyType, ptrdiff_t > > > KeysnCubes;
        inline bool is(const KeyType key, int level)
        {
            HyperCube *hypercuberef = &hypercubes[0][0];
            int next;
            unsigned int subkey;

                for(int i = 1;i < level; i++)
                {
                    subkey = (numofchilds - 1) &
                        (key >> ((numoflevels - i) * D));
                    next = hypercuberef->childs[subkey];
                    if(next == -1)
                    {
                        return false;
                    }
                    hypercuberef = &(hypercubes[i][next]);
                }
                return true;
        }

        inline Keys getSamplers(const KeyType key, int level)
        {

            //generate vectors to add
            PointT l_adjvec[ ipow< 3, D >::eval];
            PointT l_basevec[ D ];
            for(int i = 0; i < D; i++)
            {
                l_basevec[i][i] = 0x1;
                l_basevec[i][i] <<= ( sizeof(l_basevec[i]) * 8 - level - 1);
            }
            for(int i = 0; i < ipow< 3, D >::eval; i++)
            {

            }
            Keys l_surrounding;
            //std::stack< int > down;
            //KeyType key;
            KeyType b[D]; //base
            //KeyType zero =  key & (~(ipow<2,
            //            (sizeof(KeyType) * 8 / D) - level>::eval - 1));
            for(int i = 0; i < D; i++) b[i] = 0b1 << 
                ((((sizeof(KeyType) * 8 / D) - level) * D)  - i);
            //to do: snaking
            //using point symmetry to mask out distance 2
            KeyType key_back;
            //key = MetricTraitT::morton_encode(p);
           // for(int i = 0; i < D; i++) l_mortonbase[i] = 
           // MetricTraitT::morton_encode(l_basevec[i]); 
            key_back = key;
            int i;

            HyperCube *hypercuberef = &hypercubes[0][0];
            HyperCube l_c;
            int pos, next, prev;
            unsigned int subkey;

            
            //for(int i = 0; i < ipow< 3, D >::eval; i++)
            //O'Rourke, "Computing Relative Neighborhood graph in the L_1
            //and L_\infty metrics," Pattern Recognition, 1982. 
            //
            KeyType neighbourhood[ipow< ipow< 2, 2 >::eval, D >::eval];
            //KeyType sampler[ipow< 3, D >::eval - 1];
            int po[D][2];
            PointT p = MetricTraitT::morton_decode(key);
            PointT n;
            Keys samplers;
            //samplers.push_back(key);
            //L_1
            for(int i = 0; i < (D); i++)
            {

                ValueType o = 0x1 << (((sizeof(KeyType) * 8 ) / D) 
                        - (level - 1) );
                //ValueType o = 0x1 << ( sizeof(l_basevec[i]) * 8 - level - 1);
                n = p;
                n[i] += o;
                samplers.push_back(MetricTraitT::morton_encode(n));
                n = p;
                n[i] -= o;
                samplers.push_back(MetricTraitT::morton_encode(n));
                po[i][0] = -o;
                po[i][1] = o;
            }
            
            
            for(int i = 0; i < ipow< 2, D >::eval; i++)
            {
                n = p;
                for(int j = 0; j < D; j++)
                    n[j] += po[j][i & (0x1 << j) ? 1 : 0];
                samplers.push_back(MetricTraitT::morton_encode(n));
            }

            for(int k = 0; k < D; k++)
            {
                for(int i = 0; i < ipow< 2, D >::eval; i++)
                {
                    n = p;
                    for(int j = 0; j < D; j++)
                        n[(j + 1) % D] += po[(j + 1) % D][i & (0x1 << j)
                            ? 1 : 0];
                    n[k] = p[k];
                    samplers.push_back(MetricTraitT::morton_encode(n));
                }
            }
           
            return samplers;
        }
        inline Keys getAdjacencies(const KeyType key, int level)
        {

            //generate vectors to add
            PointT l_adjvec[ ipow< 3, D >::eval];
            PointT l_basevec[ D ];
            for(int i = 0; i < D; i++)
            {
                l_basevec[i][i] = 0x1;
                l_basevec[i][i] <<= ( sizeof(l_basevec[i]) * 8 - level - 1);
            }
            for(int i = 0; i < ipow< 3, D >::eval; i++)
            {

            }
            Keys l_surrounding;
            //std::stack< int > down;
            //KeyType key;
            KeyType b[D]; //base
            //KeyType zero =  key & (~(ipow<2,
            //            (sizeof(KeyType) * 8 / D) - level>::eval - 1));
            for(int i = 0; i < D; i++) b[i] = 0b1 << 
                ((((sizeof(KeyType) * 8 / D) - level) * D)  - i);
            //to do: snaking
            //using point symmetry to mask out distance 2
            KeyType key_back;
            //key = MetricTraitT::morton_encode(p);
           // for(int i = 0; i < D; i++) l_mortonbase[i] = 
           // MetricTraitT::morton_encode(l_basevec[i]); 
            key_back = key;
            int i;

            HyperCube *hypercuberef = &hypercubes[0][0];
            HyperCube l_c;
            int pos, next, prev;
            unsigned int subkey;

            
            //for(int i = 0; i < ipow< 3, D >::eval; i++)
            //O'Rourke, "Computing Relative Neighborhood graph in the L_1
            //and L_\infty metrics," Pattern Recognition, 1982. 
            //
            KeyType neighbourhood[ipow< ipow< 2, 2 >::eval, D >::eval];
            //KeyType sampler[ipow< 3, D >::eval - 1];
            int po[D][2];
            PointT p = MetricTraitT::morton_decode(key);
            PointT n;
            Keys samplers;
            //samplers.push_back(key);
            //L_1
            for(int i = 0; i < (D); i++)
            {

                ValueType o = 0x1 << (((sizeof(KeyType) * 8 ) / D) 
                        - (level - 1) );
                //ValueType o = 0x1 << ( sizeof(l_basevec[i]) * 8 - level - 1);
                n = p;
                n[i] += o;
                samplers.push_back(MetricTraitT::morton_encode(n));
                n = p;
                n[i] -= o;
                samplers.push_back(MetricTraitT::morton_encode(n));
                po[i][0] = -o;
                po[i][1] = o;
            }
            
            
            for(int i = 0; i < ipow< 2, D >::eval; i++)
            {
                n = p;
                for(int j = 0; j < D; j++)
                    n[j] += po[j][i & (0x1 << j) ? 1 : 0];
                samplers.push_back(MetricTraitT::morton_encode(n));
            }

            for(int k = 0; k < D; k++)
            {
                for(int i = 0; i < ipow< 2, D >::eval; i++)
                {
                    n = p;
                    for(int j = 0; j < D; j++)
                        n[(j + 1) % D] += po[(j + 1) % D][i & (0x1 << j)
                            ? 1 : 0];
                    n[k] = p[k];
                    samplers.push_back(MetricTraitT::morton_encode(n));
                }
            }
           
            //return samplers;
            //for(int i = 0; i < ipow< ipow< 2, 2 >::eval, D >::eval; i++)
                //\f$ 2^n-tree level distance d 
                //{2^d}^n
                //mask out distance 2 within
                //p = \infty
            //for(int i = 0; i < ipow< 3, D >::eval - 1; i++)
            for(KeyType sampler: samplers)
                if(is(sampler, level)) l_surrounding.push_back(sampler);
            return l_surrounding;

        }
        /*
        inline KeysnCubes getAdjKC(const KeyType key, int level)
        {

            //generate vectors to add
            PointT l_adjvec[ ipow< 3, D >::eval];
            PointT l_basevec[ D ];
            for(int i = 0; i < D; i++)
            {
                l_basevec[i][i] = 0x1;
                l_basevec[i][i] <<= ( sizeof(l_basevec[i]) * 8 - level - 1);
            }
            for(int i = 0; i < ipow< 3, D >::eval; i++)
            {

            }
            Keys l_surrounding;
            //std::stack< int > down;
            //KeyType key;
            KeyType b[D]; //base
            //KeyType zero =  key & (~(ipow<2,
            //            (sizeof(KeyType) * 8t / D) - level>::eval - 1));
            for(int i = 0; i < D; i++) b[i] = 0b1 << 
                ((((sizeof(KeyType) * 8 / D) - level) * D)  - i);
            //to do: snaking
            //using point symmetry to mask out distance 2
            KeyType key_back;
            //key = MetricTraitT::morton_encode(p);
           // for(int i = 0; i < D; i++) l_mortonbase[i] = 
           // MetricTraitT::morton_encode(l_basevec[i]); 
            key_back = key;
            int i;

            HyperCube *hypercuberef = &hypercubes[0][0];
            HyperCube l_c;
            int pos, next, prev;
            unsigned int subkey;

            
            //for(int i = 0; i < ipow< 3, D >::eval; i++)
            //O'Rourke, "Computing Relative Neitghborhood graph in the L_1
            //and L_\infty metrics," Pattern Recognition, 1982. 
            //
            KeyType neighbourhood[ipow< ipow< 2, 2 >::eval, D >::eval];
            //KeyType sampler[ipow< 3, D >::eval - 1];
            int po[D][2];
            PointT p = MetricTraitT::morton_decode(key);
            PointT n;
            Keys samplers;
            samplers.push_back(key);
            //L_1
            for(int i = 0; i < (D); i++)
            {

                ValueType o = 0x1 << (((sizeof(KeyType) * 8 ) / D) 
                        - (level - 1) );
                //ValueType o = 0x1 << ( sizeof(l_basevec[i]) * 8 - level - 1);
                n = p;
                n[i] += o;
                samplers.push_back(MetricTraitT::morton_encode(n));
                n = p;
                n[i] -= o;
                samplers.push_back(MetricTraitT::morton_encode(n));
                po[i][0] = -o;
                po[i][1] = o;
            }
            for(int i = 0; i < ipow< 2, D >::eval; i++)
            {
                n = p;
                for(int j = 0; j < D; j++)
                    n[j] += po[j][i & (0x1 << j) ? 1 : 0];
                samplers.push_back(MetricTraitT::morton_encode(n));
            }
            //return samplers;
            //for(int i = 0; i < ipow< ipow< 2, 2 >::eval, D >::eval; i++)
                ;//\f$ 2^n-tree level distance d 
                //{2^d}^n
                //mask out distance 2 within
                //p = \infty
            //for(int i = 0; i < ipow< 3, D >::eval - 1; i++)
            for(KeyType sampler: samplers)
                if(is(sampler, level)) l_surrounding.push_back(sampler);

            return l_surrounding;

        }
*/

        inline bool isCube(const PointT &p, int level)
        {   
            if(level == 0) return true;

            KeyType key = MetricTraitT::morton_encode(p);
            HyperCube *hypercuberef = &hypercubes[0][0];

            int next;
            unsigned int subkey;

            for(int i = 1;i < level; i++)
            {
                subkey = (numofchilds - 1) & (key >> ((numoflevels - i)
                            * (D + M)));
                next = hypercuberef->childs[subkey];
                if(next == -1)
                {
                    return false;
                }
                hypercuberef = &(hypercubes[i][next]);
            }

            return false;
        }

        inline HyperCube& operator [](KeyType key)
        {
            return getCubebyIndex(key);
        }
        inline HyperCube& getCubebyIndex(KeyType key)
        {
            HyperCube *hypercubep = &hypercubes[0][0];
            for(int i = 1; i < numofsubkeys; i++)
            {
                int subkey = (numofchilds - 1) &
                    (key >> ((numoflevels - i) * (D + M)));
                int next = hypercubep->childs[subkey];
                if(next == -1) break;
                hypercubep = &hypercubes[i][next];

            }

            return *hypercubep;
        }

        inline bool insertPoint(const PointT &p, const KeyType key)
        {

            return true;

        }

        inline ptrdiff_t* getkNN(int k, const KeyType key, ptrdiff_t *NN)
        {
            //ptrdiff_t NN[k];
            HyperCube *hypercubep = &hypercubes[0][0];
            for(int i = 1; i < numofsubkeys; i++)
            {
                int subkey = (numofchilds - 1) &
                    (key >> ((numoflevels - i) * (D + M)));
                int next = hypercubep->childs[subkey];
                if(next == -1) break;
                hypercubep = &hypercubes[i][next];


            }
            return NN;
        }

        inline int insert(const KeyType k, const ptrdiff_t id)
        {
            PointT p_b;

            KeyType key_back;
            KeyType key = k;
            key_back = key;
            int i;

            HyperCube *hypercuberef = &hypercubes[0][0];
            HyperCube l_c;
            int pos, next;
            unsigned int subkey;

            //#pragma omp parallel for
            hypercuberef->weight++;
            hypercuberef->age = 0;
            for(i = 1;i <= numoflevels; i++)
            {
                subkey = (numofchilds - 1) &
                    (key >> ((numoflevels - i) * (D + M)));
                next = hypercuberef->childs[subkey];
                if(next == -1)
                {
                    {
                        pos = counter[i];
                        l_c.key = key_back;
                        l_c.level = i;
                        l_c.vertex_id = id;
                        hypercubes[i].push_back(l_c);
                        //pos = hypercubes[i].size() - 1;
                        hypercuberef->childs[subkey] = pos;

                        next = pos;
                        counter[i]++;
                    }
                }
                hypercuberef->age = 0;
                hypercuberef = &(hypercubes[i][next]);
            }

            hypercubes[numoflevels - 1][next].isleaf = true;
            hypercubes[numoflevels - 1][next].vertex_id = id;

            return pos;
        }

        Containment< HyperCube, Allocator< HyperCube > >
            hypercubes[numoflevels + 1];
        Containment< PointT, Allocator< PointT > > points;
        Containment< KeyType, Allocator< KeyType > > keys;
        int sub_keys[numofsubkeys];
        int counter[numoflevels + 1];

};



template < int D, template< int _D > class M >
struct MortonSpace: public TopologicalSpace< D >
{
    MortonSpace()
    {
    }

    MortonSpace(const MortonSpace &s)
    {
    }
};

template < template< int D > class M >
struct MortonSpace< 1, M >: public TopologicalSpace< 1 >
{
    typename M< 4 >::KeyType v;

    MortonSpace()
    {
        v = 0;
    }
    MortonSpace(const MortonSpace &s)
    {
        v = s.v;
    }
};

//specialise the simplices for morton space
/*
template< int D, template< int _D > class M > using MortonHalfSimplex =
AbstractHalfSimplex< D, LinkType::Single, AccessScheme::Index,
    std::vector, std::allocator, MortonSpace< D + 1, M > >;
template< int D, template< int _D > class M > using MortonHalfSimplicialComplex =
AbstractHalfSimplicialComplex< D, LinkType::Single, AccessScheme::Index,
    std::vector, std::allocator, MortonSpace< D + 1, M > >;
*/
template< int D, template< int _D > class M > using MortonHalfSimplex =
AbstractHalfSimplex< D, LinkType::Single, AccessScheme::Index,
    std::vector, std::allocator, MortonSpace< 1, M > >;
template< int D, template< int _D > class M > using MortonHalfSimplicialComplex =
AbstractHalfSimplicialComplex< D, LinkType::Single, AccessScheme::Index,
    std::vector, std::allocator, MortonSpace< 1, M > >;
template< int D, template< int _D > class M > using
MortonHalfSimplicialComplexTopologyTrait = AbstractHalfSimplicialComplexTopologyTrait<
MortonHalfSimplicialComplex< D, M > >;
template< int D, template< int _D > class M > using
MortonHalfSimplicialComplexTopologicalOperations = AbstractHalfSimplicialComplexTopologicalOperations<
MortonHalfSimplicialComplex< D, M > >;
template< int D, template< int _D > class M  > using
MortonHalfSimplicialComplexIterator = AbstractHalfSimplicialComplexIterator<
MortonHalfSimplicialComplex< D, M > >;

template < int D >
struct CompressedMetricalSpace: public TopologicalSpace< D >
{
    enum { d = D };
};





template < int D, template< int _D > class _M >
struct LinearSpaceCompressed: public MetricSpace< D, _M >
{
    enum{ d = D };

    typedef _M< D > Metric;
    typedef HyperCubeTree< D, 0, _M, std::vector, std::allocator > Tree;
    typedef typename Metric::KeyType KeyType;
    typedef typename Metric::ValueType ValueType;
    typedef typename Metric::PValueType PValueType;
    typedef typename Metric::template ElementT< d > PointT;
    typedef typename ExteriorPower< 1, d, LinearSpaceCompressed >::Vector
        Vector;
    typedef typename ExteriorPower< 0, 0, LinearSpaceCompressed>::Vector
        Scalar;
   
    typedef std::vector< std::pair< KeyType, KeyType > > Edges;
    typedef std::vector< std::array< KeyType, 5 > > Tets;

    typedef MortonHalfSimplicialComplex< D, _M > SC;
    template < int _D >
        using HalfSimplex = MortonHalfSimplex< _D, _M >;
  //  typedef MortonHalfSimplex< D, _M > S;
    typedef MortonHalfSimplicialComplexIterator< d, _M > Iterator;
    template < int _D >
        using Container = typename SC::template 
        Container< _D >;
    //using Vertices = MortonHalfSimplicialComplexIterator< 0, _M >::typename Container;



    /*
    struct Iterator: MortonHalfSimplicialComplexIterator < d, _M>
    {
    };
    */
    class Vertex: public MortonHalfSimplex< 0, _M >
    {
        public:
        enum { k = 1, d = D};
        //KeyType v;
        /*
        inline PValueType operator [](ptrdiff_t n) { 
            return Metric::morton_decode(this->v)[n];
        }
        */
        //todo: check performance differences (PValueType vs. Scalar)
        inline Scalar operator [] (ptrdiff_t n) {
            Scalar s;
            s.v[0] = Metric::morton_decode(this->v)[n];
            return s;
        }
        Vertex()
        {
        }
        Vertex( std::initializer_list< Scalar > val)
        {
        }
        Vertex(const PointT &p)
        {
            this->v = Metric::morton_encode(p);
        }
        Vertex(const Vertex& p)
        {
            this->v = p.v;
        }
        inline Vertex& operator = (const Vertex& p)
        {
            this->v = p.v;
            return *this;
        }
    };

    //typedef Vertex Vector; //i guess, every vector is a 1-cell, so ...
/*
    inline Iterator insert(Vertex &v)
    {
    }

    inline Iterator insert(Vector &v)
    {
    }
    */
    //inline Iterator insert(PointT &p)
    //
    
     template < int __D > inline Container< __D >&
        getSimplices(MortonHalfSimplex< __D, _M > &mhs)
        {
            Iterator iter(&simplicial_complex);
            return iter.getSimplices(mhs);
        }
    
    inline Iterator getIterator()
    {
        Iterator iter(&simplicial_complex);
        return iter;
    }
    
    inline Container< 0 >& getVertices()
    {
        MortonHalfSimplex< 0, _M > v;
        Iterator iter(&simplicial_complex);
        return iter.getSimplices(v);
    }
    

    inline Tets & getTets()
    {
        return m_tets;
    }
    inline Edges& getEdges()
    {
        return m_edges;
    }
    inline KeyType insert(KeyType k)
    {
        
        if(access_tree.is(k, 7)) return k;
        typename Tree::Keys neighbours =
            access_tree.getAdjacencies(k, 7);
        //Iter(simplicial_complex;
        Vertex v;
        v.v = k;
        MortonHalfSimplex< 0, _M > vv;
        vv.v = k;
        MortonHalfSimplex< 1, _M > he0, he1;
        std::vector< MortonHalfSimplex< 1, _M > > halfedges;
        MortonHalfSimplex< 0, _M > tetra[4];
        Iterator iter0, iter1;
        iter0 = simplicial_complex.insert(vv);
       // iter0 = simplicial_complex.insert(v); //which vertex to use?
        access_tree.insert(k, iter0[0]);
        //if(neighbours.size() < 3) return k;
        //
        int si = neighbours.size();


        std::multimap< int, KeyType > sorted_by_dist;
       // sorted_by_dist.insert(std::make_pair(0, k));
        for(KeyType lk: neighbours)
        {
            KeyType pk = (*(static_cast< typename SC::template
                        Container< 0 > * >
                        (simplicial_complex[0])))
                [access_tree[lk].vertex_id].v;
            sorted_by_dist.insert(
                    std::make_pair(Metric::dist(
                            Metric::morton_decode(k),
                            Metric::morton_decode(pk)
                            ), pk)
                    );
        }

        std::array< KeyType, 5 > tet; int c = 0;
        tet[0] = tet[1] = tet[2] = tet[3] = k;
        /*
        for(KeyType lk: neighbours)
        {
            he0.lower = access_tree[k].vertex_id;
            he1.lower = access_tree[lk].vertex_id;
            iter0 = simplicial_complex.insert(he0);
            iter1 = simplicial_complex.insert(he1);
            
            (*(static_cast< typename SC::template Container< 1 > * >(
                    simplicial_complex[1])))[iter0[1]].opponent = iter1[1];
            (*(static_cast< typename SC::template Container< 1 > * >(
                    simplicial_complex[1])))[iter1[1]].opponent = iter0[1];
                    
            
            KeyType k0 = (*(static_cast< typename SC::template
                        Container< 0 > * >
                        (simplicial_complex[0])))[he1.lower].v;
            KeyType k1 = (*(static_cast< typename SC::template
                        Container< 0 > * >
                        (simplicial_complex[0])))[he0.lower].v;
            tet[c % 4] = k0; tet[(c+1) % 4] = k1;
            c++;

            m_edges.push_back(std::make_pair(k1, k0));

        }
        */
        /*
        if(neighbours.size() >= 1){
            KeyType k0, k1;
            typename std::multimap< int, KeyType >::iterator beg =
                sorted_by_dist.begin();
            k0 = (*(static_cast< typename SC::template
                        Container< 0 > * >
                        (simplicial_complex[0])))
                [access_tree[beg->second].vertex_id].v;
            beg++;
            
            k1 = (*(static_cast< typename SC::template
                        Container< 0 > * >
                        (simplicial_complex[0])))
                [access_tree[beg->second].vertex_id].v;
           // beg++;

            m_edges.push_back(std::make_pair(k1, k0));
        }
*/

        if(neighbours.size() > d){
            /*
            KeyType tet0[4];
            tet[0] = (*(static_cast< typename SC::template
                        Container< 0 > * >
                        (simplicial_complex[0])))
                [access_tree[neighbours[0]].vertex_id].v;
            tet[1] = (*(static_cast< typename SC::template
                        Container< 0 > * >
                        (simplicial_complex[0])))
                [access_tree[neighbours[1]].vertex_id].v;
            tet[2] = (*(static_cast< typename SC::template
                        Container< 0 > * >
                        (simplicial_complex[0])))
                [access_tree[neighbours[2]].vertex_id].v;
            tet[3] = (*(static_cast< typename SC::template
                        Container< 0 > * >
                        (simplicial_complex[0])))
                [access_tree[neighbours[3]].vertex_id].v;
            m_tets.push_back(tet);
            */

            //need to merge simplicial complexes
            //the unification process is driven by metrics;
            //we're using powerset graph topology
            /*
            void* l_simplex[ d + 1 ];
            l_simplex[0] = NULL; //dimension -1 -> empty set
            for(int i = 1; i < d; i++)
            {
                l_simplex[i] = new HalfSimplex<
            }
            */
            HalfSimplex< d > simp;
            typename std::multimap< int, KeyType >::iterator beg =
                sorted_by_dist.begin();
            tetra[0] = (*(static_cast< typename SC::template
                        Container< 0 > * >
                        (simplicial_complex[0])))
                [access_tree[beg->second].vertex_id];
            tet[0] = (*(static_cast< typename SC::template
                        Container< 0 > * >
                        (simplicial_complex[0])))
                [access_tree[beg->second].vertex_id].v;
            beg++;

            tetra[1] = (*(static_cast< typename SC::template
                        Container< 0 > * >
                        (simplicial_complex[0])))
                [access_tree[beg->second].vertex_id];
            tet[1] = (*(static_cast< typename SC::template
                        Container< 0 > * >
                        (simplicial_complex[0])))
                [access_tree[beg->second].vertex_id].v;
            beg++;

            tetra[2] = (*(static_cast< typename SC::template
                        Container< 0 > * >
                        (simplicial_complex[0])))
                [access_tree[beg->second].vertex_id];
            tet[2] = (*(static_cast< typename SC::template
                        Container< 0 > * >
                        (simplicial_complex[0])))
                [access_tree[beg->second].vertex_id].v;
            beg++;

            tetra[3] = (*(static_cast< typename SC::template
                        Container< 0 > * >
                        (simplicial_complex[0])))
                [access_tree[beg->second].vertex_id];
            tet[3] = (*(static_cast< typename SC::template
                        Container< 0 > * >
                        (simplicial_complex[0])))
                [access_tree[beg->second].vertex_id].v;
            iter0.make(tetra);
            m_tets.push_back(tet);


            for(int i = 0; i < 4; i++)
                m_edges.push_back(std::make_pair(tet[i % d], tet[(i+1)%d]));
            m_edges.push_back(std::make_pair(tet[0], tet[2]));
            m_edges.push_back(std::make_pair(tet[1], tet[3]));


        }


        return k;



        //access_tree.insert(k, -1);

        //check is()
        //if is() return;
        //insert()
        //getNeighbourhood
        //simplicial decompose depth/height-dependend
        //numbers of neighbours -> dimension for triangulation
        //retriangulate
        
        //return access_tree.insert(k, -1);
    }
    inline KeyType insert(PointT &p)
    {
        KeyType k;
        Vertex v;
        k = Metric::morton_encode(p);
        v.v = k;

        ptrdiff_t NN[D * 3];
        access_tree.getkNN(D * 3, k, NN);
        
        std::map< KeyType, ptrdiff_t > m;

        //for(int i = 0; i++; i < D * 3)
        {

           // std::cout << (simplicial_complex.getContainer<0>(0))[i].v;//[i].v;
            //typename HalfSimplicialComplex< D >::Container< 0 > c;
          //  std::cout << (*((typename HalfSimplicialComplex< D >::template Container< 0 > *)
          //      simplicial_complex.simplex_containers[0]))[i].v;
          //  std::cout << simplicial_complex[i];
            //m.insert(make_pair(((HalfSimplicialComplex::Container<0> *)(simplicial_complex[0]))[i].v, NN[i]));
          //  std::cout << simplicial_complex.bla;
        }

        //getkNN() -> get 1-Ring for each NN
        //-> intersection of all 1-Ring
        //-> proj. space
        //access_tree.insert(

        return v.v;
    }
    inline void clear()
    {
        access_tree.clear();
    }
    Tree access_tree; //access tree for 0-simplices
    Edges m_edges;
    Tets m_tets;
    MortonHalfSimplicialComplex< D, _M > simplicial_complex; //replace
    //std-containers with proper trees with metrical traits
    Vector e[D]; //basis

};
template< int D >
using EuklidianSpaceCompressedFuint16 = LinearSpaceCompressed< D,
      SimpleEuklidianMetricFuint16 >;
template< int D >
using EuklidianSpaceCompressedFuint32 = LinearSpaceCompressed< D,
      SimpleEuklidianMetricFuint32 >;
template< int D >
using EuklidianSpaceCompressedFuint64 = LinearSpaceCompressed< D,
      SimpleEuklidianMetricFuint64 >;


// vector <-> polynome / functional / function
//todo AS <-> multivectors/pseudoscalar (wedge product &
//other clifford algebra stuff) half simplices -> SO(n)
//wedge product, wedge sum ... bouquet of circles
//V = \bigoplus_{n \in \mathbb{N}} V_n <- make graded vector spaces
//projective spaces -> quotient spaces of topological spaces(?)
//-> finite fields -> elliptic curves
//Grassmannian(k, V) is a algebraic subvariety of projective space P(Λ^kV) ->
//Pluecker embeddingalComplexIterator >::doit(*this);	
//exterior algebra -> simplicial complex
//derived from power set(?) correspondency \f$P(X) \cong \{0, 1\}^X\f$
//from char. function to isomorphism.
//metrical dependency of the cell decomposition
#endif
