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

#ifndef LE_SPACE_H
#define LE_SPACE_H

#include <cstddef>
#include <cstdint>
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
struct ibinom{ 
    enum { 
        eval = ifact< N >::eval / ( ifact< K >::eval *  ifact< N - K >::eval)
    };
};
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

//quick hack integer chooser stuff
//
template< int S > struct SimpleLargerInteger{
	unsigned long long values[S/64];
};

template< int S > struct IntegerChooser
{
	typedef typename IF< S <= sizeof(unsigned char), unsigned char,
		IF< S <= sizeof(unsigned short), unsigned short,
		IF< S <= sizeof(unsigned int), unsigned int,
		IF< S <= sizeof(unsigned long long), unsigned long long,
        SimpleLargerInteger< S > > > > >::RET IntegerType;
};

template< > struct IntegerChooser< 3 >
{
    enum { size = 3 };
    typedef unsigned int IntegerType;
};

template< > struct IntegerChooser< 4 >
{
	enum { size = 4};
	typedef unsigned int IntegerType;
};


template< > struct IntegerChooser< 6 >
{
	enum { size = 6};
	typedef __attribute__ ((aligned (1))) unsigned long long IntegerType;
};

template< > struct IntegerChooser< 8 >
{
	enum { size = 8};
	typedef __attribute__ ((aligned (32))) unsigned long long IntegerType;
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
struct Set {};
template< int _Dim >
struct TopologicalSpace: Set< _Dim > {};
template< int _Dim >
struct UniformSpace: TopologicalSpace< _Dim > {};
template< int _Dim, template < int __D > class _M  >
struct MetricSpace: UniformSpace < _Dim > {};

/**
 * kinda abstract simplex(?), one example for an element
 */
template< int _Dim, typename _Trait, LinkType _LType, AccessScheme _AScheme,
    template< class U , class V > class Containment , template< class U >
    class Allocator, class _Space > class AbstractSimplicialComplex{};
template< typename _ASD > struct AbstractSimplicialComplexTopologyTrait{};
template< typename _Iterator > 
struct AbstractSimplicialComplexIteratorFunctionTrait{};
template< typename _ASD > struct AbstractSimplicialComplexIterator{};
template< int _Dim, LinkType _LType, AccessScheme _AScheme,  
    template< class U, class V > class _Containment,
    template< class U > class _Allocator, class _Space >
    struct AbstractSimplex: _Space {};

///why d + 1 you may ask. 'cause \f$ {d+1 \choose d} = d+1\f$ , u kno'. 
//thus we're terminating the recursion at dimension = -1
template< int _Dim,  template< class U, class V > class _Containment,
    template< class U > class _Allocator, template < int D > class _Space >
    struct AbstractSimplex< _Dim, LinkType::Single, AccessScheme::Index,
    _Containment, _Allocator, _Space< _Dim > >: _Space < _Dim >
{
    enum {d = _Dim};
    ptrdiff_t upper, opponent, next;
    ptrdiff_t lower[_Dim + 1]; 
};

///terminate the recursion in the empty simplex (simplicial set) with
//dimension -1
template< template< class U, class V > class _Containment, 
    template< class U > class _Allocator, 
    template < int D > class _Space >
    struct AbstractSimplex< -1, LinkType::Single,
    AccessScheme::Index, _Containment,
    _Allocator, _Space< -1 > >: _Space< -1 >
{
    enum { d = -1, };
};


///why d + 1 you may ask. 'cause \f$ {d+1 \choose d} = d+1\f$ , u kno'. thus
//we're terminating the recursion at dimension = -1
template< int _Dim,  template< class U, class V > class _Containment,
    template< class U > class _Allocator >
    struct AbstractSimplex< _Dim, LinkType::Single,
    AccessScheme::Index, _Containment, _Allocator, Set< _Dim > >: Set < _Dim >
{
    enum {d = _Dim};
    ptrdiff_t upper, opponent, next;
    ptrdiff_t lower[_Dim + 1];
};

///terminate the recursion in the empty simplex (simplicial set) with dimension -1
template< template< class U, class V > class _Containment,
    template< class U > class _Allocator >
struct AbstractSimplex< -1, LinkType::Single, AccessScheme::Index, _Containment,
    _Allocator, Set< -1 > >: Set< -1 >
{
    enum { d = -1, };
};
template< int _D, class _SC > struct ContainerFiller
{ 
    static inline void fill(_SC& s)
    {
        typedef AbstractSimplex< _D - 1, LinkType::Single, 
                AccessScheme::Index, _SC::template _Containment,
                _SC::template _Allocator,
                typename _SC::template Space< _D - 1 > > AbstractSimplexT;
        typedef typename _SC::_Trait::template 
            SimplexT< _D, AbstractSimplexT > Simplex;
        typedef typename _SC::template 
            _Containment< Simplex, 
            _SC::template _Allocator< Simplex > > Simplices;
        s.simplex_containers[_D - 1] = new Simplices;
        ContainerFiller< _D - 1, _SC >::fill(s);
    };	
};

template< class _SC > 
struct ContainerFiller< 0, _SC >{ static inline void fill(_SC&){}};

template< int _D, class _IT, class  _SC > struct IterFiller
{
    static inline void fill(_IT& i)
    {
        typedef  AbstractSimplex< _D - 1, LinkType::Single,
                 AccessScheme::Index, _SC::template _Containment,
                 _SC::template _Allocator, 
                 typename _SC::template Space < _D - 1 > > AbstractSimplexT;
        i.iterdata[_D] = new AbstractSimplexT;
        IterFiller< _D - 1, _IT, _SC >::fill(i);
    }
};

template< class _IT, class _SC > struct IterFiller< 0, _IT, _SC >{
    static inline void fill(_IT&){}
};

template< int _Dim, class _Trait,
    template< class U , class V > class _Containment,
    template< class U > class _Allocator, class _Space> 
    class AbstractSimplicialComplex< _Dim, _Trait, 
    LinkType::Single, AccessScheme::Index,
    _Containment, _Allocator, _Space >
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
template< int _Dim, class _Trait, LinkType _LType, AccessScheme _AScheme,
    template< class U , class V > class _Containment,
    template< class U > class _Allocator, class _Space> 
class AbstractSimplicialComplexIterator<
AbstractSimplicialComplex< _Dim, _Trait, _LType,
    _AScheme, _Containment, _Allocator, _Space > >
{
    public:
        typedef AbstractSimplicialComplex< _Dim, _Trait, _LType,
                _AScheme, _Containment, _Allocator, _Space > ASCT;
        typedef AbstractSimplicialComplexTopologyTrait< ASCT > ASCTopoT;
        typedef AbstractSimplicialComplexIteratorFunctionTrait<
            AbstractSimplicialComplexIterator > IterFuncT;

        AbstractSimplicialComplexIterator()
        {
            IterFiller< _Dim, AbstractSimplicialComplexIterator,
                ASCT >::fill(this);
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
                ASCTopoT::template simplexCCV<_D,
                    AbstractSimplicialComplexIterator >::doit(this);	
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

template< int _Dim, class _Trait,
    template< class U , class V > class _Containment,
    template< class U > class _Allocator, class _Space> 
class AbstractSimplicialComplexTopologyTrait<  
AbstractSimplicialComplex< _Dim, _Trait, LinkType::Single,
    AccessScheme::Index, _Containment, _Allocator, _Space > >
{
    public:

        typedef AbstractSimplicialComplex< _Dim, _Trait,
                LinkType::Single, AccessScheme::Index, _Containment,
                _Allocator, _Space > ASCT;
        typedef AbstractSimplicialComplexIterator< ASCT > IterT;
        //use const_expr for operations
        template< int _D, class _It >
            struct simplexFlip
            {
                inline bool doit(IterT &iter)
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
            struct simplexAlign
            {
                inline bool doit(IterT &iter)
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


#include <vector>

template< int D > using Simplex = AbstractSimplex< D, LinkType::Single,
    AccessScheme::Index, std::vector, std::allocator, TopologicalSpace< D > >;


//simple implementation of a point
//kinda n-vector stuff if used in a vector-space

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
			for(auto v : val) values[i++] = v;
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
    //typedef double ValueType;
    //should i pull the Tensor Definitions from the exterior algebra?
    //ok, let's do it :|
    typedef typename ExteriorPower< 0, 0, LinearSpace>::Vector Scalar;
    /*
    struct Scalar: Simplex< 0 >
    {
        enum { k = 0, d = D };

    };*/
    struct Vertex: Simplex< 1 >//derive from simplex, cause every
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

    typedef Vertex Vector; //i guess, every vector is a 1-cell, so ...

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
    template < int N > using ElementT = _Point < N >;
    typedef typename IntegerChooser< D * sizeof(ValueType) >::IntegerType 
        KeyType;
    typedef Point Vector;
    static inline Vector dist(Point p0, Point p1)
    {
        Vector d;
        for(ptrdiff_t i = 0; i < D; i++)d[i] = p1[i] - p0[i];
        return d;
    }

    //draftly putting Morton Code into the metric
    //only makes sense with integer

    static inline ValueType interleavetokey(Point p)
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
template< int D > using SimplePointFint64 = SimplePoint< D, int64_t >;
template< int D > using SimpleEuklidianMetricFdouble = EuklidianMetric< D,
    double,
    SimplePointFdouble >;
template< int D > using SimpleEuklidianMetricFint64 = EuklidianMetric< D,
    int64_t, 
    SimplePointFint64 >;

//template< int D >
//using EuklidianSpace = LinearSpace<D, EuklidianMetric>;
template< int D >
using SimpleEuklidianSpaceFdouble = LinearSpace< D,
      SimpleEuklidianMetricFdouble >;


template< int K, int D, class S > struct ExteriorPower: QuotientSpace< K, D >
{
    enum{ d = nchoosek< D, K >::eval };
    typedef typename S::Metric Metric;
    struct Vector;

    typedef typename ExteriorPower< 0, 0, S >::Vector Scalar;
    //typedef typename ExteriorPower< 0, 0, S >::Vector Scalar;

    struct Vector: Simplex< K >
    {
        typename Metric::template ElementT< d > v;
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

template< int D, int M, typename Type, template < int _D, typename _Type,
    typename ... _CArgs > class Point,
    template< int __D, typename __Type,
    template < int ___D, typename ___Type, typename ... _CArgs >
    class __Point > class MetricTrait, template< class U, class V >
    class Containment, template< class U > class Allocator  >
class	HyperCubeTree
{
	public:
		typedef typename MetricTrait< D, Type, Point >::ValueType ValueType;
		typedef typename MetricTrait< D, Type, Point >::KeyType KeyType;
		typedef typename MetricTrait< D, Type, Point>::Point PointT;
		typedef MetricTrait< D, Type, Point > MetricTraitT;
        enum {
            dimension = D, alex_dimension = M,
            numofsubkeys = (sizeof(ValueType) * 8) * D / (D + M), 
            numoflevels = (sizeof(ValueType) * 8),
            numofchilds = ipow< 2, dimension + alex_dimension >::value,
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
        /**
         *
         * quick hack simple graph
         *
         */

		struct HyperCube
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
					p.set(hc.p);
                    isleaf = hc.isleaf;
                    weight = hc.weight;
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
        inline Indices getAdjacencies(const PointT &p, int level)
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
            Indices l_surrounding;
            KeyType key;
            KeyType key_back;
            key = MetricTraitT::interleavetokey(p);
			key_back = key;
            int i;

			HyperCube *hypercuberef = &hypercubes[0][0];
            HyperCube l_c;
			int pos, next, prev;
            unsigned int subkey;

            for(i = 1;i < level; i++)
            {
                subkey = (numofchilds - 1) & (key >> ((numoflevels - i)
                            * (D + M)));
                prev = next;
                next = hypercuberef->childs[subkey];
				if(next == -1)
				{
                    return l_surrounding;
				}
				hypercuberef = &(hypercubes[i][next]);
			}

            HyperCube *parent = &(hypercubes[i][prev]);
            l_surrounding.push_back(next);

        }

        inline bool isCube(const PointT &p, int level)
        {   
            if(level == 0) return true;

            KeyType key = MetricTraitT::interleavetokey(p);
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

        inline int insertPoint(const KeyType k)
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
                        l_c.p.set(p_b);
                        l_c.level = i;
                        hypercubes[i].push_back(l_c);
                        hypercuberef->childs[subkey] = pos;

                        next = pos;
                        counter[i]++;
                    }
				}
                hypercuberef->age = 0;
				hypercuberef = &(hypercubes[i][next]);
			}

            hypercubes[numoflevels - 1][next].isleaf = true;
			return pos;
		}

	    inline int insertPoint(const PointT &p)
		{
            PointT p_b;
            p_b.set(p);
            
            KeyType key_back;
            KeyType key = MetricTraitT::interleavetokey(p);
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
                        l_c.p.set(p_b);
                        l_c.level = i;
                        hypercubes[i].push_back(l_c);
                        hypercuberef->childs[subkey] = pos;

                        next = pos;
                        counter[i]++;
                    }
				}
                hypercuberef->age = 0;
				hypercuberef = &(hypercubes[i][next]);
			}

            hypercubes[numoflevels - 1][next].isleaf = true;
			return pos;
		}


        inline void genKeys(const typename Containment< 
                PointT, Allocator< PointT > >::difference_type start,
                const typename Containment< PointT,
                Allocator< PointT > >::difference_type stop)
        {
#pragma omp parallel
            for(; start < stop; start++)
                keys[start] = MetricTraitT::interleavetokey(points[start]);

        }

        inline void genTree(const typename Containment<
                PointT, Allocator< PointT > >::difference_type start,
                const typename Containment< PointT,
                Allocator< PointT > >::difference_type stop)
        {
            genKeys(start, stop);
#pragma omp parallel
            for(; start < stop; start++)
                insertPoint(points[start], keys[start]);

        }



        typedef Containment< HyperCube, Allocator< HyperCube > > Layer;

        inline Layer& getLayer(int level)
        {
            return hypercubes[level % numoflevels];
        }

        Containment< HyperCube, Allocator< HyperCube > >
            hypercubes[numoflevels + 1];
        Containment< PointT, Allocator< PointT > > points;
        Containment< KeyType, Allocator< KeyType > > keys;
        int sub_keys[numofsubkeys];
        int counter[numoflevels + 1];
	
};
// vector <-> polynome / functional / function
//todo AS <-> multivectors/pseudoscalar (wedge product &
//other clifford algebra stuff) half simplices -> SO(n)
//wedge product, wedge sum ... bouquet of circles
//V = \bigoplus_{n \in \mathbb{N}} V_n <- make graded vector spaces
//projective spaces -> quotient spaces of topological spaces(?)
//-> finite fields -> elliptic curves
//Grassmannian(k, V) is a algebraic subvariety of projective space P(Λ^kV) ->
//Pluecker embedding
//exterior algebra -> simplicial complex
//derived from power set(?) correspondency \f$P(X) \cong \{0, 1\}^X\f$
//from char. function to isomorphism.
//metrical dependency of the cell decomposition
#endif
