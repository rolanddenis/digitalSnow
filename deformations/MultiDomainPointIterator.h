#include <iterator> // Iterator traits
#include <boost/iterator/iterator_facade.hpp>

template < typename TDomainIterator >
class MultiDomainPointIterator
  : public boost::iterator_facade <
      MultiDomainPointIterator<TDomainIterator>,
      typename std::iterator_traits<TDomainIterator>::value_type::Point const,
      std::input_iterator_tag
    >
{
public:
  using Self   = MultiDomainPointIterator<TDomainIterator>;
  using DomainIterator = typename std::remove_const< typename std::decay<TDomainIterator>::type >::type;
  using Domain = typename std::iterator_traits<DomainIterator>::value_type;
  using Point  = typename Domain::Point;
  using PointIterator  = typename Domain::ConstIterator;
  using Value  = Point;
  using Reference = Point const&;

  //TODO: check that Domain is a CDomain model.

public:

  /** Default constructor.
   *
   * Creates an uninitialized (and thus invalid) iterator.
   */
  MultiDomainPointIterator() = delete;

  /** Begin iterator.
   *
   * Creates a multi-domain iterator pointing to the first point.
   *
   * @param aDomainBeginIt  An iterator pointing to the first domain.
   * @param aDomainEndIt    An iterator pointing after the last domain.
   */
  MultiDomainPointIterator( DomainIterator const& aDomainBeginIt, DomainIterator const& aDomainEndIt )
      : myDomainIt( aDomainBeginIt )
      , myDomainEndIt( aDomainEndIt )
      , myPointIt( myDomainIt->begin() )
      , myPointEndIt( myDomainIt->end() )
    {
    }

  /** Past-the-end iterator.
   *
   * Creates a multi-domain iterator pointing after the last point.
   *
   * @param aDomainBeginIt  An iterator pointing to the first domain.
   * @param aDomainEndIt    An iterator pointing after the last domain.
   */
  MultiDomainPointIterator( DomainIterator const& aDomainBeginIt, DomainIterator const& aDomainEndIt, bool /* last */ )
      : myDomainIt( aDomainEndIt )
      , myDomainEndIt( aDomainEndIt )
      , myPointIt( aDomainBeginIt->begin() )
      , myPointEndIt( aDomainEndIt->end() )
    {
    }

  /// Default copy constructor.
  MultiDomainPointIterator( Self const& ) = default;

  /// Default move constructor.
  MultiDomainPointIterator( Self &&     ) = default;

  /// Default copy assignment operator.
  Self& operator= ( Self const& ) = default;

  /// Default move assignment operator.
  Self& operator= ( Self &&     ) = default;

  /// Default destructor.
  ~MultiDomainPointIterator() = default;

private:
  friend class boost::iterator_core_access; // Friend access to boost iterator_facade.

  /// Iterator increment (boost::iterator_facade interface).
  void increment()
    {
      if ( ++myPointIt == myPointEndIt )
        /* We have to check if we are not pointing to the past-the-end domain
         * in order to avoid invalid domain pointer dereference.
         */
        if ( ++myDomainIt != myDomainEndIt )
          {
            myPointIt     = myDomainIt->begin();
            myPointEndIt  = myDomainIt->end();
          }
    }

  /// Iterator comparison (boost::iterator_facade interface).
  inline
  bool equal( Self const& other ) const
    {
      /* Two iterators are equal either:
       * - if they point to the same current domain and point,
       * - if they both point to the past-the-end domain.
       */
      return myDomainIt == other.myDomainIt &&
        ( myDomainIt == myDomainEndIt || myPointIt == other.myPointIt );
    }

  /// Iterator dereference (boost::iterator_facade interface).
  inline
  Reference dereference() const
    {
      return *myPointIt;
    }

private:
  DomainIterator  myDomainIt;     ///< Current domain iterator.
  DomainIterator  myDomainEndIt;  ///< Past-the-end domain iterator.
  PointIterator   myPointIt;      ///< Current point iterator.
  PointIterator   myPointEndIt;   ///< Past-the-end point iterator.  

};

/* GNU coding style */
/* vim: set ts=2 sw=2 expandtab cindent cinoptions=>4,n-2,{2,^-2,:2,=2,g0,h2,p5,t0,+2,(0,u0,w1,m1 : */
