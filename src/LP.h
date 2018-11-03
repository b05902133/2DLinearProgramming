#ifndef LP_H
#define LP_H

#include <vector>
#include <tuple>
#include <utility>
#include <list>

class LP2D
{
  public:

    enum class IType
    {
      plus,
      minus
    };

    /*
     *  The position of x* relative to x_m
     */
    enum class RelativePosition
    {
      noSolution,
      left,
      right,
      equal
    };

    using Constraint = std::tuple<double,double,double>; // a, b, c

    double solve( const std::vector<Constraint> &constraints );

  private:

    using Bound           = std::pair<double,double>; // min, max
    using ConstraintList  = std::list<Constraint>;
    using Iterator        = ConstraintList::iterator;
    using RxSource        = std::tuple<IType,Iterator,Iterator>;

    static constexpr size_t mSizeSmall = 2;

    void    init   ( const std::vector<Constraint> &constraints );
    double  iterate();

    void              collectRxs( ConstraintList &I, IType iType,
                                  std::vector<double> &rxs,
                                  std::vector<RxSource> &rxSources );
    void              evalParams();
    RelativePosition  evalOptPosition();
    double            solveReduced();

    void removeConstraint( ConstraintList &constraints, Iterator &it1, Iterator &it2, bool removeIt1 );

    ConstraintList mIp; // I^+
    ConstraintList mIn; // I^-

    double mXl; // left boundary
    double mXr; // right boundary

    double  mXm;
    double  mAlpha;
    double  mBeta;
    Bound   mS;
    Bound   mT;
};

#endif
