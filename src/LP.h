#ifndef LP_H
#define LP_H

#include <vector>
#include <tuple>
#include <utility>
#include <list>

class LP2D
{
  public:

    using Bound         = std::pair<double,double>;         // min, max
    using Constradouble = std::tuple<double,double,double>; // a, b, c

    double solve( const std::vector<Constraint> &constraints );

  private:

    void    init   ( const std::vector<Constraint> &constraints );
    double  iterate();

    list<Constraint> mIp; // I^+
    list<Constraint> mIn; // I^-

    double mXl; // left boundary
    double mXr; // right boundary

    double  mXm;
    double  mAlpha;
    double  mBeta;
    Bound   mS;
    Bound   mT;
};

#endif
