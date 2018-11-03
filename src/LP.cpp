#include "LP.h"

#include <limits>
#include <cassert>
#include <algorithm>
#include <vector>
using namespace std;

#include "Select.h"

// public member functions
double LP2D::solve( const std::vector<Constraint> &constraints )
{
  init( constraints );

  return iterate();
}
// end public member functions

// private member functions
void LP2D::init( const std::vector<Constraint> &constraints )
{
  mXl = numeric_limits<double>::lowest();
  mXr = numeric_limits<double>::max();

  for( const Constraint &constraint : constraints )
  {
     double b = get<1>( constraint );

     if( b == 0 )
     {
       double a = get<0>( constraint );
       double c = get<2>( constraint );
       double x = abs( c / a );

       if( x < mXl || x > mXr ) continue;

       if( a > 0 ) // x <= abc( c / a )
         mXr = min( mXr, x );
       else if( a < 0 ) // x >= abc( c / a )
         mXl = max( mXl, x );
     }
     else if( b > 0 )
       mIp.push_back( constraint );
     else // b < 0
       mIn.push_back( constraint );
  }
}

double LP2D::iterate()
{
  if( mXr < mXl ) return numeric_limits<double>::max();

  while( !mIn.empty() && mIp.size() + mIn.size() > mSizeSmall )
  {
    vector<double>    rxs;
    vector<RxSource>  rxSources;
    Select            selectEngine;
    RelativePosition  optPosition;

    collectRxs( mIp, IType::plus,   rxs, rxSources );
    collectRxs( mIn, IType::minus,  rxs, rxSources );

    if( rxs.empty() ) continue;

    mXm = selectEngine.select( rxs, rxs.size() / 2 + 1 );

    evalParams();

    optPosition = evalOptPosition();

    if( optPosition == RelativePosition::noSolution ) return numeric_limits<double>::max();
    if( optPosition == RelativePosition::equal      )
    {
      size_t i;

      for( i = 0 ; i < rxs.size() ; ++i )
         if( rxs[i] == mXm ) break;

      const Constraint  &constraint = *get<1>( rxSources[i] );
      const double      a           = get<0>( constraint    );
      const double      b           = get<1>( constraint    );
      const double      c           = get<2>( constraint    );

      return ( c - a * mXm ) / b;
    }

    if( optPosition == RelativePosition::left )
    {
      for( size_t i = 0 ; i < rxs.size() ; ++i )
      {
         if( rxs[i] < mXm ) continue;

         Iterator it1 = get<1>( rxSources[i] );
         Iterator it2 = get<2>( rxSources[i] );
         double   a1  = get<0>( *it1 );
         double   b1  = get<1>( *it1 );
         double   a2  = get<0>( *it2 );
         double   b2  = get<1>( *it2 );

         if( get<0>( rxSources[i] ) == IType::plus )
           removeConstraint( mIp, it1, it2, -a1 / b1 < -a2 / b2 );
         else
           removeConstraint( mIn, it1, it2, -a1 / b1 > -a2 / b2 );
      }
      mXr = mXm;
    }
    else // RelativePosition::right
    {
      for( size_t i = 0 ; i < rxs.size() ; ++i )
      {
         if( rxs[i] > mXm ) continue;

         Iterator it1 = get<1>( rxSources[i] );
         Iterator it2 = get<2>( rxSources[i] );
         double   a1  = get<0>( *it1 );
         double   b1  = get<1>( *it1 );
         double   a2  = get<0>( *it2 );
         double   b2  = get<1>( *it2 );

         if( get<0>( rxSources[i] ) == IType::plus )
           removeConstraint( mIp, it1, it2, -a1 / b1 > -a2 / b2 );
         else
           removeConstraint( mIn, it1, it2, -a1 / b1 < -a2 / b2 );
      }
      mXl = mXm;
    }
  }
  return solveReduced();
}

/*
 *  a1 * x + b1 * y = c1
 *  a2 * x + b2 * y = c2
 *
 *  => x = ( b2 * c1 - b1 * c2 ) / ( a1 * b2 - a2 * b1 )
 */
void LP2D::collectRxs(  ConstraintList &I, IType iType,
                        std::vector<double> &rxs,
                        std::vector<RxSource> &rxSources )
{
  Iterator it1;
  Iterator it2;

  for(  it1 = I.begin(), it2 = ++I.begin() ;
        it1 != I.end() && it2 != I.end() ;
        ++( ++it1 ), ++( ++it2 ) )
  {
     double a1 = get<0>( *it1 );
     double b1 = get<1>( *it1 );
     double c1 = get<2>( *it1 );
     double a2 = get<0>( *it2 );
     double b2 = get<1>( *it2 );
     double c2 = get<2>( *it2 );
     double rx;

     if( a1 == a2 && b1 == b2 ) // parallel
     {
       removeConstraint(  I, it1, it2,
                          ( ( iType == IType::plus   && c1 > c2 ) ||
                            ( iType == IType::minus  && c1 < c2 ) ) );
       continue;
     }

     rx = ( b2 * c1 - b1 * c2 ) / ( a1 * b2 - a2 * b1 );

     if( rx < mXl )
     {
       removeConstraint(  I, it1, it2,
                          ( ( iType == IType::plus   && ( -a1 / b1 > -a2 / b2 ) ) ||
                            ( iType == IType::minus  && ( -a1 / b1 < -a2 / b2 ) ) ) );
       continue;
     }

     if( rx > mXr )
     {
       removeConstraint(  I, it1, it2,
                          ( ( iType == IType::plus   && ( -a1 / b1 < -a2 / b2 ) ) ||
                            ( iType == IType::minus  && ( -a1 / b1 > -a2 / b2 ) ) ) );
       continue;
     }

     rxs.push_back( rx );
     rxSources.push_back( make_tuple( iType, it1, it2 ) );
  }
}

void LP2D::evalParams()
{
  mAlpha        = numeric_limits<double>::lowest();
  mBeta         = numeric_limits<double>::max();
  get<0>( mS )  = numeric_limits<double>::max();
  get<1>( mS )  = numeric_limits<double>::lowest();
  get<0>( mT )  = numeric_limits<double>::max();
  get<1>( mT )  = numeric_limits<double>::lowest();

  // evaluate alpha y
  for( const Constraint &constraint : mIn )
  {
     double a = get<0>( constraint );
     double b = get<1>( constraint );
     double c = get<2>( constraint );

     mAlpha = max( mAlpha, ( -a * mXm + c ) / b );
  }
  // end evaluate alpha y

  // evaluate beta y
  for( const Constraint &constraint : mIp )
  {
     double a = get<0>( constraint );
     double b = get<1>( constraint );
     double c = get<2>( constraint );

     mBeta = min( mBeta, ( -a * mXm + c ) / b );
  }
  // end evaluate beta y

  // evaluate s_max, s_min
  for( const Constraint &constraint : mIn )
  {
     double a     = get<0>( constraint );
     double b     = get<1>( constraint );
     double c     = get<2>( constraint );
     double slope = -a / b;

     if( a * mXm + b * mAlpha != c ) continue;

     get<0>( mS ) = min( get<0>( mS ), slope );
     get<1>( mS ) = max( get<1>( mS ), slope );
  }
  // end evaluate s_max, s_min

  // evaluate t_max, t_min
  for( const Constraint &constraint : mIn )
  {
     double a     = get<0>( constraint );
     double b     = get<1>( constraint );
     double c     = get<2>( constraint );
     double slope = -a / b;

     if( a * mXm + b * mAlpha != c ) continue;

     get<0>( mT ) = min( get<0>( mT ), slope );
     get<1>( mT ) = max( get<1>( mT ), slope );
  }
  // end evaluate t_max, t_min
}

LP2D::RelativePosition LP2D::evalOptPosition()
{
  if( mAlpha <= mBeta )
  {
    if( get<1>( mS ) < 0 ) return RelativePosition::right;
    if( get<0>( mS ) > 0 ) return RelativePosition::left;
    return RelativePosition::equal;
  }
  else
  {
    if( get<1>( mS ) < get<0>( mT ) ) return RelativePosition::right;
    if( get<0>( mS ) > get<1>( mT ) ) return RelativePosition::left;
    return RelativePosition::noSolution;
  }
}

/*
 *  a1 * x + b1 * y = c1
 *  a2 * x + b2 * y = c2
 *
 *  => y = ( a1 * c2 - a2 * c1 ) / ( a1 * b2 - a2 * b1 )
 */
double LP2D::solveReduced()
{
  assert( mIn.empty() || mIn.size() + mIp.size() <= mSizeSmall ); // precondition

  if( mIn.size() < mSizeSmall ) return numeric_limits<double>::lowest();

  Iterator  it1 = mIn.begin();
  Iterator  it2 = ++mIn.begin();
  double    a1  = get<0>( *it1 );
  double    b1  = get<1>( *it1 );
  double    c1  = get<2>( *it1 );
  double    a2  = get<0>( *it2 );
  double    b2  = get<1>( *it2 );
  double    c2  = get<2>( *it2 );

  return ( a1 * c2 - a2 * c1 ) / ( a1 * b2 - a2 * b1 );
}

void LP2D::removeConstraint( ConstraintList &constraints, Iterator &it1, Iterator &it2, bool removeIt1 )
{
  Iterator it;

  if( removeIt1 )
  {
    it = it1;
    --it1;
  }
  else
  {
    it = it2;
    --it1;
    --it2;
  }
  constraints.erase( it );
}
// end private member functions
