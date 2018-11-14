#include "LP.h"

#include <limits>
#include <cassert>
#include <algorithm>
#include <vector>
#include <cmath>

#ifndef NDEBUG
#include <iostream>
#endif
using namespace std;

#include "Select.h"

namespace LinearProgramming
{

// public member functions
/*!
 *  \return
 *
 *    - numeric_limits<double>::min():  -infinity
 *    - numeric_limits<double>::max():  no solution
 *    - else value:                     minimun y value
 */
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
     if( constraint.b() == 0 )
     {
       double x = constraint.xIntercept();

       if( x < mXl || x > mXr ) continue;

       if( constraint.a() > 0 ) // x <= abc( c / a )
         mXr = min( mXr, x );
       else if( constraint.a() < 0 ) // x >= abc( c / a )
         mXl = max( mXl, x );
     }
     else if( constraint.b() > 0 )
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
#ifndef NDEBUG
    clog << "I+ size: " << mIp.size() << "\t I- size:" << mIn.size() << "\n";
#endif

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

      return get<1>( rxSources[i] )->y( mXm );
    }

    if( optPosition == RelativePosition::left )
    {
      for( size_t i = 0 ; i < rxs.size() ; ++i )
      {
         if( rxs[i] < mXm ) continue;

         Iterator it1 = get<1>( rxSources[i] );
         Iterator it2 = get<2>( rxSources[i] );

         if( get<0>( rxSources[i] ) == IType::plus )
           removeConstraint( mIp, it1, it2, it1->slope() < it2->slope() );
         else
           removeConstraint( mIn, it1, it2, it1->slope() > it2->slope() );
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

         if( get<0>( rxSources[i] ) == IType::plus )
           removeConstraint( mIp, it1, it2, it1->slope() > it2->slope() );
         else
           removeConstraint( mIn, it1, it2, it1->slope() < it2->slope() );
      }
      mXl = mXm;
    }
    if( mXr < mXl ) return numeric_limits<double>::max();
  }
  return solveReduced();
}

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
     double rx;

     if( isParallel( *it1, *it2 ) ) // parallel
     {
       double yShift1 = it1->yIntercept();
       double yShift2 = it2->yIntercept();

       removeConstraint(  I, it1, it2,
                          ( ( iType == IType::plus   && yShift1 > yShift2 ) ||
                            ( iType == IType::minus  && yShift1 < yShift2 ) ) );
       continue;
     }

     rx = intersectX( *it1, *it2 );

     if( rx < mXl )
     {
       removeConstraint(  I, it1, it2,
                          ( ( iType == IType::plus   && ( it1->slope() > it2->slope() ) ) ||
                            ( iType == IType::minus  && ( it1->slope() < it2->slope() ) ) ) );
       continue;
     }

     if( rx > mXr )
     {
       removeConstraint(  I, it1, it2,
                          ( ( iType == IType::plus   && ( it1->slope() < it2->slope() ) ) ||
                            ( iType == IType::minus  && ( it1->slope() > it2->slope() ) ) ) );
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
     mAlpha = max( mAlpha, constraint.y( mXm ) );
  // end evaluate alpha y

  // evaluate beta y
  for( const Constraint &constraint : mIp )
     mBeta = min( mBeta, constraint.y( mXm ) );
  // end evaluate beta y

  // evaluate s_max, s_min
  for( const Constraint &constraint : mIn )
  {
     double slope = constraint.slope();

     if( round( constraint.a() * mXm + constraint.b() * mAlpha ) != constraint.c() ) continue;

     get<0>( mS ) = min( get<0>( mS ), slope );
     get<1>( mS ) = max( get<1>( mS ), slope );
  }
  // end evaluate s_max, s_min

  // evaluate t_max, t_min
  for( const Constraint &constraint : mIp )
  {
     double slope = constraint.slope();

     if( round( constraint.a() * mXm + constraint.b() * mBeta ) != constraint.c() ) continue;

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

double LP2D::solveReduced()
{
  assert( mIn.empty() || mIn.size() + mIp.size() <= mSizeSmall ); // precondition

#ifndef NDEBUG
  clog << "solve reduced probem\n";
#endif

  if( mIn.empty() ) return numeric_limits<double>::lowest();

  ConstraintList constraints;

  constraints.insert( constraints.end(), mIp.begin(), mIp.end() );
  constraints.insert( constraints.end(), mIn.begin(), mIn.end() );

  Iterator it1 = constraints.begin();
  Iterator it2 = ++constraints.begin();

  // make sure slope of it1 > slope of it2
  if( it1->slope() < it2->slope() )
    swap( it1, it2 );

  double slope1 = it1->slope();
  double slope2 = it2->slope();

  if( slope1 == slope2 )
  {
    double yShift1 = it1->yIntercept();
    double yShift2 = it2->yIntercept();

    if      ( it1->b() > 0 && it2->b() < 0 )
    {
      if( yShift2 > yShift1 ) return numeric_limits<double>::max();

      return ( slope2 > 0 ) ? it2->y( mXl ) : it2->y( mXr );
    }
    else if ( it1->b() < 0 && it2->b() > 0 )
    {
      if( yShift1 > yShift2 ) return numeric_limits<double>::max();

      return ( slope1 > 0 ) ? it1->y( mXl ): it1->y( mXr );
    }
    else // b1 < 0 && b2 < 0
    {
      if( yShift1 > yShift2 )
        return ( slope1 > 0 ) ? it1->y( mXl ): it1->y( mXr );
      else
        return ( slope2 > 0 ) ? it2->y( mXl ): it2->y( mXr );
    }
  }

  if      ( slope2 > 0 )  // all slope > 0
  {
    if      ( it1->b() > 0 && it2->b() < 0 )  // quadrant I
    {
      double x = intersectX( *it1, *it2 );

      if( x > mXr ) return numeric_limits<double>::max();
      if( x < mXl ) return it2->y( mXl );

      return intersectY( *it1, *it2 );
    }
    else if ( it1->b() < 0 && it2->b() > 0 )  // quadrant III
    {
      if( mXl != numeric_limits<double>::lowest() ) return it1->y( mXl );

      return numeric_limits<double>::lowest();
    }
    else                          // quadrant II
    {
      if( mXl != numeric_limits<double>::lowest() ) return it2->y( mXl );

      return numeric_limits<double>::lowest();
    }
  }
  else if ( slope1 < 0 )  // all slope < 0
  {
    if      ( it1->b() < 0 && it2->b() > 0 )  // quadrant II
    {
      double x = intersectX( *it1, *it2 );

      if( x < mXl ) return numeric_limits<double>::max();
      if( x > mXr ) return it1->y( mXr );

      return intersectY( *it1, *it2 );
    }
    else if ( it1->b() > 0 && it2->b() < 0 )  // quadrant IV
    {
      if( mXr != numeric_limits<double>::max() ) return it2->y( mXr );

      return numeric_limits<double>::lowest();
    }
    else                          // quadrant III
    {
      if( mXr != numeric_limits<double>::max() ) return it1->y( mXr );

      return numeric_limits<double>::lowest();
    }
  }
  else
  {
    if      ( it1->b() > 0 && it2->b() < 0 )  // right
    {
      if( mXr != numeric_limits<double>::max() ) return it2->y( mXr );

      return numeric_limits<double>::lowest();
    }
    else if ( it1->b() < 0 && it2->b() > 0 )  // left
    {
      if( mXl != numeric_limits<double>::lowest() ) return it1->y( mXl );

      return numeric_limits<double>::lowest();
    }
    else                          // up
    {
      double x = intersectX( *it1, *it2 );

      if( x < mXl ) return it1->y( mXl );
      if( x > mXr ) return it2->y( mXr );

      return intersectY( *it1, *it2 );
    }
  }
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

}