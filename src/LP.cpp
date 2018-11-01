#include "LP.h"

#include <limits>
#include <cassert>
#include <algorithm>
using namespace std;

// public member functions
double LP::solve( const std::vector<Constraint> &constraints )
{
  init( constraints );

  return iterate();
}
// end public member functions

// private member functions
void LP::init( const std::vector<Constraint> &constraints )
{
  mXl = numeric_limits<double>::min();
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

double LP::iterate()
{
}
// end private member functions
