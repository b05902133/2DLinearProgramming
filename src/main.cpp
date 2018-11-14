#include <iostream>
#include <cmath>
#include <limits>
using namespace std;

#include "LP.h"
using namespace LinearProgramming;

int main()
{
  LP2D                solver;
  vector<Constraint>  constraints;
  int                 constraintNum;
  double              yMin;

  cin >> constraintNum;

  constraints.reserve( constraintNum );

  for( int i = 0 ; i < constraintNum ; ++i )
  {
     int a, b, c;

     cin >> a >> b >> c;

     constraints.push_back( Constraint( a, b, c ) );
  }

  yMin = round( solver.solve( constraints ) );

  if( yMin == numeric_limits<double>::max() ) // nosolution
    cout << "NA\n";
  else if( yMin == numeric_limits<double>::lowest() ) // -infinity
    cout << "-INF\n";
  else
    cout << yMin << "\n";

  return 0;
}
