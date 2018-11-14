#include "Constraint.h"

namespace LinearProgramming
{

Constraint::Constraint( double a, double b, double c )
{
  setA( a );
  setB( b );
  setC( c );
}

}