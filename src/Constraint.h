#ifndef CONSTRAINT_H
#define CONSTRAINT_H

#include <tuple>

namespace LinearProgramming
{

using Data = double;

/*!
  Linear programming constraint data structure.

  Represent constraint a * x + b * y <= c.
*/
class Constraint : public std::tuple<double,double,double>
{
  public:

    Constraint() = default;
    Constraint( double a, double b, double c );

    inline double a() const;
    inline double b() const;
    inline double c() const;

    inline void setA( double a );
    inline void setB( double b );
    inline void setC( double c );

    inline double y( double x ) const;
    inline double x( double y ) const;

    inline double yIntercept() const;
    inline double xIntercept() const;

    inline double slope() const;
};

inline bool   isParallel( const Constraint &c1, const Constraint &c2 );
inline double intersectX( const Constraint &c1, const Constraint &c2 );
inline double intersectY( const Constraint &c1, const Constraint &c2 );

inline double Constraint::a() const { return std::get<0>( *this ); }
inline double Constraint::b() const { return std::get<1>( *this ); }
inline double Constraint::c() const { return std::get<2>( *this ); }

inline void Constraint::setA( double a ) { std::get<0>( *this ) = a; }
inline void Constraint::setB( double b ) { std::get<1>( *this ) = b; }
inline void Constraint::setC( double c ) { std::get<2>( *this ) = c; }

inline double Constraint::y( double x ) const { return ( c() - a() * x ) / b(); }
inline double Constraint::x( double y ) const { return ( c() - b() * y ) / a(); }

inline double Constraint::yIntercept() const { return c() / b(); }
inline double Constraint::xIntercept() const { return c() / a(); }

inline double Constraint::slope() const { return -a() / b(); }

inline bool isParallel( const Constraint &c1, const Constraint &c2 )
{ return ( c1.a() == c2.a() && c1.b() == c2.b() ); }

/*!
 *  a1 * x + b1 * y = c1
 *  a2 * x + b2 * y = c2
 *
 *  => x = ( b2 * c1 - b1 * c2 ) / ( a1 * b2 - a2 * b1 )
 */
inline double intersectX( const Constraint &c1, const Constraint &c2 )
{ return ( c2.b() * c1.c() - c1.b() * c2.c() ) / ( c1.a() * c2.b() - c2.a() * c1.b() ); }

/*!
 *  a1 * x + b1 * y = c1
 *  a2 * x + b2 * y = c2
 *
 *  => y = ( a1 * c2 - a2 * c1 ) / ( a1 * b2 - a2 * b1 )
 */
inline double intersectY( const Constraint &c1, const Constraint &c2 )
{ return ( c1.a() * c2.c() - c2.a() * c1.c() ) / ( c1.a() * c2.b() - c2.a() * c1.b() ); }

}

#endif