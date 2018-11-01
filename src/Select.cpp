#include "Select.h"

#include <cmath>
#include <algorithm>
#include <cassert>
#include <limits>
#include <iostream>
using namespace std;

Select::Data Select::select( std::vector<Data> datas, std::size_t k )
{
  assert( 0 < k && k <= datas.size() ); // precondition

  // solution space is small
  if( datas.size() <= mR )
  {
    sort( datas.begin(), datas.end() );

    return datas[k-1];
  }
  // end solution space is small

  const int groupNum = ceil( static_cast<double>( datas.size() ) / mR );

  vector<Data> groupMediums;
  vector<Data> lessThanMediums;
  vector<Data> greaterThanMediums;

  groupMediums.reserve( groupNum );
  lessThanMediums.reserve( datas.size() / 2 + 1 );
  greaterThanMediums.reserve( datas.size() / 2 + 1 );

  while( datas.size() % mR != 0 )
    datas.push_back( numeric_limits<Data>::max() );

  for( int i = 0 ; i < groupNum ; ++i )
  {
     auto itBegin = datas.begin() + i * mR;
     auto itEnd   = itBegin + mR;

     sort( itBegin, itEnd );

     groupMediums.push_back( *( itBegin + mR / 2 ) );
  }

  Data medium = ( groupNum == 1 ) ? groupMediums.front() :
                                    select( groupMediums, groupMediums.size() / 2 );

  for( const Data &data : datas )
  {
     if     ( data < medium ) lessThanMediums.push_back( data );
     else if( data > medium ) greaterThanMediums.push_back( data );
  }

  if      ( k <= lessThanMediums.size()      )
    return select( lessThanMediums, k );
  else if ( k == lessThanMediums.size() + 1 )
    return medium;
  else
    return select( greaterThanMediums, k - lessThanMediums.size() - 1 );
}
