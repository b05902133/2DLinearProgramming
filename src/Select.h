#ifndef SELECT_H
#define SELECT_H

#include <vector>

class Select
{
  public:

    using Data = int;

    Data select( std::vector<Data> datas, std::size_t k );

  private:

    constexpr static int mR = 5;

    std::vector<Data> mDatas;
};

#endif
