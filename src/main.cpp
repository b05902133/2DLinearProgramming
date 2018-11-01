#include <iostream>
using namespace std;

#include "Select.h"

int main()
{
  Select engine;

  cout << engine.select( { 1, 4, 2, 3, 5, 7, 6, 8, 10, 9 }, 3 ) << "\n";

  return 0;
}
