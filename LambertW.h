#ifndef _utl_LambertW_h_
#define _utl_LambertW_h_

#include <limits>


namespace utl {

  template<int branch>
  double LambertW(const double x);


  inline
  double
  LambertW(const int branch, const double x)
  {
    switch (branch) {
    case -1: return LambertW<-1>(x);
    case  0: return LambertW<0>(x);
    default: return std::numeric_limits<double>::quiet_NaN();
    }
  }

}


#endif
