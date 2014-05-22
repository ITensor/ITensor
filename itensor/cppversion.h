#ifndef __ITENSOR_CPPVERSION_H_
#define __ITENSOR_CPPVERSION_H_

#ifndef USE_CPP11

#include "boost/foreach.hpp"
#define Foreach BOOST_FOREACH

#include "boost/array.hpp"

namespace itensor {
using boost::array;
};

#include "boost/make_shared.hpp"
namespace itensor {
using boost::shared_ptr;
using boost::make_shared;
};


#else // USE_CPP11 is defined

#define Foreach(X,Y) for(X : Y)

#include <array>
namespace itensor {
using std::array;
using std::shared_ptr;
using std::make_shared;
};

#include <memory>
namespace itensor {
using std::shared_ptr;
using std::make_shared;
};

#endif

#endif
