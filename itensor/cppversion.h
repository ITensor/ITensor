#ifndef __ITENSOR_CPPVERSION_H_
#define __ITENSOR_CPPVERSION_H_

#ifdef USE_CPP11

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

#include <random>
namespace itensor {
using std::mt19937;
using std::uniform_real_distribution;
};

#include <functional>
namespace itensor {
using std::function;
};

#else // USE_CPP11 is not defined

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

#include "boost/random/mersenne_twister.hpp"
#include "boost/random/uniform_real_distribution.hpp"
namespace itensor {
using boost::random::mt19937;
using boost::random::uniform_real_distribution;
};

#include "boost/function.hpp"
namespace itensor {
using boost::function;
};

#endif

#endif
