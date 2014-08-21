#ifndef __ITENSOR_TYPES_H_
#define __ITENSOR_TYPES_H_

namespace itensor {

using Real = double;

};

#include <array>
namespace itensor {
using std::array;
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
using std::uniform_int_distribution;
};

#include <functional>
namespace itensor {
using std::function;
};

#endif
