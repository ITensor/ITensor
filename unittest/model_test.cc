#include "test.h"
#include <boost/test/unit_test.hpp>
#include "model/spinhalf.h"
#include "model/spinone.h"
#include "model/hubbard.h"
#include "model/spinless.h"
#include "model/tj.h"

using namespace itensor;

struct ModelTestDefaults
    {
    static const int N = 10;

    ModelTestDefaults()
        {
        }

    ~ModelTestDefaults() { }

    };

BOOST_FIXTURE_TEST_SUITE(ModelTest,ModelTestDefaults)

TEST(SpinHalfModel)
    {
    SpinHalf model(N);

    model.op("Sz",2); 
    model.op("S+",2); 
    model.op("S-",2); 
    model.op("Sp",2); 
    model.op("Sm",2); 
    model.op("Sx",2); 
    model.op("Sy",2); 
    model.op("ISy",2); 
    }

TEST(SpinOneModel)
    {
    SpinOne model(N);

    model.op("Sz",2); 
    model.op("S+",2); 
    model.op("S-",2); 
    model.op("Sp",2); 
    model.op("Sm",2); 
    model.op("Sx",2); 
    //model.op("Sy",2); 
    model.op("ISy",2); 
    }

TEST(HubbardModel)
    {
    Hubbard model(N);

    model.op("Nup",2); 
    model.op("Ndn",2); 
    model.op("Nupdn",2); 
    model.op("Ntot",2); 
    model.op("Sz",2); 
    model.op("Cup",2); 
    model.op("Cdn",2); 
    model.op("Aup",2); 
    model.op("Adn",2); 
    model.op("F",2); 
    }

TEST(SpinlessModel)
    {
    Spinless model(N);

    model.op("N",2); 
    model.op("A",2); 
    model.op("Adag",2); 
    model.op("F",2); 
    }

TEST(tJModel)
    {
    tJ model(N);

    model.op("Nup",2); 
    model.op("Ndn",2); 
    model.op("Aup",2); 
    model.op("Adn",2); 
    model.op("F",2); 
    }


BOOST_AUTO_TEST_SUITE_END()
