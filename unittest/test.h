#ifndef __TEST_H
#define __TEST_H

#define CHECK_CLOSE BOOST_CHECK_CLOSE 
#define CHECK_EQUAL BOOST_CHECK_EQUAL 
#define CHECK       BOOST_CHECK

#define TEST   BOOST_AUTO_TEST_CASE

#include <iostream>

inline double 
Func(double x)
    {
    return x*x;
    }

inline int 
Func(int x)
    {
    return x*x;
    }

class Functor
    {
    public:

    double
    operator()(double x) const
        {
        return x*x;
        }

    int
    operator()(int x) const
        {
        return x*x;
        }

    };

#endif
