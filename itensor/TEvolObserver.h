//
// Distributed under the ITensor Library License, Version 1.0.
//    (See accompanying LICENSE file.)
//
#ifndef __ITENSOR_TEVOLOBSERVER_H
#define __ITENSOR_TEVOLOBSERVER_H
#include "observer.h"

#define Cout std::cout
#define Endl std::endl
#define Format boost::format

//
// Class for monitoring time evolution calculations.
// The measure and checkDone methods are virtual
// so that behavior can be customized in a
// derived class.
//

class TEvolObserver : public Observer
    {
    public:
    
    TEvolObserver(const OptSet& opts = Global::opts());

    virtual ~TEvolObserver() { }

    void virtual
    measure(const OptSet& opts = Global::opts());
    
    bool virtual
    checkDone(const OptSet& opts = Global::opts());

    private:

    /////////////
    //
    // Data Members

    bool done_,
         show_percent_;

    //
    /////////////

    }; // class TEvolObserver

inline TEvolObserver::
TEvolObserver(const OptSet& opts) 
    : 
    done_(false),
    show_percent_(opts.getBool("ShowPercent",true))
    { 
    }


void inline TEvolObserver::
measure(const OptSet& opts)
    {
    const Real t = opts.getReal("Time");
    if(show_percent_)
        {
        const Real ttotal = opts.getReal("TotalTime");
        Real percentdone = (100.*t)/ttotal;
        if(percentdone < 99.5 || (fabs(t-ttotal) < 1E-10))
            {
            Cout << Format("\b\b\b%2.f%%") % percentdone;
            Cout.flush();
            }
        }
    }


bool inline TEvolObserver::
checkDone(const OptSet& opts)
    {
    const Real t = opts.getReal("Time");
    if(fileExists("STOP_TEVOL"))
        {
        Cout << "File STOP_TEVOL found: stopping this time evolution run at time " << t << Endl;
        system("rm -f STOP_TEVOL");
        return true;
        }

    //Set done_ flag to true so any outer callers using this Observer will also terminate.
    if(fileExists("STOP_TEVOL_ALL"))
        {
        Cout << "File STOP_TEVOL_ALL found: stopping this time evolution at time " << t << Endl;
        system("rm -f STOP_TEVOL_ALL");
        done_ = true;
        return done_;
        }
    
    return done_;
    }

#undef Cout
#undef Endl
#undef Format

#endif // __ITENSOR_TEVOLOBSERVER_H
