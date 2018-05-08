//
// Distributed under the ITensor Library License, Version 1.2
//    (See accompanying LICENSE file.)
//
#ifndef __ITENSOR_TEVOLOBSERVER_H
#define __ITENSOR_TEVOLOBSERVER_H
#include "itensor/util/readwrite.h"
#include "itensor/mps/observer.h"

namespace itensor {

//
// Class for monitoring time evolution calculations.
// The measure and checkDone methods are virtual
// so that behavior can be customized in a
// derived class.
//

class TEvolObserver : public Observer
    {
    public:
    
    TEvolObserver(const Args& args = Global::args());

    virtual ~TEvolObserver() { }

    void virtual
    measure(const Args& args = Global::args());
    
    bool virtual
    checkDone(const Args& args = Global::args());

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
TEvolObserver(const Args& args) 
    : 
    done_(false),
    show_percent_(args.getBool("ShowPercent",true))
    { 
    }


void inline TEvolObserver::
measure(const Args& args)
    {
    const Real t = args.getReal("Time");
    if(show_percent_)
        {
        const Real ttotal = args.getReal("TotalTime");
        Real percentdone = (100.*t)/ttotal;
        if(percentdone < 99.5 || (std::fabs(t-ttotal) < 1E-10))
            {
            printf("\b\b\b%2.f%%",percentdone);
            std::cout.flush();
            }
        }
    }


bool inline TEvolObserver::
checkDone(const Args& args)
    {
    const Real t = args.getReal("Time");
    if(fileExists("STOP_TEVOL"))
        {
        println("File STOP_TEVOL found: stopping this time evolution run at time ",t);
        system("rm -f STOP_TEVOL");
        return true;
        }

    //Set done_ flag to true so any outer callers using this Observer will also terminate.
    if(fileExists("STOP_TEVOL_ALL"))
        {
        println("File STOP_TEVOL_ALL found: stopping this time evolution at time ",t);
        system("rm -f STOP_TEVOL_ALL");
        done_ = true;
        return done_;
        }
    
    return done_;
    }

} //namespace itensor


#endif // __ITENSOR_TEVOLOBSERVER_H
