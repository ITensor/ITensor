//
// Copyright 2018 The Simons Foundation, Inc. - All Rights Reserved.
//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
//    http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.
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
    
    TEvolObserver(Args const& args = Args::global());

    virtual ~TEvolObserver() { }

    void virtual
    measure(Args const& args = Args::global());
    
    bool virtual
    checkDone(Args const& args = Args::global());

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
        std::remove("STOP_TEVOL");
        return true;
        }

    //Set done_ flag to true so any outer callers using this Observer will also terminate.
    if(fileExists("STOP_TEVOL_ALL"))
        {
        println("File STOP_TEVOL_ALL found: stopping this time evolution at time ",t);
        std::remove("STOP_TEVOL_ALL");
        done_ = true;
        return done_;
        }
    
    return done_;
    }

} //namespace itensor


#endif // __ITENSOR_TEVOLOBSERVER_H
