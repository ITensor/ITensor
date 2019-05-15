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
#include "itensor/global.h"

namespace itensor {

bool&
Global::checkArrows()
    {
    static bool checkArrows_ = true;
    return checkArrows_;
    }
bool&
Global::debug1()
    {
    static bool debug1_ = false;
    return debug1_;
    }
bool&
Global::debug2()
    {
    static bool debug2_ = false;
    return debug2_;
    }
bool&
Global::debug3()
    {
    static bool debug3_ = false;
    return debug3_;
    }
bool&
Global::debug4()
    {
    static bool debug4_ = false;
    return debug4_;
    }
bool&
Global::printdat()
    {
    static bool printdat_ = false;
    return printdat_;
    }
Real&
Global::printScale()
    {
    static Real printScale_ = 1E-10;
    return printScale_;
    }
bool&
Global::showIDs()
    {
    static bool showIDs_ = true;
    return showIDs_;
    }
Real
Global::random(int seed)
    {
    using Generator = std::mt19937;
    using Distribution = std::uniform_real_distribution<Real>;

    static Generator rng(std::time(NULL)+getpid());
    static Distribution dist(0,1);

    if(seed != 0)  //reseed rng
        {
        rng = Generator(seed);
        }

    return dist(rng);
    }
void
Global::warnDeprecated(const std::string& message)
    {
    static int depcount = 1;
    if(depcount <= 10)
        {
        println("\n\n",message,"\n");
        ++depcount;
        }
    }
bool&
Global::read32BitIDs()
    {
    static bool read32_ = false;
    return read32_;
    }
} //namespace itensor
