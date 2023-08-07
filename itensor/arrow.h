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
#ifndef __ITENSOR_ARROW_H
#define __ITENSOR_ARROW_H

#include <iostream>
#include "util/error.h"

namespace itensor {

/*
*
* The Arrow enum is used to label how indices
* transform under a particular symmetry group. 
* Indices with an Out Arrow transform as vectors
* (kets) and with an In Arrow as dual vectors (bras).
*
* Conventions regarding arrows:
*
* * Arrows point In or Out, never right/left/up/down.
*
* * The Site indices of an MPS representing a ket point Out.
*
* * Conjugation switches arrow directions.
*
* * All arrows flow Out from the ortho center of an MPS 
*   (assuming it's a ket - In if it's a bra).
*
* * IQMPOs are created with the same arrow structure as if they are 
*   orthogonalized to site 1, but this is just a default since they 
*   aren't actually ortho. If position is called on an IQMPO it follows 
*   the same convention as for an MPS except Site indices point In and 
*   Site' indices point Out.
*
* * Local site operators have two IQIndices, one unprimed and pointing In, 
*   the other primed and pointing Out.
*
*/

enum Arrow { In = -1, Out = 1, Neither = 0 };

Arrow inline
toArrow(int i)
    {
    int In_int = static_cast<int>(In);
    if(In_int == i) return In;
    return Out;
    }

Arrow inline
operator-(Arrow dir)
    {
#ifdef DEBUG
    if(dir == Neither)
        error("Cannot reverse Arrow direction 'Neither'");
#endif
    return (dir == In ? Out : In);
    }

inline std::ostream& 
operator<<(std::ostream& s, Arrow D)
    { 
    switch(D)
        {
        case In:
            s << "In";
            return s;
        case Out:
            s << "Out";
            return s;
        case Neither:
            s << "Neither";
            return s;
        default:
            error("Missing Arrow case");
        }
    return s; 
    }

struct ArrowError : ITError
    {
    ArrowError(const std::string& message) 
        : ITError(message)
        { }
    };

} // namespace itensor

#endif

