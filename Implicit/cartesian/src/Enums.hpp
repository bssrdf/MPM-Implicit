// *********************************************************************************
// File: Enums.hpp
//
// Details: Enumerator for describing the particle type tags
//
// Author: Samila Bandara, University of Cambridge
//
// Version: 1.0
// *********************************************************************************

#ifndef CARTESIAN_SRC_ENUMS_H_
#define CARTESIAN_SRC_ENUMS_H_

#include <vector>
#include <string>
#include <iostream>
#include <algorithm>
#include "cartesian/src/MpmItems.hpp"

namespace mpm {
enum ParticleTypeTag {
    SOIL,
    WATER
};
}

#endif  // CARTESIAN_SRC_ENUMS_H_


