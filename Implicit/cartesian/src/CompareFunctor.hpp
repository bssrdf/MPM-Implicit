// *********************************************************************************
// File: CompareFunctor.hpp
//
// Details: Comparison of functor for pointers
//
// Dependency: NULL
//
// Author: Samila Bandara and Kosala Bandara, University of Cambridge
//
// Version / Revision: 1.0 (2009)
// *********************************************************************************

#ifndef CARTESIAN_SRC_COMPAREFUNCTOR_H_
#define CARTESIAN_SRC_COMPAREFUNCTOR_H_

namespace mpm {
template<typename T>
class CompareFunctor;
}

// comparison of functor for pointers
template<typename T>
class mpm::CompareFunctor {
public:
    bool operator()(const T& lhs, const T& rhs) const {
        return &(*lhs) < &(*rhs);
    }
};

#endif  // CARTESIAN_SRC_COMPAREFUNCTOR_H_

