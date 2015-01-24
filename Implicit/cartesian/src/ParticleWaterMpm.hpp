// *********************************************************************************
// File: ParticleWaterMpm.hpp
//
// Details: Base class for Water particles
//
// Dependency: MPM_ITEMS contains all template parameters
//
// Author: Samila Bandara, University of Cambridge
//
// Version: 1.0
// *********************************************************************************

#ifndef CARTESIAN_SRC_PARTICLEWATERMPM_H_
#define CARTESIAN_SRC_PARTICLEWATERMPM_H_

#include <cassert>
#include <set>
#include <iostream>

#include <boost/array.hpp>
#include <boost/numeric/ublas/vector.hpp>

#include "cartesian/src/MpmItems.hpp"
#include "cartesian/src/ParticleMpm.hpp"
#include "cartesian/src/CompareFunctor.hpp"
#include "cartesian/src/Enums.hpp"

namespace mpm {
template<typename MPM_ITEMS>
class ParticleWaterMpm;
namespace ublas = boost::numeric::ublas;
}


// \brief class for water particles
// \details Includes basic calculations for particles
// \tparam MPM_ITEMS   the main template file which includes all the classes

template<typename MPM_ITEMS>
class mpm::ParticleWaterMpm : public mpm::ParticleMpm<MPM_ITEMS> {
public:
    typedef MPM_ITEMS MpmItems;
    static const unsigned dof = MpmItems::dof;
    static const unsigned dim = MpmItems::dim;
    static const unsigned numNodes = MpmItems::nNodes;

    static const mpm::ParticleTypeTag myType = mpm::WATER;

    typedef typename MpmItems::NodeHandle NodeType;
    typedef typename MpmItems::ElementHandle ElementType;
    typedef typename MpmItems::ParticleHandle ParticleType;

    typedef typename MpmItems::MaterialHandle MaterialType;

    typedef typename MpmItems::NodePtr NodePtr;
    typedef typename MpmItems::ElementPtr ElementPtr;
    typedef typename MpmItems::ParticlePtr ParticlePtr;

    typedef typename MpmItems::MaterialPtr MaterialPtr;

    typedef ublas::bounded_vector<double, MpmItems::dim> VecDim;
    typedef ublas::bounded_vector<double, numNodes > VecNV;
    typedef ublas::bounded_matrix<double, MpmItems::dim, numNodes> MatDimNV;
    typedef ublas::bounded_matrix<double, MpmItems::dim, MpmItems::dim> MatDimDim;

protected:
    typedef std::vector<ElementPtr> ElementVec_;
    // typedef mpm::CompareFunctor<ElementPtr> ElementCompare;
    // typedef std::set<ElementPtr, ElementCompare> ElementSet;

public:
    // Constructor with id number
    ParticleWaterMpm(const unsigned  & id) : ParticleType(id) {}

    // Constructor with id and coordinates
    ParticleWaterMpm(const unsigned & id, const VecDim & cv) : ParticleType(id, cv) {}

    // ~ParticleMpm () { dofIndices_.clear();  particleForces_.clear(); } // ???


private:
    // typedef std::vector<std::pair<unsigned,double>>           ParticleForceVec_;
    // typedef boost::array<unsigned, dof>                       IndexArray_;


public:
    // cache material to particle and initialise material properties of particle
    // void cacheMaterial(MaterialPtr  m);
    void cacheMaterial(MaterialPtr  m) {
        material_ = m;
        std::cout << "water" << "\n";
        return;
    }

    // Calculate mass of particle
    double massParticle() {
        return mass_;
    }

protected:
    //  calculate mass of particle (constant value)
    void computeMass_();

protected:
    unsigned id_;              // identifiying number
    double mass_;              // mass of particle (constant)
    double density_;           // density of particle (solid density)
    double porosity_;          // porosity of particle
    MaterialPtr material_;     // Pointer to material object
};

#include "ParticleWaterMpm.ipp"

#endif  // CARTESIAN_SRC_PARTICLEWATERMPM_H_
