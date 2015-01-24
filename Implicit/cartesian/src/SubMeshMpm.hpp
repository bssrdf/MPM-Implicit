// *********************************************************************************
// File: SubMeshMpm.hpp
//
// Details: Base class for sub-meshes
//
// Dependency: MPM_ITEMS contains all template parameters
//
// Author: Krishna Kumar and Samila Bandara, University of Cambridge
//
// Version: 1.0
// *********************************************************************************

#ifndef CARTESIAN_SRC_SUBMESHMPM_H_
#define CARTESIAN_SRC_SUBMESHMPM_H_

#include <algorithm>
#include <iostream>
#include <set>
#include <string>
#include <vector>

#include "cartesian/src/MpmItems.hpp"

namespace mpm {
template<typename MPM_ITEMS>
class SubMeshMpm;
namespace ublas = boost::numeric::ublas;

}


// \brief container for sub-meshes
// \details Has containers with pointers to the nodes and the elements
// \tparam MPM_ITEMS   the main template file which includes all the classes

template<typename MPM_ITEMS>
class mpm::SubMeshMpm {
public:
    typedef MPM_ITEMS MpmItems;

    static const unsigned dim = MpmItems::dim;
    static const unsigned numNodes = MpmItems::nNodes;

    typedef typename MpmItems::NodeHandle NodeType;
    typedef typename MpmItems::ElementHandle ElementType;
    typedef typename MpmItems::ParticleHandle ParticleType;
    typedef typename MpmItems::SubMeshHandle SubMeshType;
    typedef typename MpmItems::MeshHandle MeshType;

    typedef typename MpmItems::NodePtr NodePtr;
    typedef typename MpmItems::ElementPtr ElementPtr;
    typedef typename MpmItems::ParticlePtr ParticlePtr;
    typedef typename MpmItems::ParticleSoilPtr ParticleSoilPtr;
    typedef typename MpmItems::SubMeshPtr SubMeshPtr;
    typedef typename MpmItems::MeshPtr MeshPtr;

    typedef mpm::CompareFunctor<ElementPtr> ElementCompare;
    typedef std::set<ElementPtr, ElementCompare> ElementSet;
    typedef mpm::CompareFunctor<ParticleSoilPtr> ParticleCompare;
    typedef std::set<ParticleSoilPtr, ParticleCompare> ParticleSet;

    typedef ublas::bounded_vector<unsigned, numNodes> VecNN;
    typedef ublas::bounded_vector<double, dim> VecDim;
    typedef ublas::bounded_vector<unsigned, dim> VecDimUnsign;

protected:
    typedef std::vector<NodePtr> NodeVec_;
    typedef std::vector<ElementPtr> ElementVec_;
    typedef typename ElementSet::iterator ElementItr_;

public:
    SubMeshMpm(MeshPtr mPtr) : myMeshPtr_(mPtr) {
        elementSet_.clear();
        elementsOfParticlesSet_.clear();
        particleSet_.clear();
    }

    ~SubMeshMpm() {
        elementSet_.clear();
        elementsOfParticlesSet_.clear();
        particleSet_.clear();
        myMeshPtr_ = NULL;
    }

    // Read sub mesh details
    std::istream& readSubMesh(std::istream& subms);

    // Give slope angle for sloped sub-mesh
    const double giveSlopeAngle() {
        return slopeAngle_;
    }

    // add particles to elements and submesh, add element to particle
    bool addParticleToSubMesh(ParticleSoilPtr pPtr);

    // Give element set that contains particles
    ElementSet giveElementsOfParticles() {
        return elementsOfParticlesSet_;
    }

    // Give particle set in the submesh
    ParticleSet giveParticles() {
        return particleSet_;
    }

    void printCheckSMesh() {
        ElementItr_ eItr;
        for (eItr = elementsOfParticlesSet_.begin(); eItr != elementsOfParticlesSet_.end(); eItr++)
            (*eItr) -> printCheckElem();
    }

    // write particle data (id, 2, position, velocity, stress, strain)
    void writeParticleData(std::ostream& out);

private:
    // check whether a particle is inside the submesh
    bool isParticleIncluded(ParticleSoilPtr pPtr);

    // assign particle and element pointers to correspoding locations
    void setParticle(ParticleSoilPtr pPtr);

    // add pointers to all the elements inside sub-mesh
    void addElement(ElementPtr ePtr) {
        elementSet_.insert(ePtr);
    }

    // add pointer to an element where each particle is located
    void addElementOfParticle(ElementPtr ePtr) {
        elementsOfParticlesSet_.insert(ePtr);
    }

    // add pointers of particles inside submesh into a set
    void addParticle(ParticleSoilPtr pPtr) {
        particleSet_.insert(pPtr);
    }

protected:
    MeshPtr myMeshPtr_;                  // pointer to mesh
    unsigned meshType_;                  // type of the mesh (0- horizontal/ 1- slope)
    double slopeAngle_;                  // slope angle w.r.t. 0 direction
    unsigned elemIdStart_;               // id of first element in the sub mesh
    VecDimUnsign numElem_;               // number of elements in the submesh
    VecDim dxElem_;                      // lengths of element
    VecDim nFirstESCoord_;               // coordinates of first node of submesh
    VecDim nLastEECoord_;                // coordinates of last node of submesh
    ElementSet elementSet_;              // elements inside sub-mesh
    ElementSet elementsOfParticlesSet_;  // elements contain particles
    ParticleSet particleSet_;            // particle pointer inside submesh
};

#include "SubMeshMpm.ipp"


#endif  // CARTESIAN_SRC_SUBMESHMPM_H_
