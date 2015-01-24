// *********************************************************************************
// File: ParticleMpm.hpp
//
// Details: Base class for particles
//
// Dependency: MPM_ITEMS contains all template parameters
//
// Author: Krishna Kumar and Samila Bandara, University of Cambridge
//
// Version: 1.0
// *********************************************************************************

#ifndef CARTESIAN_SRC_PARTICLEMPM_H_
#define CARTESIAN_SRC_PARTICLEMPM_H_

#include <algorithm>
#include <cassert>
#include <iostream>
#include <set>
#include <utility>
#include <vector>

#include <boost/array.hpp>
#include <boost/numeric/ublas/vector.hpp>

#include "cartesian/src/MpmItems.hpp"
#include "cartesian/src/CompareFunctor.hpp"

namespace mpm {
template<typename MPM_ITEMS>
class ParticleMpm;
namespace ublas = boost::numeric::ublas;
}


// \brief Base class for particles
// \details Includes basic calculations for particles
// \tparam MPM_ITEMS   the main template file which includes all the classes

template<typename MPM_ITEMS>
class mpm::ParticleMpm {
public:
    typedef MPM_ITEMS MpmItems;
    static const unsigned dof = MpmItems::dof;
    static const unsigned dim = MpmItems::dim;
    static const unsigned numNodes = MpmItems::nNodes;

    typedef typename MpmItems::NodeHandle NodeType;
    typedef typename MpmItems::ElementHandle ElementType;
    typedef typename MpmItems::ParticleHandle ParticleType;

    typedef typename MpmItems::NodePtr NodePtr;
    typedef typename MpmItems::ElementPtr ElementPtr;
    typedef typename MpmItems::ParticlePtr ParticlePtr;

    typedef ublas::bounded_vector<double, dim> VecDim;
    typedef ublas::bounded_vector<double, numNodes> VecNN;
    typedef ublas::bounded_matrix<double, dim, numNodes> MatDimNN;
    typedef ublas::bounded_matrix<double, dim, dim> MatDimDim;

protected:
    typedef std::vector<ElementPtr> ElementVec_;
    typedef std::vector<NodePtr> NodeVec_;

public:
    // Constructor with id number
    ParticleMpm(const unsigned& id) : id_(id) {
        coord_.clear();
        ePtrOfP_ = NULL;
        meshType_ = 0;
    }

    // Constructor with id and coordinates
    ParticleMpm(const unsigned& id, const VecDim& cv) : id_(id), coord_(cv) { }
    // ~ParticleMpm () { dofIndices_.clear();  particleForces_.clear(); }  // ?

private:
    typedef std::vector<std::pair<unsigned, double>> ParticleForceVec_;
    typedef boost::array<unsigned, dof>  IndexArray_;

public:
    // Set coordinates \details
    // \param[in] c  Vector containing the coordinates
    void setParticleCoordinates(const VecDim& c) {
        coord_ = c;
    }

    // Give the coordinates \details
    // \retval coord_ The coordinates of this particle
    VecDim giveParticleCoordinates() const {
        return coord_;
    }

 //-----------------------------------------------------------------------------
                            // Give Shape functions 

   // Give shape functions of particle
    VecNN giveParicleSfun() const {
        return sfunc_;
    }

    // Give shape function of particle for a given node
    double giveParticleSfunAtNode(unsigned& id) {
        return sfunc_(id);
    }

//-----------------------------------------------------------------------------
                          // Give gradient shape functions 

   // Give gradient shape functions of particle
    MatDimNN giveParicleGradShepeFunctions() const {
        return dsfuncDX_;
    }

    double giveParticleGradSfunAtNode(unsigned i, unsigned& j) {
        return dsfuncDX_(i,j);
    }

//------------------------------------------------------------------------------

    // Add element pointer and node pointers vector relevant to particle location
    void addElementAndNodes(ElementPtr ePtr, NodeVec_ nVec);

    // Give element pointer which contains the particle
    ElementPtr giveElement() {
        return ePtrOfP_;
    }

    // WRITE: Store particle external forces
    void storeParticleForce(const unsigned& component, const double& value);

    // WRITE: return constraints
    void giveParticleForces(ParticleForceVec_& GlobalParticleForces) const;

    // add mesh type and slope angle if it is a sloped mesh
    void addMeshTypeAndSlope(const unsigned& meshType, const double& slopeAngle);

    // Insert coordinates to given insert iterator
    template<typename IT> IT passParticleCoordinates(IT backIter) const;

    // Return particle's id
    unsigned giveParticleId() const {
        return id_;
    }

    // evaluate shape functions and gradients of shape functions for horizontal mesh
    void evaluateShpfunAndGradShpfunHorMesh();

    // evaluate shape functions and gradients of shape functions for mixed mesh
    void evaluateShpfunAndGradShpfunMixedMesh();

    // clear pointers of nodes and elements
    void clearMeshDetails();

    // NOT USING: Compute stiffness component from particle
    void computeStiffness(const VecDim& xi); // THIS SHOULD BE IN DERIVED CLASS

    // Set Material Type for Particle
    void setMaterial(const unsigned materialID) {
        matID_ = materialID;
    }

    // Write id to a stream
    std::ostream& writeParticleId(std::ostream& out) const;
    // Write coordinates to a stream
    std::ostream& writeParticleCoordinates(std::ostream& out) const;
    // Read coordinates from a stream
    std::istream& readParticleSelf(std::istream& inp);

protected:
    // send sign
    double sign_(double var);

    // Evaluation of shape function and gradient
    // Calculate local coordinate for each particle location
    void localCoordinates_();

    // Calculate local coordinate for each particle location for sloped mesh
    void localCoordinatesSlope_();

    // Return values of shape function for nodes
    void sfun_();

    // Return gradient of shape function at xi
    void sfunGrad_(MatDimNN& dsfuncDxi);

    //  Compute global derivatives for horizontal / slope (2D/3D)
    void globalDerivatives_();

    // Compute simple global derivatives only for regular mesh (horizontal)
    void globalDerivativesHorMesh_();

    //  Compute global derivatives for any mesh type (higher computing time)
    void globalDerivativesExtended_();

    // Give inverse of Jacobi matrix for horizontal / slope mesh
    void jacobianInverse_(MatDimDim& Ji);

    // Give Jacobi matrix for horizontal / slope mesh
    void jacobian_(MatDimDim& J);

protected:
    unsigned id_;                       // identifiying number
    unsigned matID_;                    // Material ID
    VecDim coord_;                      // coordinates
    IndexArray_ dofIndices_;            // indices of the degrees of freedom
    unsigned meshType_;                 // type of the mesh (0- horizontal/ 1- slope)
    double  slopeAngle_;                // slope angle w.r.t. 0 direction
    ElementPtr ePtrOfP_;                // pointer to the element which has particle
    NodeVec_ nodesOfParticle_;          // nodes of element in which particle is in
    VecDim xi_;                         // local coordinates
    VecNN  sfunc_;                      // shape functions
    MatDimNN dsfuncDX_;                 // gradient shape function
    MatDimNN dsfuncDXCenter_;           // gradient shape fucntion at the center of the element
    ParticleForceVec_ particleForces_;  // <dof,value>-pairs for prescribed dofs
};

#include "ParticleMpm.ipp"

#endif  // CARTESIAN_SRC_PARTICLEMPM_H_
