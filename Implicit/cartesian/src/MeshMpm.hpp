// *********************************************************************************
// File: MeshMpm.hpp
//
// Details: Container for the MPM nodes, elements and sub-meshes
//          Has containers with pointers to the nodes, elements and pointers
//
// Dependency: MPM_ITEMS contains all template parameters
//
// Author: Krishna Kumar and Samila Bandara, University of Cambridge
//
// Version: 1.0
// *********************************************************************************

#ifndef CARTESIAN_SRC_MESHMPM_H_
#define CARTESIAN_SRC_MESHMPM_H_

#include <vector>
#include <string>
#include <iostream>
#include <algorithm>
#include <boost/bind.hpp>

#include "cartesian/src/MpmItems.hpp"
#include "corlib/Mesh.hpp"
#include "cartesian/src/Enums.hpp"

namespace mpm {
template<typename MPM_ITEMS>
class MeshMpm;
}

// \brief container for the MPM nodes, elements and sub-meshes
// \details Has containers with pointers to the nodes, elements and pointers
// \tparam MPM_ITEMS  Contains all the template parameters

template<typename MPM_ITEMS>
class mpm::MeshMpm : public corlib::Mesh<typename MPM_ITEMS::NodeHandle, typename MPM_ITEMS::ElementHandle> {
public:
    typedef MPM_ITEMS  MpmItems;

    static const unsigned dim = MpmItems::dim;

    typedef typename MpmItems::NodeHandle NodeType;
    typedef typename MpmItems::ElementHandle ElementType;
    typedef typename MpmItems::ParticleHandle  ParticleType;
    // typedef typename MpmItems::ParticleSoilHandle ParticleSoilType;
    // typedef typename MpmItems::ParticleWaterHandle ParticleWaterType;
    typedef typename MpmItems::ParticleCloudSoilHandle ParticleCloudSoilType;
    // typedef typename MpmItems::ParticleCloudWaterHandle ParticleCloudWaterType;
    typedef typename MpmItems::SubMeshHandle SubMeshType;
    typedef typename MpmItems::MatrixHandle MatrixType;

    typedef typename MpmItems::NodePtr NodePtr;
    typedef typename MpmItems::ElementPtr ElementPtr;
    typedef typename MpmItems::ParticlePtr ParticlePtr;
    typedef typename MpmItems::ParticleSoilPtr ParticleSoilPtr;
    // typedef typename MpmItems::ParticleWaterPtr ParticleWaterPtr;
    typedef typename MpmItems::ParticleCloudSoilPtr ParticleCloudSoilPtr;
    // typedef typename MpmItems::ParticleCloudWaterPtr ParticleCloudWaterPtr;
    typedef typename MpmItems::SubMeshPtr  SubMeshPtr;
    typedef typename MpmItems::MatrixPtr MatrixPtr;


protected:
    typedef corlib::Mesh<NodeType, ElementType> Mesh_;
    typedef std::vector<NodePtr> NodeVec_;
    typedef std::vector<ElementPtr> ElementVec_;
    typedef std::vector<SubMeshPtr> SubMeshVec_;
    typedef std::vector<ParticleSoilPtr> ParticleSoilVec_;
    typedef std::vector<ParticlePtr> ParticleVec_;
    // typedef std::vector<ParticleWaterPtr> ParticleWaterVec_;



public:
    typedef typename NodeVec_::const_iterator NodeConstIterator;
    typedef typename ElementVec_::const_iterator ElementConstIterator;
    typedef typename SubMeshVec_::const_iterator SubMeshConstIterator;

    typedef mpm::CompareFunctor<ElementPtr> ElementCompare;
    typedef std::set<ElementPtr, ElementCompare> ElementSet;
    typedef mpm::CompareFunctor<NodePtr> NodeCompare;
    typedef std::set<NodePtr, NodeCompare> NodeSet;

    typedef ublas::bounded_vector<double, dim> VecDim;

protected:
    typedef typename ElementSet::iterator ElementSetIterator_;
    typedef typename NodeSet::iterator NodeSetIterator_;
//-----------------------------------------------------------------------------

public:
// Constructor with SMF file, & submesh information file
    MeshMpm(std::istream& smf, std::istream& subms);

// Destructor which clears the containers
    ~MeshMpm();

//------------------------------------------------------------------------------
                       // Nodal constraints functions

// Read constraints from stream for explicit MPM only
    void readConstraintsExplicit(std::istream& cstr);

// Read constraints from stream for explicit MPM for mixed mesh with sloped boundary
    void readConstraintsMixedMeshExplicit(std::istream& cstr);
    
// Set the constraints given access to tuples
    // template<typename IT> void setConstraints(IT iter, IT end);

//------------------------------------------------------------------------------

                  // locate particles, elements and nodes

// locate particles and elements (templated by soil/water particle cloud type)
    void locateParticlesAndElements(ParticleCloudSoilPtr pCloudSoilPtr);

// CHECK:print the particles in elements
    void printParticlesInElement();

// set of nodes connected to each single node which has paricles
    void addNodesToNodes();

// template <typename P_CLOUD_TYPE>
    // void locateParticlesAndElements(P_CLOUD_TYPE* pCloudType);

// add pointer to elements and nodes which contain particles
    void addElementAndNodesOfSoilParticle(ElementPtr ePtr, NodeVec_ nVec) {
        elementsOfSoilPSet_.insert(ePtr);
        for (unsigned i = 0; i < nVec.size(); i++)
            nodesOfSoilPSet_.insert(nVec[i]);
    }

//-----------------------------------------------------------------------------

// intialise data inside mesh, elements, nodes
    void initialise();

//-----------------------------------------------------------------------------
        // Generic application of functors to node and element containers

// Iterate over elements which contains soil particles
    template<typename OP> OP iterateOverElementsOfSoilP(OP op);

// Iterate over nodes which contains soil particles
    template<typename OP> OP iterateOverNodesOfSoilP(OP op);

// Iterate over submeshes
    template<typename OP> OP iterateOverSubMeshes(OP op);

// Iterate over submeshes with predicate
    template<typename OP, typename PRED>
    OP iterateOverSubMeshesWithPredicate(OP op, PRED p);

//-----------------------------------------------------------------------------
             // Const access to node and element containers

// Give const access to begin of subMesh container
    SubMeshConstIterator SubMeshesBegin() const {
        return subMeshes_.begin();
    }

// Give const access to end of subMesh container
    SubMeshConstIterator SubMeshesEnd() const {
        return subMeshes_.end();
    }

// Give const access to specific SubMesh
    SubMeshPtr SubMeshIterator(const unsigned numSubMesh) const;

//-----------------------------------------------------------------------------
                          // return data

// Give element pointer
    ElementPtr giveElement(unsigned id) {
        return elements_[id];
    }

// Give pointers to nodes which contain particles
    NodeSet giveNodesOfParticles() {
        return nodesOfSoilPSet_;
    }

// Give pointers to elements which contain particles
    ElementSet giveElementsOfParticles() {
        return elementsOfSoilPSet_;
    }

//------------------------------------------------------------------------------
                      // Size queries

// Return number of elements
    unsigned numSubMeshes() const {
        return subMeshes_.size();
    }

//------------------------------------------------------------------------------
protected:
    using corlib::Mesh<NodeType, ElementType>::nodes_;    // Vector with node pointers
    using corlib::Mesh<NodeType, ElementType>::elements_; // Vector with element pointers
    SubMeshVec_ subMeshes_;                               // Vector with sub mesh pointers
    ElementSet elementsOfSoilPSet_;                       // Elements which contains soil particles
    NodeSet nodesOfSoilPSet_;                             // nodes which contains soil particles
    // ParticleSoilVec_     soilParticles_;               // Vector with soil particle pointers
    // ParticleWaterVec_     waterParticles_;             // Vector with water particle pointers
    // ParticleCloudSoilPtr  pSoilCloudPtr_;
    // ParticleCloudWaterPtr pWaterCloudPtr_;
};

#include "MeshMpm.ipp"

#endif  // CARTESIAN_SRC_MESHMPM_H_
