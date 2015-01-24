// *********************************************************************************
// File: MeshMpm.ipp
//
// Details: Container for the MPM nodes, elements and sub-meshes
//          Has containers with pointers to the nodes, elements and pointers
//
// Dependency: MPM_ITEMS contains all template parameters
//
// Author: Krishna Kumar, Samila Bandara | University of Cambridge
//
// Version: 1.0
// *********************************************************************************

#include <boost/tuple/tuple.hpp>

// Constructor given a stream to a simple mesh file. This file must have the
// following format:

//     NN   NE \n
//     x1  y1  z1  \n
//     ....        \n
//     xNN yNN zNN \n
//     v1_1  v1_2 ... v1_NPE   \n
//     ...                     \n
//     vNE_1 vNE_2 ... vNE_NPE \n

// where NN is the number of nodes, NE the number of elements, and NPE refers to
// the number of nodes per element. Note that the latter number is fixed at
// compile time. (Refer the image)
// \image html numbering.jpg
// \param smf Input stream

template<typename MPM_ITEMS>
mpm::MeshMpm<MPM_ITEMS>::MeshMpm(std::istream& smf, std::istream& subms) :Mesh_(smf) {

    this->iterateOverElements(boost::bind(&ElementType::computeCtrCoordAndLength, _1));
    // for xi calculations( for constant mesh only, not applicable if mesh changes with time )

    unsigned numSubMeshes;

    subms >> numSubMeshes;
    subMeshes_.reserve(numSubMeshes);

    for (unsigned s = 0; s < numSubMeshes; s ++) {
        SubMeshPtr theSubMesh = new SubMeshType(this);
        theSubMesh->readSubMesh(subms);
        subMeshes_.push_back(theSubMesh);
    }

    elementsOfSoilPSet_.clear();
    nodesOfSoilPSet_.clear();

    // typedef ublas::bounded_matrix< double, 3, 3 >    Mat3X3_;
    // Mat3X3_ AA;
    // AA(0,0) = -1.3027756; AA(0,1) = 1.732; AA(0,2) = 0.;
    // AA(1,0) = 1.732; AA(1,1) = -2.30277; AA(1,2) = 0.;
    // AA(2,0) = 0.; AA(2,1) = 0.; AA(2,2) = -1.3027756;
    // const double val = corlib::inverse( AA );
    // std::cout<<"AA"<<AA<<"\n";

    return;
}

//-----------------------------------------------------------------------------
                     // Free all dynamically allocated memory

template<typename MPM_ITEMS>
mpm::MeshMpm<MPM_ITEMS>::~MeshMpm() {
    this->iterateOverSubMeshes(corlib::deleteFunctor());
    subMeshes_.clear();
    return;
}

//-----------------------------------------------------------------------------
                 // Intialise data inside mesh, elements, nodes

template<typename MPM_ITEMS>
void mpm::MeshMpm<MPM_ITEMS>::initialise() {
    elementsOfSoilPSet_.clear();
    nodesOfSoilPSet_.clear();
    this->iterateOverNodes(boost::bind(&NodeType::initialise, _1));
    this->iterateOverElements(boost::bind(&ElementType::initialise, _1));
    return;
}

//-----------------------------------------------------------------------------
                // locate particles and elements
                // assign particle pointers into corresponding element
                // assign element pointer into corresponding particle
                // only applicable for soil particles

template<typename MPM_ITEMS>
void mpm::MeshMpm<MPM_ITEMS>::locateParticlesAndElements(ParticleCloudSoilPtr pSoilCloudPtr) {
    unsigned numSoilParticles = pSoilCloudPtr->numParticles();
    bool status[numSoilParticles];

    #pragma omp parallel sections num_threads(2)
    {
        #pragma omp section
        {
            for (unsigned i = 0; i < numSoilParticles; i ++) {
                status[i] = 0;
                for (unsigned j = 0; j < subMeshes_.size(); j ++) {
                    ParticleSoilPtr pSoilPtr = pSoilCloudPtr->giveParticlePtr(i);
                    SubMeshPtr sMeshPtr = subMeshes_[j];
                    if (sMeshPtr->addParticleToSubMesh(pSoilPtr)) {
                        status[i] = 1;
                        break;
                    }
                }
                FSI_VERIFY(status[i]);
            }
        }
    }

    #pragma omp barrier
    return;
}

//------------------------------------------------------------------------------
/*
                 // print the coordinates of particels in each element

template<typename MPM_ITEMS>
void mpm::MeshMpm<MPM_ITEMS>::printParticlesInElement() {

    ParticleVec_ particlesInElement;
    unsigned i = 0;
//    std::cout << "Number of elements: " << elementsOfSoilPSet_.size() << "\n";

    for (ElementSetIterator_ iterE = elementsOfSoilPSet_.begin(); iterE != elementsOfSoilPSet_.end(); ++iterE) {
        ElementPtr ptrE = *iterE;
        particlesInElement = ptrE->giveParticles();
        std::cout << "Element " << i << ": " << "\n";
        for (unsigned j = 0; j < particlesInElement.size(); ++j) {
            ParticlePtr ptrP = particlesInElement.at(j);
            VecDim coord = ptrP->giveParticleCoordinates();
            std::cout << coord << "\n";
        }
        i++;
        std::cout << "\n";
    }

}
*/
//-----------------------------------------------------------------------------
               // add elements (only which contains particles) connected 
               // to each node (only which containes particles)

template<typename MPM_ITEMS>
void mpm::MeshMpm<MPM_ITEMS>::addNodesToNodes() {

    for (ElementSetIterator_ iterE = elementsOfSoilPSet_.begin(); iterE != elementsOfSoilPSet_.end(); ++iterE) {
        ElementPtr ptrE = *iterE;
        NodeVec_ nodesOfElem = ptrE->giveNodes();
        for (unsigned i = 0; i < nodesOfElem.size(); ++i) {
            NodePtr ptrN = nodesOfElem.at(i);
            ptrN->addElement(ptrE);
            for (unsigned j = 0; j < nodesOfElem.size(); ++j) {
                ptrN->addNode(nodesOfElem.at(j));
            }
        }
    }
/*
// print the connected nodes of each nodes only which contain particles
    for (NodeSetIterator_ iterN = nodesOfSoilPSet_.begin(); iterN != nodesOfSoilPSet_.end(); ++iterN) {
        NodePtr ptrNode = *iterN;
        unsigned id = ptrNode->giveId();
        NodeSet setOfNodes = ptrNode->giveNodeSet();
        std::cout << "Node " << id << ": Connected Nodes" << std::endl;
     for (NodeSetIterator_ iterNode = setOfNodes.begin(); iterNode != setOfNodes.end(); ++iterNode) {
         NodePtr nPtr = *iterNode;
         unsigned nId = nPtr->giveId();
         std::cout << nId << " ";
     }
     std::cout << std::endl;
    }
*/
}

//-----------------------------------------------------------------------------
            // locate particles (soil/water) and elements
            // assign particle pointers into corresponding element
            // assign element pointer into corresponding particle
            // check what is the particle type and use the corresponding vectors (SOIL/WATER)

// template<typename MPM_ITEMS>
// template<typename P_CLOUD_TYPE>
// void mpm::MeshMpm<MPM_ITEMS>::locateParticlesAndElements(P_CLOUD_TYPE* pCloudPtr)
// {
//     if(pCloudPtr->myType == mpm::SOIL) {
// 	pSoilCloudPtr_ = dynamic_cast<ParticleCloudSoilPtr>(pCloudPtr); //avoid type casting as much as possible
// 	soilParticles_ =  pSoilCloudPtr_->giveParticles();
// 	for (unsigned i = 0; i < soilParticles_.size(); i++) {
// 	    ParticleSoilPtr pSoilPtr = soilParticles_[i];
// 	    for (unsigned j = 0; j < subMeshes_.size(); j++) {
// 		SubMeshPtr sMeshPtr = subMeshes_[j];
// 		if (sMeshPtr->addParticleToSubMesh(pSoilPtr)) break;
// 		std::cout << "soil particle ID= "<< i << " , submesh ID= " << j << ", particle not added " <<  std::endl;
// 	    }
// 	}
//     }else if(pCloudPtr->myType == mpm::WATER) {
// 	pWaterCloudPtr_ = dynamic_cast<ParticleCloudWaterPtr>(pCloudPtr);
//     	waterParticles_ =  pWaterCloudPtr_->giveParticles();
//     	for (unsigned i = 0; i < waterParticles_.size(); i++) {
//     	    ParticleWaterPtr pWaterPtr = waterParticles_[i];
//     	    for (unsigned j = 0; j < subMeshes_.size(); j++) {
//     		SubMeshPtr sMeshPtr = subMeshes_[j];
//     		if (sMeshPtr -> addParticleToSubMesh(pWaterPtr)) break;
//     		std::cout << "water particle ID= "<<i << " , submesh ID= " <<j<<", particle not added "<<  std::endl;
//     	    }
//     	}
//     }

// }

//------------------------------------------------------------------------------
                  // Apply a functor to all elements that contains soil particles using std::for_each
                // \tparam OP         Type of the element operation
                // \param[in,out] op  Specific operation applied to all elements
                // \retval OP         for_each returns the operator

template<typename MPM_ITEMS>
template<typename OP>
OP mpm::MeshMpm<MPM_ITEMS>::iterateOverElementsOfSoilP(OP op) {
    ElementSetIterator_ begin = elementsOfSoilPSet_.begin();
    ElementSetIterator_ end   = elementsOfSoilPSet_.end();
    return std::for_each(begin, end, op);
}

//------------------------------------------------------------------------------
              // Apply a functor to all nodes inside elements which contains soil particles
              // using std::for_each
              // \tparam OP         Type of the node operation
              // \param[in,out] op  Specific operation applied to all nodes
              // \retval OP         for_each returns the operator

template<typename MPM_ITEMS>
template<typename OP>
OP mpm::MeshMpm<MPM_ITEMS>::iterateOverNodesOfSoilP(OP op) {
    NodeSetIterator_ begin = nodesOfSoilPSet_.begin();
    NodeSetIterator_ end   = nodesOfSoilPSet_.end();
    return std::for_each(begin, end, op);
}

//------------------------------------------------------------------------------
               // Apply a given functor to all sub-meshes making use of std::for_each
               // \tparam OP         Type of the submesh operation
               // \param[in,out] op  Specific operation applied to all submeshes
               // \retval OP         for_each returns the operator

template<typename MPM_ITEMS>
template<typename OP>
OP mpm::MeshMpm<MPM_ITEMS>::iterateOverSubMeshes(OP op) {
    typename SubMeshVec_::iterator begin = subMeshes_.begin();
    typename SubMeshVec_::iterator end   = subMeshes_.end();
    return std::for_each(begin, end, op);
}

//------------------------------------------------------------------------------
// Apply a given functor only to those subMeshes which have a given predicate.
// Given the predicate 'p', the functor 'op' will only be applied to those
// subMeshes for which holds 'p(X) = true'.
// \tparam OP          Type of functor to be applied to subMeshes
// \tparam PRED        Type of predicate functor
// \param[in,out] op   The functor to be applied
// \param[in] p        The predicate
// \retval op          Returns the given functor

template<typename MPM_ITEMS>
template<typename OP, typename PRED>
OP mpm::MeshMpm<MPM_ITEMS>::iterateOverSubMeshesWithPredicate(OP op, PRED p) {
    typename SubMeshVec_::iterator begin = subMeshes_.begin();
    typename SubMeshVec_::iterator end   = subMeshes_.end();
    return corlib::for_each_if(begin, end, op, p);
}

//------------------------------------------------------------------------------
               // Given an number numNode, return the correponding node pointer
               // \param[in] numNode  Number of the requested node
               // \retval    NodePtr_ Pointer to the requested node

template<typename MPM_ITEMS>
typename MPM_ITEMS::SubMeshPtr mpm::MeshMpm<MPM_ITEMS>::SubMeshIterator(const unsigned numSubMesh) const {
    return subMeshes_.at(numSubMesh);
}

//------------------------------------------------------------------------------
              // Read from a stream the nodal constraints. 
              // The file must have the format

//   NGC     NFC    \n
//   n_1     d_1    \n
//   ...            \n
//   n_NGC   d_NGC  \n
//   n_1     dn_1     sdn_1     \n
//   ...                        \n
//   n_NFC   dn_NFC   sdn_NFC   \n

//   where NGC is the number of general constraints in this file, which are used
//   to set velocity and acceleration components in the d_i direction to
//   be zero. n_i is the node number, d_i is the direction number (0|1|2)
//   NFC is the number of friction constraints, which are used to limit the
//   acceleration in the tangential direction to boundary.
//   n_i is the node number, dn_i is the perpendicular direction to boundary* (0|1|2)
//   sdn_i is the sign (+1/-1) of dn_i w.r.t coordinate system ( i.e. =-1 if bottom of
//   mesh boundary ).
//   \param[in] cstr Input stream containing the constraints

template<typename MPM_ITEMS>
void mpm::MeshMpm<MPM_ITEMS>::readConstraintsExplicit(std::istream& cstr) {
    // get number of constraints
    unsigned numGenConstraints, numFrictionConstraints;
    cstr >> numGenConstraints >> numFrictionConstraints;
    // go through general constraints
    for (unsigned c = 0; c < numGenConstraints; c++) {
        unsigned nodenum, direction;
        cstr >> nodenum >> direction;
        nodes_.at(nodenum)->storeGeneralConstraint(direction);
    }
    // go through friction constraints
    for (unsigned c = 0; c < numFrictionConstraints; c++) {
        unsigned nodenum, normDirection;
        int signNormDirection;
        cstr >> nodenum >> normDirection >> signNormDirection;
        nodes_.at(nodenum)->storeFrictionConstraint(normDirection, signNormDirection);
    }
}

//------------------------------------------------------------------------------
                  // Read from a stream the nodal constraints for a mesh
                  //  with mixed mesh ( sloped and horizontal boundaries ).
                  //   The file must have the format:

//   NGC       NFC        NGCS        NFCS  \n
//   n_1       d_1        \n
//   ...                  \n
//   n_NGC     d_NGC      \n
//   n_1       dn_1       sdn_1       \n
//   ...                              \n
//   n_NFC     dn_NFC     sdn_NFC     \n
//   n_1       sm_1       \n
//   ...                  \n
//   ns_NGCS   sm_NGCS    \n
//   n_1       sm_1       sdns_1      \n
//   ...                              \n
//   n_NFCS    sm_NFCS    sdns_NFCS   \n

//   where NGC is the number of general constraints, NFC is the number of friction
//   constraints, NGCS is the number of general constraints in slope boundary and
//   NFCS is the number of friction constraints in slope boundary. General constraints
//   are used to set velocity and acceleration normal to the boundary to zero. Friction
//   constraints are used to limit acceleration tangential to boundary according to
//   Coulomb's friction law.

//   n_i is the node number, d_i is the direction number (0|1|2).
//   dn_i is the perpendicular direction to boundary* (0|1|2), sdn_i is the sign
//   (+1/-1) of dn_i w.r.t coordinate system ( i.e. =-1 if bottom of mesh boundary ).
//   sm_i is the id of submesh which includes the node. sdns_i is the sign of the
//   normal force w.r.t. boundary ( i.e. =-1 if bottom of the boundary ).
//   \param[in] cstr Input stream containing the constraints

template<typename MPM_ITEMS>
void mpm::MeshMpm<MPM_ITEMS>::readConstraintsMixedMeshExplicit(std::istream& cstr) {
    // get number of constraints
    unsigned numGenCons, numFrictionCon, numGenConsSlope, numFrictionConSlope;
    cstr >> numGenCons >> numFrictionCon >> numGenConsSlope >> numFrictionConSlope;

    // go through general constraints for 0/1/2 directions
    for (unsigned c = 0; c < numGenCons; c++) {
        unsigned nodenum, direction;
        cstr >> nodenum >> direction;
        nodes_.at(nodenum) -> storeGeneralConstraint(direction);
    }

    // go through friction constraints for 0/1/2 directions
    for (unsigned c = 0; c < numFrictionCon; c++) {
        unsigned nodenum, normDirection;
        int signNormDirection;
        cstr >> nodenum >> normDirection >> signNormDirection;
        nodes_.at(nodenum)->storeFrictionConstraint(normDirection, signNormDirection);
    }

    // go through general constraints for slope boundaries
    for (unsigned c = 0; c < numGenConsSlope; c++) {
        unsigned nodenum, subMeshId;
        double slopeAngle;
        cstr >> nodenum >> subMeshId;
        slopeAngle = subMeshes_[subMeshId]->giveSlopeAngle();
        nodes_.at(nodenum)->storeGeneralConstraintSlopeBn(slopeAngle);
    }

    // go through friction constraints for slope boundaries
    for (unsigned c = 0; c < numFrictionConSlope; c++) {
        unsigned nodenum, subMeshId;
        int signNormDirection;
        double slopeAngle;
        cstr >> nodenum >> subMeshId >> signNormDirection;
        slopeAngle = subMeshes_[subMeshId] -> giveSlopeAngle();
        nodes_.at(nodenum) ->
        storeFrictionConstraintSlopeBn(slopeAngle, signNormDirection);
    }
}
//------------------------------------------------------------------------------
     // * Set the nodal constraints given iterators to a container of constraint
//   triples: <n,d,v> with the node number \e n, the constraint direction \e d,
//   and the value of the constraint \e v.
//   \tparam IT       Type of iterator to triple container
//   \param[in] iter  Begin of the triple container
//   \param[in] end   End   of the triple container

// template<typename MPM_ITEMS>
// template<typename IT>
// void mpm::MeshMpm<MPM_ITEMS>::setConstraints(IT iter, IT end) {
//     for ( ; iter != end; ++ iter ) {
//       boost::tuple<unsigned,unsigned,double> t = *iter;
//       const unsigned nodenum   = t.get<0>();
//       const unsigned component = t.get<1>();
//       const double   value     = t.get<2>();
//       nodes_.at(nodenum)->storeConstraint(component, value);
//     }
// }

//------------------------------------------------------------------------------

