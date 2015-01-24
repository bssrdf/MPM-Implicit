//-----------------------------------------------------------------------------
                   // free all dynamically allocate memory

template<typename MPM_ITEMS>
mpm::SolverMpm<MPM_ITEMS>::~SolverMpm() {

    delete theMatrix;
}

//-----------------------------------------------------------------------------
                              // initialise
template<typename MPM_ITEMS>
void mpm::SolverMpm<MPM_ITEMS>::initialise() {

    nodesOfParticles_.clear();
    theMatrix->initialise();
}

//-----------------------------------------------------------------------------
                          // Matrix computation

template<typename MPM_ITEMS>
void mpm::SolverMpm<MPM_ITEMS>::evaluateMatrix(MeshPtr& mesh, double& dt) {

    this->setPtrsToNodesAndElementsOfP(mesh);
// 
    unsigned locationOfM_Matrix = 0;
    unsigned locationOfK_Matrix = 0;
    unsigned numOfActiveNodes = nodesOfParticles_.size();

    for (unsigned i = 0; i < nodesOfParticles_.size(); ++i) {
        this->evaluateNonZeroEntries(nodesOfParticles_.at(i),i);
    }
    // theMatrix->printnonZeroColumnPositions();

    theMatrix->resizeMatrix(numOfActiveNodes);

    for (unsigned i = 0; i < nodesOfParticles_.size(); ++i) {
        this->computeNonZeroM_Entries(nodesOfParticles_.at(i), locationOfM_Matrix);
        this->computeNonZeroK_Entries(nodesOfParticles_.at(i), locationOfK_Matrix, dt);
        this->computeNonZeroF_Entries(nodesOfParticles_.at(i), dt, i);
    }

    // theMatrix->printMmatrix();
    // theMatrix->printKmatrix();
    // theMatrix->printFmatrix();


}

//-----------------------------------------------------------------------------
                     // evaluate the "nodesOfParticles_" vector
           // which contains the pointers to all nodes which contain particles

template<typename MPM_ITEMS>
void mpm::SolverMpm<MPM_ITEMS>::setPtrsToNodesAndElementsOfP(MeshPtr ptrMesh) {

    NodeSet nodesOfP_ = ptrMesh->giveNodesOfParticles();
    ElementSet elementsOfP_ = ptrMesh->giveElementsOfParticles();
    for (NodeSetIterator_ iterN = nodesOfP_.begin(); iterN != nodesOfP_.end(); ++iterN) {
        NodePtr ptrNode = *iterN;
        nodesOfParticles_.push_back(ptrNode);
        // unsigned id = ptrNode->giveId();
        // std::cout << id << "\n";
    }

    for (ElementSetIterator_ iterE = elementsOfP_.begin(); iterE != elementsOfP_.end(); ++iterE) {
        ElementPtr ptrElem = *iterE;
        elementsOfParticles_.push_back(ptrElem);
        // unsigned id = ptrElem->giveElementId();
        // std::cout << id << "\n";
    }
}

//------------------------------------------------------------------------------
                     // evaluate the vecor ("columnMatrix" in Matrix file)
           // which contains the column position of nonzero entries in matrices

template<typename MPM_ITEMS>
void mpm::SolverMpm<MPM_ITEMS>::evaluateNonZeroEntries(NodePtr rowId, unsigned row) {

    NodeSet nodesOfNode = rowId->giveNodeSet();
    for (NodeSetIterator_ iterN = nodesOfNode.begin(); iterN != nodesOfNode.end(); ++iterN) {
        NodePtr ptrN = *iterN;
        unsigned id = ptrN->giveId();
        // std::cout << id << " ";
        for (unsigned i = 0; i < nodesOfParticles_.size(); ++i) {
            if (id == nodesOfParticles_.at(i)->giveId()) {
                theMatrix->nonZeroColumnPositions(i);
                theMatrix->nonZeroRowPositions(row);
                break;
            } 
        }
    }
    // std::cout << std::endl;
}

//-----------------------------------------------------------------------------
                     // evaluate the nonzero entries in M matrix

template<typename MPM_ITEMS>
void mpm::SolverMpm<MPM_ITEMS>::computeNonZeroM_Entries(NodePtr rowId, unsigned& n) {

    ElementSet elementsOfNode = rowId->giveElements();
    NodeSet nodesOfNode = rowId->giveNodeSet();
    unsigned num = 0;

    for (NodeSetIterator_ iterN = nodesOfNode.begin(); iterN != nodesOfNode.end(); ++iterN) {
         NodePtr nPtr = *iterN;
         double MEntry = 0;
         double rho = 1000.0;
         // double massOfP = rho*0.4*0.4*0.4;


        for (ElementSetIterator_ iterE = elementsOfNode.begin(); iterE != elementsOfNode.end(); ++iterE) {
             ElementPtr ptrE = *iterE;
             unsigned nodeRow = ptrE->giveNodeLocation(rowId);
             // std::cout << "Node Row: " << nodeRow << " ";
             ParticleSoilVec_ particlesOfE = ptrE->giveParticles();

            if (ptrE->isNodeIncluded(nPtr)) {
                unsigned nodeCol = ptrE->giveNodeLocation(nPtr);
                // std::cout << "Node Col: " << nodeCol<< std::endl;
                for (unsigned i = 0; i < particlesOfE.size(); ++i) {
                    double massOfP = particlesOfE.at(i)->massParticle();
                    double sFunRow = particlesOfE.at(i)->giveParticleSfunAtNode(nodeRow);
                    double sFunCol = particlesOfE.at(i)->giveParticleSfunAtNode(nodeCol);
                    MEntry += massOfP*sFunRow*sFunCol;
                }
            }
        }
        unsigned position = num + n;
        // std::cout << "position: "<<  position << std::endl;
        theMatrix->setNonZeroMEntries(MEntry,(position));
         ++num;
         // std::cout << "MEntry: " << MEntry << std::endl;
         // std::cout << std::endl;
    }
    n += nodesOfNode.size();
    // std::cout << "n: " << n << std::endl;

}

//------------------------------------------------------------------------------
                    // evaluate the nonzero entries in all K matrices

template<typename MPM_ITEMS>
void mpm::SolverMpm<MPM_ITEMS>::computeNonZeroK_Entries(NodePtr rowId_, unsigned& n_, double dt_) {

    ElementSet elementsOfNode_ = rowId_->giveElements();
    NodeSet nodesOfNode_ = rowId_->giveNodeSet();
    unsigned num_ = 0;

    for (NodeSetIterator_ iterN = nodesOfNode_.begin(); iterN != nodesOfNode_.end(); ++iterN) {
         NodePtr nPtr_ = *iterN;
         double Kxx_Entry = 0;
         double Kyy_Entry = 0;
         double Kxy_Entry = 0;
         double mu_ = 0.001;
         double rho_ = 1000.;
         // double massOfP_ = rho_*0.4*0.4*0.4;

        for (ElementSetIterator_ iterE = elementsOfNode_.begin(); iterE != elementsOfNode_.end(); ++iterE) {
             ElementPtr ptrE_ = *iterE;
             unsigned nodeRow_ = ptrE_->giveNodeLocation(rowId_);
             // std::cout << "Node Row: " << nodeRow << " ";
             ParticleSoilVec_ particlesOfE_ = ptrE_->giveParticles();

            if (ptrE_->isNodeIncluded(nPtr_)) {
                unsigned nodeCol_ = ptrE_->giveNodeLocation(nPtr_);
                // std::cout << "Node Col: " << nodeCol<< std::endl;
                for (unsigned i = 0; i < particlesOfE_.size(); ++i) {
                    double massOfP_ = particlesOfE_.at(i)->massParticle();
                    double dXsFunRow = particlesOfE_.at(i)->giveParticleGradSfunAtNode(0,nodeRow_);
                    double dYsFunRow = particlesOfE_.at(i)->giveParticleGradSfunAtNode(1,nodeRow_);
                    double dXsFunCol = particlesOfE_.at(i)->giveParticleGradSfunAtNode(0,nodeCol_);
                    double dYsFunCol = particlesOfE_.at(i)->giveParticleGradSfunAtNode(1,nodeCol_);
                    Kxx_Entry += (massOfP_*mu_*dt_/rho_)*(2*dXsFunRow*dXsFunCol + dYsFunRow*dYsFunCol);
                    Kyy_Entry += (massOfP_*mu_*dt_/rho_)*(dXsFunRow*dXsFunCol + 2*dYsFunRow*dYsFunCol);
                    Kxy_Entry += (massOfP_*mu_*dt_/rho_)*(dYsFunRow*dXsFunCol);
                }
            }
        }
        unsigned position_ = num_ + n_;
        // std::cout << "position: "<<  position << std::endl;
        theMatrix->setNonZeroKEntries(Kxx_Entry, Kyy_Entry, Kxy_Entry, (position_));
         ++num_;
         // std::cout << "MEntry: " << MEntry << std::endl;
         // std::cout << std::endl;
    }
    n_ += nodesOfNode_.size();
    // std::cout << "n: " << n << std::endl;

}

//------------------------------------------------------------------------------
                    // evaluate the entries of known F matrix

template<typename MPM_ITEMS>
void mpm::SolverMpm<MPM_ITEMS>::computeNonZeroF_Entries(NodePtr row, double dt, unsigned position) {
    VecDof velocityTerm = (row->giveLumpedMassAtNodes())*(row->giveSoilVelocity());
    VecDof pressureTerm = row->givePressureTermAtNodes();
    VecDof bodyForceTerm = row->giveExtForceTermAtNodes();

    VecDof F_Entry = velocityTerm + dt*pressureTerm + dt*bodyForceTerm;
    // std::cout << "F_Entry: " << F_Entry << std::endl;

    theMatrix->setNonZeroFEntries(F_Entry,position);
}

//------------------------------------------------------------------------------
