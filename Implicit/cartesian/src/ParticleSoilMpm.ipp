// *********************************************************************************
// File: ParticleSoilMpm.ipp
//
// Details: Implementation of base class for the MPM Soil Particles
//
// Dependency: MPM_ITEMS contains all template parameters
//
// Author: Krishna Kumar, Samila Bandara, University of Cambridge
//
// Version: 1.0
// *********************************************************************************

// cache material to particle and initialise material properties of particle
template<typename MPM_ITEMS>
void mpm::ParticleSoilMpm<MPM_ITEMS>::cacheMaterial(std::vector<MaterialPtr> materialPtrs, const std::vector<double> pSpacing) {
    material_ = materialPtrs.at(matID_);
    double grainDensity = material_->density();
    porosity_ = material_ -> porosity();
    density_ = grainDensity * (1. - porosity_);
    phi_ = material_ -> phi();
    this->computeMass_(pSpacing);
    return;
}


// calculate mass of particle (constant across the calculation)
template<typename MPM_ITEMS>
void mpm::ParticleSoilMpm<MPM_ITEMS>::computeMass_(const std::vector<double> pSpacing) {
    mass_ = pSpacing.at(0) * pSpacing.at(1) * pSpacing.at(2) * density_;
    return;
}


// Assign mass to relevant nodes
template<typename MPM_ITEMS>
void mpm::ParticleSoilMpm<MPM_ITEMS>::assignMassToNodes() {
    VecNN nMass = mass_ * sfunc_;
    for (unsigned i = 0; i < numNodes; i ++)
        nodesOfParticle_[i]->assignSoilMass(nMass(i));
    return;
}


// Assign velocity to relevant nodes
template<typename MPM_ITEMS>
void mpm::ParticleSoilMpm<MPM_ITEMS>::assignMomentumToNodes() {
    VecDof nMomentum;
    for (unsigned i = 0; i < numNodes; i ++) {
        nMomentum = mass_ * velocity_ * sfunc_(i);
        nodesOfParticle_[i]->assignSoilMomentum(nMomentum);
    }
    return;
}


// Compute strain vector of soil particle { \f$ \epsilon_{xx}, \epsilon_{yy},
// \epsilon_{zz}, \gamma_{xy}, \gamma_{yz}, \gamma_{xz} \f$ } using particle
// velocity and B matrix. (Note:\f$ \gamma_{xy} = 2\epsilon_{xy} \f$)

template<typename MPM_ITEMS>
void mpm::ParticleSoilMpm<MPM_ITEMS>::computeStrain(const double dt) {
    VecNComp strainRate, tempSR;
    VecNNVec nVelocity;

    strainRate.clear();
    dStrain_.clear();

    for (unsigned i = 0; i < numNodes; i ++)
        nVelocity(i)  = nodesOfParticle_[i]->giveSoilVelocity();

    this->computeBMatrix_();
    for (unsigned i = 0; i < numNodes; i ++) {
        ublas::noalias(tempSR) = ublas::prod(B_(i), nVelocity(i));
        strainRate += tempSR;             // tempSR is used for efficiency
    }
    if (dim == 2) {
        dStrain_(0) = strainRate(0) * dt;
        dStrain_(1) = strainRate(1) * dt;
        dStrain_(3) = strainRate(2) * dt;
    } else if (dim == 3)
        dStrain_ = strainRate * dt;

    strain_ += dStrain_;

    // Alternative approach
    // MatDimDim velocityGradient, strainRate;
    // VecDof nodeVelocity;
    // MatNNDim nAllVelocity;

    // for (unsigned i = 0; i < numNodes; i ++) {
    //   nodeVelocity = nodesOfParticle_[i]->giveSoilVelocity();
    //   for (unsigned j = 0; j < dim; j ++)
    // 	nAllVelocity(i, j) = nodeVelocity(j);
    // }

    // velocityGradient = 0.5 * (ublas::prod(dsfuncDX_, nAllVelocity));
    // strainRate = velocityGradient + ublas::trans(velocityGradient);

    return;

}


// Compute strain vector of soil particle { \f$ \epsilon_{xx}, \epsilon_{yy},
// \epsilon_{zz}, \gamma_{xy}, \gamma_{yz}, \gamma_{xz} \f$ } using particle
// velocity and B matrix. (Note:\f$ \gamma_{xy} = 2\epsilon_{xy} \f$)

template<typename MPM_ITEMS>
void mpm::ParticleSoilMpm<MPM_ITEMS>::computeStrainBBar(const double dt) {
    VecNComp strainRateBBar, tempSRBBar;
    VecNNVec nVelocityBBar;

    strainRateBBar.clear();
    dStrainBBar_.clear();

    for (unsigned i = 0; i < numNodes; i ++)
        nVelocityBBar(i)  = nodesOfParticle_[i]->giveSoilVelocity();

    this->computeBBarMatrix_();

    for (unsigned i = 0; i < numNodes; i ++) {
        ublas::noalias(tempSRBBar) = ublas::prod(BBar_(i), nVelocityBBar(i));
        strainRateBBar += tempSRBBar;             // tempSR is used for efficiency
    }

    if (dim == 2) {
        dStrainBBar_(0) = strainRateBBar(0) * dt;
        dStrainBBar_(1) = strainRateBBar(1) * dt;
        dStrainBBar_(3) = strainRateBBar(2) * dt;
        pressure_ = pressure_ + 1.67E6*(dStrainBBar_(0) + dStrainBBar_(1));
    }
    else if (dim == 3) {
        dStrainBBar_ = strainRateBBar * dt;
        pressure_ = pressure_ + 1.67E6*(dStrainBBar_(0) + dStrainBBar_(1) + dStrainBBar_(2));
    }

    strain_ += dStrainBBar_;

    return;

}


// Compute Volumetric Strain Increament at the center of the element
template<typename MPM_ITEMS>
void mpm::ParticleSoilMpm<MPM_ITEMS>::computeVolStrainIncCenter(const double dt) {
    VecNComp strainRateCenter, tempSRCenter;
    VecNNVec nVelocityCenter;

    strainRateCenter.clear();
    dStrainCenter_.clear();

    for (unsigned i = 0; i < numNodes; i ++)
        nVelocityCenter(i)  = nodesOfParticle_[i]->giveSoilVelocity();

    this->computeBMatrixCenter_();

    for (unsigned i = 0; i < numNodes; i ++) {
        ublas::noalias(tempSRCenter) = ublas::prod(BCenter_(i), nVelocityCenter(i));
        strainRateCenter += tempSRCenter;             // tempSR is used for efficiency
    }

    if (dim == 2) {
        dStrainCenter_(0) = strainRateCenter(0) * dt;
        dStrainCenter_(1) = strainRateCenter(1) * dt;
        dStrainCenter_(2) = strainRateCenter(2) * dt;
        dStrainCenter_(3) = strainRateCenter(3) * dt;
    }
    else if (dim == 3) {
        dStrainCenter_ = strainRateCenter * dt;
    }

}


// Compute stress vector of soil particle { \f$ \sigma_{xx}, \sigma_{yy},
// \sigma_{zz}, \sigma_{xy}, \sigma_{yz}, \sigma_{xz} \f$ } using a soil model.
// (Note:\f$ \sigma_{xy} = 2\tau_{xy} \f$)

template<typename MPM_ITEMS>
void mpm::ParticleSoilMpm<MPM_ITEMS>::computeStress(ParticleCloudSoilPtr pCloudSoilPtr) {
    MatDimDim velocityGradient;
    VecDof nodeVelocity;
    MatNNDim nAllVelocity;

    for (unsigned i = 0; i < numNodes; i ++) {
        nodeVelocity = nodesOfParticle_[i]->giveSoilVelocity();
        for (unsigned j = 0; j < dim; j ++)
            nAllVelocity(i, j) = nodeVelocity(j);
    }
    RateOfStrainI2_ = 0.;
    velocityGradient = 0.5 * (ublas::prod(dsfuncDX_, nAllVelocity));
    strainRate_ = velocityGradient + ublas::trans(velocityGradient);

    if (dim == 2)
        material_->computeStress(dStrainBBar_, pressure_, stress_ , pCloudSoilPtr, id_);
    else
        material_->computeStress3D(dStrainBBar_, pressure_, stress_ , pCloudSoilPtr, id_);

    return;
}


// Compute strain matrix (B matrix).

template<typename MPM_ITEMS>
void mpm::ParticleSoilMpm<MPM_ITEMS>::computeBMatrix_() {

    MatNCompDim Bn;
    B_.clear();

    if (dim == 2) {
        for (unsigned i = 0; i < numNodes; i ++) {
            Bn(0,0)= dsfuncDX_(0,i);
            Bn(0,1)= 0;
            Bn(1,0)= 0;
            Bn(1,1)= dsfuncDX_(1,i);
            Bn(2,0)= dsfuncDX_(1,i);
            Bn(2,1)= dsfuncDX_(0,i);
            B_(i) = Bn;
        }
    }
    else if (dim ==3) {
        for (unsigned i = 0; i < numNodes; i ++) {
            Bn(0,0)= dsfuncDX_(0,i);
            Bn(0,1)= 0;
            Bn(0,2)= 0;
            Bn(1,0)= 0;
            Bn(1,1)= dsfuncDX_(1,i);
            Bn(1,2)= 0;
            Bn(2,0)= 0;
            Bn(2,1)= 0;
            Bn(2,2)= dsfuncDX_(2,i);
            Bn(3,0)= dsfuncDX_(1,i);
            Bn(3,1)= dsfuncDX_(0,i);
            Bn(3,2)= 0;
            Bn(4,0)= 0;
            Bn(4,1)= dsfuncDX_(2,i);
            Bn(4,2)= dsfuncDX_(1,i);
            Bn(5,0)= dsfuncDX_(2,i);
            Bn(5,1)= 0;
            Bn(5,2)= dsfuncDX_(0,i);
            B_(i) = Bn;
        }
    }

    return;

}


// Compute strain matrix (BBar matrix). Using \={B}

template<typename MPM_ITEMS>
void mpm::ParticleSoilMpm<MPM_ITEMS>::computeBBarMatrix_() {

    MatNCompDim BBarN;
    BBar_.clear();
    double BBar00, BBar01, BBar10, BBar11;
    BBar00 = BBar01 = BBar10 = BBar11 = 0.;
    const double OneThird = 1.0/3.0;
    if (dim == 2) {
        for (unsigned i = 0; i < numNodes; i ++) {
            BBar00 = OneThird * (2*dsfuncDX_(0,i)+dsfuncDXCenter_(0,i));
            BBar01 = OneThird * (dsfuncDXCenter_(1,i) - dsfuncDX_(1,i));
            BBar10 = OneThird * (dsfuncDXCenter_(0,i) - dsfuncDX_(0,i));
            BBar11 = OneThird * (2*dsfuncDX_(1,i)+dsfuncDXCenter_(1,i));

            BBarN(0,0)= BBar00;
            BBarN(0,1)= BBar01;
            BBarN(1,0)= BBar10;
            BBarN(1,1)= BBar11;


            BBarN(2,0)= dsfuncDX_(1,i);
            BBarN(2,1)= dsfuncDX_(0,i);
            BBar_(i) = BBarN;
        }
    }
    else if (dim ==3) {
        for (unsigned i = 0; i < numNodes; i ++) {
            BBarN(0,0)= OneThird*(2*dsfuncDX_(0,i)+dsfuncDXCenter_(0,i));
            BBarN(0,1)= OneThird*(dsfuncDXCenter_(1,i) - dsfuncDX_(1,i));
            BBarN(0,2)= OneThird*(dsfuncDXCenter_(2,i) - dsfuncDX_(2,i));
            BBarN(1,0)= OneThird*(dsfuncDXCenter_(0,i) - dsfuncDX_(0,i));
            BBarN(1,1)= OneThird*(2* dsfuncDX_(1,i)+dsfuncDXCenter_(1,i));
            BBarN(1,2)= OneThird*(dsfuncDXCenter_(2,i) - dsfuncDX_(2,i));
            BBarN(2,0)= OneThird*(dsfuncDXCenter_(0,i) - dsfuncDX_(0,i));
            BBarN(2,1)= OneThird*(dsfuncDXCenter_(1,i) - dsfuncDX_(1,i));
            BBarN(2,2)= OneThird*(2*dsfuncDX_(2,i)+dsfuncDXCenter_(2,i));


            BBarN(3,0)= dsfuncDX_(1,i);
            BBarN(3,1)= dsfuncDX_(0,i);
            BBarN(3,2)= 0.;
            BBarN(4,0)= 0.;
            BBarN(4,1)= dsfuncDX_(2,i);
            BBarN(4,2)= dsfuncDX_(1,i);
            BBarN(5,0)= dsfuncDX_(2,i);
            BBarN(5,1)= 0.;
            BBarN(5,2)= dsfuncDX_(0,i);
            BBar_(i) = BBarN;
        }
    }

    return;
}


// Compute strain matrix (B matrix) at the element center.

template<typename MPM_ITEMS>
void mpm::ParticleSoilMpm<MPM_ITEMS>::computeBMatrixCenter_() {

    MatNCompDim BnCenter;
    BCenter_.clear();

    for (unsigned i = 0; i < numNodes; i ++) {
        if (dim == 2) {
            BnCenter(0,0)= dsfuncDXCenter_(0,i);
            BnCenter(0,1)= 0;
            BnCenter(1,0)= 0;
            BnCenter(1,1)= dsfuncDXCenter_(1,i);
            BnCenter(2,0)= dsfuncDXCenter_(1,i);
            BnCenter(2,1)= dsfuncDXCenter_(0,i);
        }
        else if (dim ==3) {
            BnCenter(0,0)= dsfuncDXCenter_(0,i);
            BnCenter(0,1)= 0;
            BnCenter(0,2)= 0;
            BnCenter(1,0)= 0;
            BnCenter(1,1)= dsfuncDXCenter_(1,i);
            BnCenter(1,2)= 0;
            BnCenter(2,0)= 0;
            BnCenter(2,1)= 0;
            BnCenter(2,2)= dsfuncDXCenter_(2,i);
            BnCenter(3,0)= dsfuncDXCenter_(1,i);
            BnCenter(3,1)= dsfuncDXCenter_(0,i);
            BnCenter(3,2)= 0;
            BnCenter(4,0)= 0;
            BnCenter(4,1)= dsfuncDXCenter_(2,i);
            BnCenter(4,2)= dsfuncDXCenter_(1,i);
            BnCenter(5,0)= dsfuncDXCenter_(2,i);
            BnCenter(5,1)= 0;
            BnCenter(5,2)= dsfuncDXCenter_(0,i);
        }
        BCenter_(i) = BnCenter;
    }

    return;

}


// Store Traction on Particles
template<typename MPM_ITEMS>
void mpm::ParticleSoilMpm<MPM_ITEMS>::storeTraction(unsigned direction, double tractionPressure) {
    tractionBn_.push_back(std::make_pair(direction, tractionPressure));
    return;
}

//------------------------------------------------------------------------------
                  // Assign velocity to relevant nodes

template<typename MPM_ITEMS>
void mpm::ParticleSoilMpm<MPM_ITEMS>::assignBodyForceToNodes(const VecDim G) {
    VecDof nBodyForce;

    for (unsigned i = 0; i < numNodes; i ++) {
        nBodyForce = mass_ * G * sfunc_(i);
        nodesOfParticle_[i]->assignExternalForce(nBodyForce);
        // std::cout << "nBodyForce : " << nBodyForce << std::endl;
    }
    return;
}

//------------------------------------------------------------------------------
                  // Assign velocity to relevant nodes

template<typename MPM_ITEMS>
void mpm::ParticleSoilMpm<MPM_ITEMS>::assignPressureToNodes() {
    VecDof nPressure;

    for (unsigned i = 0; i < numNodes; i++) {
        for (unsigned j = 0; j < dim; j++) {
            nPressure(j) = (mass_/density_)*dsfuncDX_(j,i)*pressure_;
        }
        // std::cout << "nPressure: " << nPressure << std::endl;
        nodesOfParticle_[i]->assignPressureForce(nPressure);
    }
}

//------------------------------------------------------------------------------

// Compute internal forces for relevant nodes using stress
template<typename MPM_ITEMS>
void mpm::ParticleSoilMpm<MPM_ITEMS>::computeInternalForceAtNodes() {
    VecDof nIntForce;

    if (dim == 2) {
        VecNComp tmpStress;
        tmpStress(0) = stress_(0);
        tmpStress(1) = stress_(1);
        tmpStress(2) = stress_(3);

        for (unsigned i = 0; i < numNodes; i ++) {
            nIntForce = (mass_ / density_) * (ublas::prod(ublas::trans(B_(i)), tmpStress));
            nodesOfParticle_[i]->assignInternalForce(nIntForce);
        }
    }
    else if (dim == 3) {
        for (unsigned i = 0; i < numNodes; i ++) {
            nIntForce = (mass_ / density_) * (ublas::prod(ublas::trans(B_(i)), stress_));
            nodesOfParticle_[i]->assignInternalForce(nIntForce);
        }
    }
    return;
}


// Assign tractions to nodes
template<typename MPM_ITEMS>
void mpm::ParticleSoilMpm<MPM_ITEMS>::assignTractionsToNodes(const std::vector<double> pSpacing) {
    VecDof nTractionF;
    nTractionF.clear();

    TractionVec_::const_iterator cIter = tractionBn_.begin();
    TractionVec_::const_iterator cEnd  = tractionBn_.end();

    for (; cIter != cEnd; cIter ++) {
        for (unsigned i = 0; i < numNodes; i ++) {
            nTractionF(cIter->first) = (mass_ / density_) *
                                       ((cIter->second) / pSpacing.at(0)) * sfunc_(i);
            nodesOfParticle_[i]->assignExternalForce(nTractionF);
        }
    }

    return;
}


// Update velocity and position
template<typename MPM_ITEMS>
void mpm::ParticleSoilMpm<MPM_ITEMS>::updateSoilParticle(const double dt) {
    VecDof nVelocity, nAcceleration, tempVal, tempVelocity, tempCoord;

    // Setting temporary velocity to zero
    ublas::noalias(tempVelocity) = 0.* nodesOfParticle_[0]->giveSoilVelocity();
    ublas::noalias(tempCoord)    = 0.* nodesOfParticle_[0]->giveSoilVelocity();

    // Both Velocity and Coord are multiplied by dt at the end of this function
    for (unsigned i = 0; i < numNodes; i ++) {
        nVelocity  = nodesOfParticle_[i]->giveSoilVelocity();
        ublas::noalias(tempVal) = nVelocity * sfunc_(i);
        ublas::noalias(tempCoord) += tempVal;

        //      nAcceleration  = nodesOfParticle_[i]->giveSoilAcceleration();
        //      ublas::noalias(tempVal) = nAcceleration * sfunc_(i);
        //      ublas::noalias(tempVelocity) += tempVal;
    }


    // Setting a tolerance value especially for y-direction
    // Due to the sequence of summation, there is always a small error in y
    for (unsigned i = 0; i < tempVelocity.size(); i++) {
        if (std::fabs(tempCoord[i]) < 1.0e-16) {
            tempCoord[i] = 0.;
        }
        if (std::fabs(tempVelocity[i]) < 1.0e-16) {
            tempVelocity[i] = 0.;
        }

    }

//    velocity_ += tempVelocity * dt;
    velocity_ += tempCoord;
    coord_ += tempCoord * dt;

    return;
}


// Update phase averaged density of soil (\f$ =(1-n)*\rho_{grain}\f$)
// using mass balance equation  \f$ \frac{d\rho}{dt}+\rho\nabla\cdot{v_s}=0 \f$

template<typename MPM_ITEMS>
void mpm::ParticleSoilMpm<MPM_ITEMS>::updateSoilDensity() {
    double volStrain = 0;

    for (unsigned i = 0; i < dim; i ++)
        volStrain += dStrainBBar_(i);

    density_ = density_ / (1. + volStrain);

    return;
}


// Write particle data for VTK file (continuation from particleCloud's
// writeParticleCloudData function)
// \param[out] out  output file

template<typename MPM_ITEMS>
void mpm::ParticleSoilMpm<MPM_ITEMS>::writeParticleData(std::ostream& out) {
    if (dim == 2)
        out << coord_(0) << " " << coord_(1) << " 0." << std::endl;
    else if (dim == 3)
        out << coord_(0) << " " << coord_(1) << " "<< coord_(2) << std::endl;
    return;
}


// Write particle data for VTK file (continuation from particleCloud's
// writeParticleCloudData function)
// \param[out] out  output file

template<typename MPM_ITEMS>
void mpm::ParticleSoilMpm<MPM_ITEMS>::writeParticleVelocityVTK(std::ostream& out) {
    if (dim == 2)
        out << velocity_(0) << " " << velocity_(1) << " 0." << std::endl;
    else if (dim == 3)
        out << velocity_(0) << " " << velocity_(1) << " " << velocity_(2) << std::endl;
    return;
}


template<typename MPM_ITEMS>
void mpm::ParticleSoilMpm<MPM_ITEMS>::writeParticleStressVTK(std::ostream& out) {
    out << stress_(0) << " " << stress_(1) << " " << stress_(2) << std::endl;
    // out << stress_(0) << " " << stress_(3) << " " << stress_(5) << std::endl
    //     << stress_(3) << " " << stress_(1) << " " << stress_(4) << std::endl
    //     << stress_(5) << " " << stress_(4) << " " << stress_(2) << std::endl << std::endl;
    return;
}


template<typename MPM_ITEMS>
void mpm::ParticleSoilMpm<MPM_ITEMS>::writeParticleStrainVTK(std::ostream& out) {
    if (dim == 2)
        out << strain_(0) << " " << strain_(1) << " 0." << std::endl;
    else if (dim == 3)
        out << strain_(0) << " " << strain_(1) << " " << strain_(2) << std::endl;
    // out << strain_(0) << " " << strain_(3) << " " << strain_(5) << std::endl
    //     << strain_(3) << " " << strain_(1) << " " << strain_(4) << std::endl
    //     << strain_(5) << " " << strain_(4) << " " << strain_(2) << std::endl << std::endl;
    return;
}


template<typename MPM_ITEMS>
void mpm::ParticleSoilMpm<MPM_ITEMS>::writeParticleRateOfStrainI2VTK(std::ostream& out) {
    out  << pressure_ << std::endl;
    return;
}


// store the constraints
//
// template<typename MPM_ITEMS>
// void mpm::ParticleMpm<MPM_ITEMS>::computeStiffness(const VecDim& xi) {
//     VecNV phi; this->sfun_ (xi, phi);
//     //std::cout << "sfun:phi = "<< phi[1] << std::endl;
//     return;
// }


// void setInitialStress(const double hOverburden, const VecDim G);
// TAKES VERY LONGE TIME
// calculate initial stress of particle using k0
//
// template<typename MPM_ITEMS>
// void mpm::ParticleSoilMpm<MPM_ITEMS>::setInitialStress(const double hOverburden, const VecDim G) {
//     double PI = std::atan(1.0) * 4;
//
//     if (dim == 2) {
//       stress_(1) = density_ * G(1) * hOverburden;
//       stress_(0) = (1 - sin(phi_ * PI / 180)) * stress_(1);
//       stress_(2) = stress_(0);
//       std::cout<<id_<<" "<<stress_(0)<<" "<<stress_(1)<<" "<<stress_(2)<<" 0 0 0"<<"\n";
//     }
//     else if (dim == 3) {
//       stress_(2) = density_ * G(2) * hOverburden;
//       stress_(0) = (1 - sin(phi_ * PI / 180)) * stress_(2);
//       stress_(1) = stress_(0);
//     }
//     return;
// }
