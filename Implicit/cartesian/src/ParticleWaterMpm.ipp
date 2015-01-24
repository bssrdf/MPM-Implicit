// *********************************************************************************
// File: ParticleWaterMpm.ipp
//
// Details: Implementation of base class for the MPM Water Particles
//
// Dependency: MPM_ITEMS contains all template parameters
//
// Author: Samila Bandara, University of Cambridge
//
// Version: 1.0
// *********************************************************************************


// Write only the ID of this node to a given output stream
// \param[in,out]  out Output stream
// \retval         out Output stream


// cache material to particle and initialise material properties of particle
// template<typename MPM_ITEMS>
// void mpm::ParticleSoilMpm<MPM_ITEMS>::cacheMaterial(MaterialPtr m) {
//     material_ = m;
//     // density_ = material_->density();
//     // this -> computeMass_();
//     // porosity_ = material_->porosity();
//     return;
// }


// calculate mass of particle (constant across the calculation)
// when there is no soil particles around porosity_ should be zero
template<typename MPM_ITEMS>
void mpm::ParticleWaterMpm<MPM_ITEMS>::computeMass_() {
    double pSpacing;
    pSpacing = this->ParticleType::givePartSpacing_();
    if (dim == 2)
        mass_ = pSpacing * pSpacing * 1. * density_ * (porosity_);
    else if (dim == 3)
        mass_ = pSpacing * pSpacing * pSpacing * density_ * (porosity_);
    return;
}

