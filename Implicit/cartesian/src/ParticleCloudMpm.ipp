// *********************************************************************************
// File: ParticleCloudMpm.ipp
//
// Details: Implementation of container class for the MPM particles
//
// Dependency: MPM_ITEMS contains all template parameters
//
// Author: Krishna Kumar and Samila Bandara, University of Cambridge
//
// Version: 2.0 revised by Krishna
// *********************************************************************************

#include <limits>
#include <boost/tuple/tuple.hpp>

#include "cartesian/src/MpmItems.hpp"


// Constructor given a stream to a simple particle file. This file must have the
// following format:
//    NP dX \n
//    x1  y1  z1  \n
//    ....        \n

// where NP is the number of particles, dX is the size of the particle (spacing).
// x, y, z are the coordinates of each particle
// \param smf Input stream

template<typename MPM_ITEMS, typename PARTICLE_TYPE>
mpm::ParticleCloudMpm<MPM_ITEMS, PARTICLE_TYPE>::ParticleCloudMpm(std::istream& prcld) {
    unsigned numParticles;
    // read number of particles
    prcld >> numParticles;
    prcld.ignore(std::numeric_limits<std::streamsize>::max(), '\n');

    // read particle spacing and
    particles_.reserve(numParticles);
    for (unsigned n = 0; n < numParticles; n++) {
        ParticlePtr_ theParticle = new ParticleType(n);
        theParticle->readParticleSelf(prcld);
        particles_.push_back(theParticle);
    }
    return;
}

// Free all dynamically allocated memory
template<typename MPM_ITEMS, typename PARTICLE_TYPE>
mpm::ParticleCloudMpm<MPM_ITEMS, PARTICLE_TYPE>::~ParticleCloudMpm() {
    this->iterateOverParticles(corlib::deleteFunctor());
    particles_.clear();
    return;
}


// Read from a stream of Material Types. The file must have the following format
// NP  NM mID
// pid_1 mID_1
// pid_2 mID_1
// pid_3 mID_1
// pid_4 mID_2
// ...
//
// or if there is only one material type for all the particles you can use the
// following format
//
// NP NM mID
//
// where NP is number of particles, pid is the particle id, mID is the material ID
// as defined in the Material file.

template<typename MPM_ITEMS, typename PARTICLE_TYPE>
void mpm::ParticleCloudMpm<MPM_ITEMS, PARTICLE_TYPE>::associateParticleToMaterial(std::istream& matParticleReader) {
    unsigned MaxNumberOfParticles, particle;
    unsigned MaterialTypes, material;

    matParticleReader >> MaxNumberOfParticles >> MaterialTypes;
    matParticleReader >> material;

    for (unsigned i = 0; i < MaxNumberOfParticles; i++) {
        if (MaterialTypes > 1)
            matParticleReader >> particle >> material;
        else
            particle = i;
        particles_.at(particle)->setMaterial(material);
    }
}


// Read from a stream the initial stresses. The file must have the format

//  NP
//  pid_1    s_1   s_2   s_3   s_4   s_5   s_6   \n
//  ...                                          \n
//  pid_N1   s_N1  s_N2  s_N3  s_N4  s_N5  s_N6  \n
//
//  where NP is the number of particles in this file, pid is the particle id,
//  s_1 to s_6 are the initial stresses. Each particle should contain 6 stress
//  components.
//  \param[in] initStress Input stream containing the initial stresses

template<typename MPM_ITEMS, typename PARTICLE_TYPE>
void mpm::ParticleCloudMpm<MPM_ITEMS, PARTICLE_TYPE>::setInitialStressParticles(std::istream& initStress) {
    unsigned particleNum;
    Vec6x1 iStress;

    // get number of particles with initial stress
    unsigned numParticleInStrs = 0;
    initStress >> numParticleInStrs;

    // go through initial stresses
    for (unsigned c = 0; c < numParticleInStrs; c++) {
        initStress >> particleNum;
        for (unsigned i = 0; i < 6; i++)
            initStress >> iStress(i);
        particles_.at(particleNum)->setInitialStress(iStress);
    }
    return;
}


// Read from a stream the traction boundary. The file must have the format
//
// NTP            \n
// p_1     d_1    tp_1    \n
// ...                    \n
// p_NTP   d_NTP  tp_NTP  \n
//
// where NTP is the number of traction particles in this file, which are used
// to give traction forces in the d_i direction. p_i is the particle id,
// d_i is the direction number (0|1|2) of traction force. tp_i is the traction
// pressure. NTP is the number of particles with traction pressures.
// \param[in] trbn Input stream containing the traction pressures

template<typename MPM_ITEMS, typename PARTICLE_TYPE>
void mpm::ParticleCloudMpm<MPM_ITEMS, PARTICLE_TYPE>::readTractionBoundary(std::istream& trbn) {
    // get number of traction particles
    unsigned numTractionP;
    trbn >> numTractionP;
    // go through particles
    for (unsigned c = 0; c < numTractionP; c++) {
        unsigned pNum, direction;
        double tractionPressure;
        trbn >> pNum >> direction >> tractionPressure;
        particles_.at(pNum)->storeTraction(direction, tractionPressure);
    }

    return;
}


// Apply a given functor to all particles making use of std::for_each
// \tparam OP         Type of the particle operation
// \param[in,out] op  Specific operation applied to all particles
// \retval OP         for_each returns the operator

template<typename MPM_ITEMS, typename PARTICLE_TYPE>
template<typename OP>
OP mpm::ParticleCloudMpm<MPM_ITEMS, PARTICLE_TYPE>::iterateOverParticles(OP op) {
    typename ParticleVec::iterator begin = particles_.begin();
    typename ParticleVec::iterator end   = particles_.end();
    return std::for_each(begin, end, op);
}

// Parallel Version - GNU/Parallel
template<typename MPM_ITEMS, typename PARTICLE_TYPE>
template<typename OP>
OP mpm::ParticleCloudMpm<MPM_ITEMS, PARTICLE_TYPE>::iterateOverParticlesParallel(OP op) {
    typename ParticleVec::iterator begin = particles_.begin();
    typename ParticleVec::iterator end   = particles_.end();
    return __gnu_parallel::for_each(begin, end, op);
}


// Apply a given functor only to those particles which have a given predicate.
// Given the predicate 'p', the functor 'op' will only be applied to those
// particles for which holds 'p(X) = true'.
// \tparam OP          Type of functor to be applied to particles
// \tparam PRED        Type of predicate functor
// \param[in,out] op   The functor to be applied
// \param[in] p        The predicate
// \retval op          Returns the given functor

template<typename MPM_ITEMS, typename PARTICLE_TYPE>
template<typename OP, typename PRED>
OP mpm::ParticleCloudMpm<MPM_ITEMS, PARTICLE_TYPE>::iterateOverParticlesWithPredicate(OP op, PRED p) {
    typename ParticleVec::iterator begin = particles_.begin();
    typename ParticleVec::iterator end   = particles_.end();
    return corlib::for_each_if (begin, end, op, p);
}


// Given an number numParticle, return the correponding particle pointer
// \param[in] numParticle  Number of the requested particle
// \retval    ParticlePtr_ Pointer to the requested particle

template<typename MPM_ITEMS, typename PARTICLE_TYPE>
PARTICLE_TYPE* mpm::ParticleCloudMpm<MPM_ITEMS, PARTICLE_TYPE>::
particleIterator(const unsigned numParticle) const {
    return particles_.at(numParticle);
}


// Write particle data for VTK file
// \param[out] out  output file

template<typename MPM_ITEMS, typename PARTICLE_TYPE>
void mpm::ParticleCloudMpm<MPM_ITEMS, PARTICLE_TYPE>::writeParticleCloudVelocityData(const unsigned int step, std::ostream& out) {

    unsigned numParticles = particles_.size();
    out << "# vtk DataFile Version 2.0" << std::endl;
    out << "test grain" << std::endl << "ASCII" << std::endl;
    out << "DATASET UNSTRUCTURED_GRID" << std::endl;

    out << "POINTS " << numParticles << " float" << std::endl;
    for (unsigned i = 0; i < numParticles; i++)
        particles_[i]->writeParticleData(out);

    out << "CELLS " << numParticles << " " << 2 * numParticles << std::endl;
    for (unsigned i = 0; i < numParticles; i++)
        out << "1 " << i << std::endl;

    out << "CELL_TYPES " << numParticles << std::endl;
    for (unsigned i = 0; i < numParticles; i++)
        out << "1 " << std::endl;

    out << "POINT_DATA " << numParticles << std::endl;
    out << "VECTORS Velocity float" << std::endl;
    for (unsigned i = 0; i < numParticles; i++)
        particles_[i]->writeParticleVelocityVTK(out);

    return;
}

template<typename MPM_ITEMS, typename PARTICLE_TYPE>
void mpm::ParticleCloudMpm<MPM_ITEMS, PARTICLE_TYPE>::writeParticleCloudStressData(const unsigned int step, std::ostream& out) {

    unsigned numParticles = particles_.size();
    out << "# vtk DataFile Version 2.0" << std::endl;
    out << "test grain" << std::endl << "ASCII" << std::endl;
    out << "DATASET UNSTRUCTURED_GRID" << std::endl;

    out << "POINTS " << numParticles << " float" << std::endl;
    for (unsigned i = 0; i < numParticles; i++)
        particles_[i]->writeParticleData(out);

    out << "CELLS " << numParticles << " " << 2 * numParticles << std::endl;
    for (unsigned i = 0; i < numParticles; i++)
        out << "1 " << i << std::endl;

    out << "CELL_TYPES " << numParticles << std::endl;
    for (unsigned i = 0; i < numParticles; i++)
        out << "1 " << std::endl;

    out << "POINT_DATA " << numParticles << std::endl;
    out << "VECTORS Stress float" << std::endl;
    for (unsigned i = 0; i < numParticles; i++)
        particles_[i]->writeParticleStressVTK(out);

    return;
}

template<typename MPM_ITEMS, typename PARTICLE_TYPE>
void mpm::ParticleCloudMpm<MPM_ITEMS, PARTICLE_TYPE>::writeParticleCloudStrainData(const unsigned int step, std::ostream& out) {

    unsigned numParticles = particles_.size();
    out << "# vtk DataFile Version 2.0" << std::endl;
    out << "test grain" << std::endl << "ASCII" << std::endl;
    out << "DATASET UNSTRUCTURED_GRID" << std::endl;

    out << "POINTS " << numParticles << " float" << std::endl;
    for (unsigned i = 0; i < numParticles; i++)
        particles_[i]->writeParticleData(out);

    out << "CELLS " << numParticles << " " << 2 * numParticles << std::endl;
    for (unsigned i = 0; i < numParticles; i++)
        out << "1 " << i << std::endl;

    out << "CELL_TYPES " << numParticles << std::endl;
    for (unsigned i = 0; i < numParticles; i++)
        out << "1 " << std::endl;

    out << "POINT_DATA " << numParticles << std::endl;
    out << "VECTORS Strain float" << std::endl;
    for (unsigned i = 0; i < numParticles; i++)
        particles_[i]->writeParticleStrainVTK(out);

    return;
}


template<typename MPM_ITEMS, typename PARTICLE_TYPE>
void mpm::ParticleCloudMpm<MPM_ITEMS, PARTICLE_TYPE>::writeParticleCloudRateOfStrainI2Data(const unsigned int step, std::ostream& out) {

    unsigned numParticles = particles_.size();
    out << "# vtk DataFile Version 2.0" << std::endl;
    out << "test grain" << std::endl << "ASCII" << std::endl;
    out << "DATASET UNSTRUCTURED_GRID" << std::endl;

    out << "POINTS " << numParticles << " float" << std::endl;
    for (unsigned i = 0; i < numParticles; i++)
        particles_[i]->writeParticleData(out);

    out << "CELLS " << numParticles << " " << 2 * numParticles << std::endl;
    for (unsigned i = 0; i < numParticles; i++)
        out << "1 " << i << std::endl;

    out << "CELL_TYPES " << numParticles << std::endl;
    for (unsigned i = 0; i < numParticles; i++)
        out << "1 " << std::endl;

    out << "POINT_DATA " << numParticles << std::endl;
    out << "SCALARS RateOfStrainI2 float" << std::endl;
    out << "LOOKUP_TABLE default" << std::endl;
    for (unsigned i = 0; i < numParticles; i++)
        particles_[i]->writeParticleRateOfStrainI2VTK(out);
    return;
}


// Read from a stream the nodal constraints. The file must have the format
//
//   NC
//   n_1    d_1    v_1  \n
//   ...                \n
//  n_NC   d_NC   v_NC \n
//
// where NC is the number of constraints in this file, n_i is the node number,
// d_i is the direction number (0|1|2), and v_i is the value of the i-th
// constraint.
// \param[in] cstr Input stream containing the constraints

template< typename MPM_ITEMS, typename PARTICLE_TYPE>
void mpm::ParticleCloudMpm< MPM_ITEMS, PARTICLE_TYPE>::readParticleForces(std::istream& pfrcs) {
    // get number of constraints
    unsigned numParticleForces = 0;
    pfrcs >> numParticleForces;
    // go through contraints
    for (unsigned c = 0; c < numParticleForces; c++) {
        unsigned particleNum, component;
        double value;
        pfrcs >> particleNum >> component >> value;
        particles_.at(particleNum)->storeParticleForce(component, value);
    }
}


// Set the particle constraints given iterators to a container of constraint
// triples: <n,d,v> with the particle number \e n, the constraint direction \e d,
// and the value of the constraint \e v.
// \tparam IT       Type of iterator to triple container
// \param[in] iter  Begin of the triple container
// \param[in] end   End   of the triple container

template<typename MPM_ITEMS, typename PARTICLE_TYPE>
template<typename IT>
void mpm::ParticleCloudMpm<MPM_ITEMS, PARTICLE_TYPE>::setParticleForces(IT iter, IT end) {
    for (; iter != end; ++iter) {
        boost::tuple<unsigned, unsigned, double> t = *iter;
        const unsigned particleNum   = t.get<0>();
        const unsigned component     = t.get<1>();
        const double   value         = t.get<2>();
        particles_.at(particleNum)->storeParticleForce(component, value);
    }
}


// TAKES VERY LONG TIME (so included as an input file)
// void setInitialStressParticles(const double pSpacing, const VecDim G);

// evaluate initial stress of particles using k0

// template<typename MPM_ITEMS, typename PARTICLE_TYPE>
// void mpm::ParticleCloudMpm<MPM_ITEMS, PARTICLE_TYPE>::setInitialStressParticles(const double pSpacing, const VecDim G) {
//     unsigned numParticles = particles_.size();
//     VecDim X, Y;
//     double hMax = 0.0;
//     double hOverburden;

//     for (unsigned i = 0; i < numParticles; i++) {
//       X = particles_[i]->giveParticleCoordinates();
//       for (unsigned j = 0; j < numParticles; j++) {
// 	Y = particles_[j]->giveParticleCoordinates();
// 	if ((Y(0) > X(0) - pSpacing)&&(Y(0) < X(0) + pSpacing) && (Y(1) >= X(1)) && (Y(1) > hMax)) {
// 	  hMax = Y(1);
// 	}
//       }
//       hOverburden = hMax - X(1) + pSpacing;
//       particles_[i]->setInitialStress(hOverburden, G);
//       hMax = 0.0;
//     }
//     return;
// }

