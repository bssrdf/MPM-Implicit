// *********************************************************************************
// File: ParticleCloudMpm.hpp
//
// Details: Base class container for the MPM particles
//          Has containers with pointers to the particles
//
// Dependency: MPM_ITEMS contains all template parameters
//
// Author: Krishna Kumar and Samila Bandara, University of Cambridge
//
// Version: 2.0 revised by Krishna
// *********************************************************************************

#ifndef CARTESIAN_SRC_PARTICLECLOUDMPM_H_
#define CARTESIAN_SRC_PARTICLECLOUDMPM_H_

#include <vector>
#include <string>
#include <iostream>
#include <algorithm>

#include "cartesian/src/MpmItems.hpp"
#include "cartesian/src/Enums.hpp"


namespace mpm {
template<typename MPM_ITEMS, typename PARTICLE_TYPE>
class ParticleCloudMpm;
namespace ublas = boost::numeric::ublas;
}


// \brief Base class container for the MPM particles
// \details Has containers with pointers to the particles
// \tparam MPM_ITEMS   the main template file which includes all the classes
// \tparam PARTICLE_TYPE  the type of particle (soil/water) to be used

template<typename MPM_ITEMS, typename PARTICLE_TYPE>
class mpm::ParticleCloudMpm {
public:
    typedef MPM_ITEMS MpmItems;
    typedef PARTICLE_TYPE ParticleType;

    static const mpm::ParticleTypeTag myType = ParticleType::myType;

    static const unsigned dof = MpmItems::dof;
    static const unsigned dim = MpmItems::dim;

private:
    typedef ParticleType* ParticlePtr_;

public:
    typedef std::vector<ParticlePtr_> ParticleVec;
    typedef typename ParticleVec::const_iterator ParticleConstIterator;
    typedef ublas::bounded_vector<double, dof> VecDof;
    typedef ublas::bounded_vector<double, dim> VecDim;
    typedef ublas::bounded_vector<double, 6> Vec6x1;

    // Constructor with SMF file
    ParticleCloudMpm(std :: istream& prcld);

    // Destructor which clears the containers
    ~ParticleCloudMpm();

    // give particle vector
    ParticleVec giveParticles() {
        return particles_;
    }

    // Give particle pointer for given index
    ParticlePtr_ giveParticlePtr(const unsigned ind) const {
        return particles_[ind];
    }

    // evaluate initial stress of particles using k0
    void setInitialStressParticles(std::istream& initStress);

    // Assign Material Type to particles
    void associateParticleToMaterial(std::istream& matType);

    // write soil particle data into an output file
    void writeParticleCloudVelocityData(const unsigned int step, std::ostream& out);
    void writeParticleCloudStressData(const unsigned int step, std::ostream& out);
    void writeParticleCloudStrainData(const unsigned int step, std::ostream& out);
    void writeParticleCloudRateOfStrainI2Data(const unsigned int step, std::ostream& out);

    // Generic application of functors to particle containers
    // Read particle traction pressures
    void readTractionBoundary(std :: istream& trbn);
    // Read particle external forces
    void readParticleForces(std :: istream& pfrcs);                 // NOT IN USE
    // Set the constraints given access to tuples
    template<typename IT> void setParticleForces(IT iter, IT end);  // NOT IN USE

    // Generic application of functors to particle containers
    // Iterate over particles
    template<typename OP>  OP iterateOverParticles(OP op);
    template<typename OP>  OP iterateOverParticlesParallel(OP op);
    // Iterate over particle with predicate
    template<typename OP, typename PRED>
    OP iterateOverParticlesWithPredicate(OP op, PRED p);


    // Const access to particle containers
    // Give const access to begin of particale container
    ParticleConstIterator particlesBegin() const {
        return particles_.begin();
    }
    // Give const access to end of particle container
    ParticleConstIterator particlesEnd() const {
        return particles_.end();
    }
    // Give const access to specific particle
    ParticlePtr_ particleIterator(const unsigned numParticle) const;
    // Give access to a particle pointer
    ParticlePtr_ particleIterator(const unsigned numParticle) {
        return particles_[numParticle];
    }

    // Size queries
    // Return number of particles
    unsigned numParticles() const {
        return particles_.size();
    }
    // Dummy function for composite clouds
    unsigned numParticleCoulds() const {
        return 1.;
    }

protected:
    ParticleVec particles_;  // vector of particle pointers
};

#include "ParticleCloudMpm.ipp"


#endif  // CARTESIAN_SRC_PARTICLECLOUDMPM_H_
