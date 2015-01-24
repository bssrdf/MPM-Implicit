// *********************************************************************************
// File: Bingham.hpp
//
// Details: Bingham Material Model
//
// Dependency: MaterialBase, ParticleCloudMpm and MpmItems
//
// Author: Krishna Kumar, University of Cambridge
//
// Version: 1.0
// *********************************************************************************

#ifndef CARTESIAN_MATERIAL_BINGHAM_BINGHAM_H_
#define CARTESIAN_MATERIAL_BINGHAM_BINGHAM_H_

#include <boost/numeric/ublas/matrix.hpp>

#include "cartesian/material/MaterialBaseMpm.hpp"
#include "cartesian/src/MpmItems.hpp"
#include "cartesian/src/ParticleCloudMpm.hpp"

namespace mpm {
namespace material {
class Newtonian;
namespace ublas = boost::numeric::ublas;
}
}


class mpm::material::Newtonian : public mpm::material::MaterialBaseMpm {
private:
    Newtonian();
    Newtonian(const Newtonian& other);
    // Assignement operator
    Newtonian& operator = (const Newtonian& other);

private:
    typedef ublas::bounded_matrix<double, 6, 6> Mat6X6_;
    typedef ublas::bounded_vector<double, 6> Vec6x1_;
    typedef ublas::bounded_vector<double, 3> Vec3x1_;
    typedef ublas::bounded_matrix<double, 2, 2> Mat2x2_;
    typedef ublas::bounded_matrix<double, 3, 3> Mat3x3_;

public:
    ~Newtonian();

    static MaterialBaseMpm * createMaterial() {
        return new Newtonian();
    }

    // To Create seperate material instances for each particle: Don't use
    Newtonian* getInstance() {
        // Copy initial values
        Newtonian* pNewtonian = new Newtonian;
        pNewtonian->youngModulus_ = this->youngModulus_;
        pNewtonian->poissonRatio_ = this->poissonRatio_;
        pNewtonian->density_ = this->density_;
        //    pBingham->phi_ = this->phi_;
        //    pBingham->porosity_ = this->porosity_;
        pNewtonian->tau0_ = this->tau0_;
        pNewtonian->mu_ = this->mu_;
        return pNewtonian;
    }

    // return density
    double density() {
        return density_;
    }

    // return internal friction angle of soil (phi)
    double phi() {
        return 0.;
    }

    // return porosity of soil (n)
    double porosity() {
        return 0.;
    }

    // return elasticity tensor
    void elasticityTensor(Mat6x6_& De) const;

    // compute stress
    void computeStress(const Vec6x1_& F, double& P, Vec6x1_& S, const ParticleCloudSoilPtr pCloudSoilPtr, const unsigned id);

    // compute stress 3D
    void computeStress3D(const Vec6x1_& F, double& P, Vec6x1_& S, const ParticleCloudSoilPtr pCloudSoilPtr, const unsigned id);

private:
    // send sign
    double sign_(double var);

private:
    double youngModulus_;
    double poissonRatio_;
    double density_;
    //  double phi_;
    //  double porosity_;
    double tau0_;
    double mu_;
    double strainCutoff_;

};

#include "Newtonian.ipp"

#endif  // CARTESIAN_MATERIAL_BINGHAM_BINGHAM_H_
