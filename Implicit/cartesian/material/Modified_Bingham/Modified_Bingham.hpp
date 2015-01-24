// *********************************************************************************
// File: Modified_Bingham.hpp
//
// Details: Modified Bingham Material Model
//
// Dependency: MaterialBase, ParticleCloudMpm and MpmItems
//
// Author: Krishna Kumar, University of Cambridge
//
// Version: 1.0
// *********************************************************************************

#ifndef CARTESIAN_MATERIAL_MODIFIED_BINGHAM_MODIFIED_BINGHAM_H_
#define CARTESIAN_MATERIAL_MODIFIED_BINGHAM_MODIFIED_BINGHAM_H_

#include <boost/numeric/ublas/matrix.hpp>

#include "cartesian/material/MaterialBaseMpm.hpp"

namespace mpm {
namespace material {
class Modified_Bingham;
namespace ublas = boost::numeric::ublas;
}
}


class mpm::material::Modified_Bingham : public mpm::material::MaterialBaseMpm {
private:
    Modified_Bingham();
    Modified_Bingham(const Modified_Bingham & other);
    Modified_Bingham& operator = (const Modified_Bingham& other);

private:
    typedef ublas::bounded_matrix<double, 6, 6> Mat6X6_;
    typedef ublas::bounded_vector<double, 6> Vec6x1_;
    typedef ublas::bounded_vector<double, 3> Vec3x1_;
    typedef ublas::bounded_matrix<double, 2, 2> Mat2x2_;
    typedef ublas::bounded_matrix<double, 3, 3> Mat3x3_;


public:
    ~Modified_Bingham();

    static MaterialBaseMpm* createMaterial() {
        return new Modified_Bingham();
    }

    // To Create seperate material instances for each particle: Don't use
    Modified_Bingham* getInstance() {

        Modified_Bingham* pModified_Bingham = new Modified_Bingham;
        pModified_Bingham->youngModulus_ = this->youngModulus_;
        pModified_Bingham->poissonRatio_ = this->poissonRatio_;
        pModified_Bingham->density_ = this->density_;
        pModified_Bingham->tau0_ = this->tau0_;
        pModified_Bingham->mu_ = this->mu_;

        //    pModified_Bingham->phi_ = this->phi_;
        //    pModified_Bingham->porosity_ = this->porosity_;

        return pModified_Bingham;
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

    // compute stress for 3D
    void computeStress3D(const Vec6x1_& F, double& P, Vec6x1_& S, const ParticleCloudSoilPtr pCloudSoilPtr, const unsigned id);

private:
    // send sign
    double sign_(double var) {}

private:
    double youngModulus_;
    double poissonRatio_;
    double density_;
    double tau0_;
    double mu_;
    double M_;
    double strainCutoff_;
};

#include "Modified_Bingham.ipp"

#endif  // CARTESIAN_MATERIAL_MODIFIED_BINGHAM_MODIFIED_BINGHAM_H_
