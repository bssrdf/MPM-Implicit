// *********************************************************************************
// File: MaterialBaseMpm.hpp
//
// Details: Base class for Materials
//
//
// Author: Fehmi Cirak, California Institute of Technology
//         Krishna Kumar and Samila Bandara, University of Cambridge
//
// Version: 2.0 - 2013; 1.0 - 2005
// *********************************************************************************

#ifndef CARTESIAN_MATERIAL_MATERIALBASEMPM_H_
#define CARTESIAN_MATERIAL_MATERIALBASEMPM_H_

#include <iosfwd>
#include <string>

#include <boost/numeric/ublas/matrix.hpp>

#include <corlib/PropertiesParser.hpp>

namespace mpm {
namespace material {
class MaterialBaseMpm;
namespace ublas = boost::numeric::ublas;
}
}

#include "cartesian/src/MpmItems.hpp"


class mpm::material::MaterialBaseMpm {
public:
    MaterialBaseMpm();
    virtual ~MaterialBaseMpm();

protected:
    typedef ublas::bounded_matrix<double, 6, 6> Mat6x6_;
    typedef boost::numeric::ublas::bounded_vector<double, 6> Vec6x1_;
    typedef ublas::bounded_matrix<double, 2, 2> Mat2x2_;
    typedef ublas::bounded_matrix<double, 3, 3> Mat3x3_;

public:
    // Return instances of each material to each particle
    // virtual MaterialBaseMpm* getInstance() = 0;

    typedef mpm::MpmItems<mpm::DIM, mpm::DOF, mpm::NUMNODES, mpm::SHAPEFN> MpmItemsType;

    typedef mpm::ParticleSoilMpm<MpmItemsType>  SoilParticleHandle;
    typedef mpm::ParticleCloudMpm<MpmItemsType, SoilParticleHandle> ParticleCloudSoilHandle;
    typedef SoilParticleHandle* ParticleSoilPtr;
    typedef ParticleCloudSoilHandle* ParticleCloudSoilPtr;

    // return value of mass density
    virtual double density() = 0;

    // return value of phi (internal friction angle of soil)
    virtual double phi() = 0;

    // return value of porosity
    virtual double porosity() = 0;

    // compute elasticity tensor
    virtual void elasticityTensor(Mat6x6_& De) const = 0;

    // compute stress
    virtual void computeStress(const Vec6x1_& F, double& P, Vec6x1_& S, const ParticleCloudSoilPtr pCloudSoilPtr, const unsigned id) = 0;

    // compute stress for 3D
    virtual void computeStress3D(const Vec6x1_& F, double& P, Vec6x1_& S, const ParticleCloudSoilPtr pCloudSoilPtr, const unsigned id) = 0;

    // I/O
    void readMaterialParams(std::istream& is);
    void printMaterialParams(std::ostream& os);

protected:
    template <typename T>
    void registerVariable(const std::string& name, T& variable) {
        parser_->registerPropertiesVar(name, variable);
    }

private:
    corlib::PropertiesParser* parser_;

};

#include <cartesian/material/MaterialBaseMpm.ipp>

#endif  // CARTESIAN_MATERIAL_MATERIALBASEMPM_H_
