// *********************************************************************************
// File: MpmItems.hpp
//
// Details: Contains all template parameters.
//          Items in mpm data structure for entire mpm class
//          All the classes in mpm:: are templated by this
//
// Author: Krishna Kumar and Samila Bandara, University of Cambridge
//
// Version: 2.0 Revised by Krishna, 1.0 - Samila Bandara
// *********************************************************************************

#ifndef CARTESIAN_SRC_MPMITEMS_H_
#define CARTESIAN_SRC_MPMITEMS_H_

 #define _MPM2D_ //Toggle between MPM2D and 3D Analysis

#ifdef _MPM2D_
namespace mpm {
const unsigned DIM          = 2;
const unsigned DOF          = 2;
const unsigned NUMNODES     = 4;
const corlib::shape SHAPEFN = corlib::QUADRILATERAL;
}
#else
namespace mpm {
const unsigned DIM          = 3;
const unsigned DOF          = 3;
const unsigned NUMNODES     = 8;
const corlib::shape SHAPEFN = corlib::HEXAHEDRON;
}
#endif

#include "corlib/Shape.hpp"

#include "cartesian/src/ElementMpm.hpp"
#include "cartesian/src/NodeMpm.hpp"
#include "cartesian/src/MeshMpm.hpp"
#include "cartesian/src/SubMeshMpm.hpp"
#include "cartesian/src/ParticleMpm.hpp"
#include "cartesian/src/ParticleCloudMpm.hpp"
#include "cartesian/src/ParticleSoilMpm.hpp"
#include "cartesian/src/ParticleWaterMpm.hpp"
#include "cartesian/src/MatrixMpm.hpp"
#include "cartesian/src/SolverMpm.hpp"
#include "cartesian/material/MaterialBaseMpm.hpp"
#include "cartesian/src/Enums.hpp"

namespace mpm {
template<unsigned DIM, unsigned DOF, unsigned NNODES, corlib::shape SHAPE>
class MpmItems;
}

// \brief Items in mpm data structure for entire mpm class
// \details All the classes in mpm:: are templated by this
// \tparam DIM  the dimension of the problem
// \tparam DOF  number of degrees of freedom
// \tparam NNODES  number of nodes (=4 since rectanglular mesh)

template<unsigned DIM, unsigned DOF, unsigned NNODES, corlib::shape SHAPE>
class mpm::MpmItems {
public:
    static const unsigned dim = DIM;
    static const unsigned dof = DOF;
    static const unsigned nNodes = NNODES;
    static const corlib::shape myShape = SHAPE;

    typedef typename mpm::MpmItems<dim, dof, nNodes, myShape> MpmItemsHandle;

    typedef typename mpm::NodeMpm<MpmItemsHandle> NodeHandle;
    typedef typename mpm::ElementMpm<MpmItemsHandle> ElementHandle;
    typedef typename mpm::ParticleMpm<MpmItemsHandle> ParticleHandle;
    typedef typename mpm::ParticleSoilMpm<MpmItemsHandle> ParticleSoilHandle;
    typedef typename mpm::ParticleWaterMpm<MpmItemsHandle> ParticleWaterHandle;
    typedef typename mpm::SubMeshMpm<MpmItemsHandle> SubMeshHandle;
    typedef typename mpm::MeshMpm<MpmItemsHandle> MeshHandle;
    typedef typename mpm::ParticleCloudMpm<MpmItemsHandle, ParticleSoilHandle> ParticleCloudSoilHandle;
    typedef typename mpm::ParticleCloudMpm<MpmItemsHandle, ParticleWaterHandle> ParticleCloudWaterHandle;
    typedef typename mpm::MatrixMpm<MpmItemsHandle> MatrixHandle;
    typedef typename mpm::SolverMpm<MpmItemsHandle> SolverHandle;

    typedef typename mpm::material::MaterialBaseMpm MaterialHandle;

    typedef NodeHandle* NodePtr;
    typedef ElementHandle* ElementPtr;
    typedef ParticleHandle* ParticlePtr;
    typedef ParticleSoilHandle* ParticleSoilPtr;
    typedef ParticleWaterHandle* ParticleWaterPtr;
    typedef SubMeshHandle* SubMeshPtr;
    typedef MeshHandle* MeshPtr;
    typedef ParticleCloudSoilHandle* ParticleCloudSoilPtr;
    typedef ParticleCloudWaterHandle* ParticleCloudWaterPtr;
    typedef MatrixHandle* MatrixPtr;
    typedef MaterialHandle* MaterialPtr;
    typedef SolverHandle* SolverPtr;
};

#endif  // CARTESIAN_SRC_MPMITEMS_H_

