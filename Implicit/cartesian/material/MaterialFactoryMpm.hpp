// *********************************************************************************
// File: MaterialBaseMpm.hpp
//
// Details: Base class for Materials
//
//
// Author: Fehmi Cirak, California Institute of Technology
//
// Version: 1.0 - 2005
// *********************************************************************************

#ifndef CARTESIAN_MATERIAL_MATERIALFACTORYMPM_H_
#define CARTESIAN_MATERIAL_MATERIALFACTORYMPM_H_

#include <map>
#include <cassert>
#include <string>


namespace mpm {
namespace material {
class MaterialBaseMpm;
class MaterialFactoryMpm;
}
}

class mpm::material::MaterialFactoryMpm {
public:
    typedef MaterialBaseMpm* (*CreateMaterialCallBack)();

private:
    typedef std::map<std::string, CreateMaterialCallBack> CallBackMap_;

public:
    bool registerMaterial(std::string materialType, CreateMaterialCallBack cb);
    MaterialBaseMpm* createMaterial(std::string& materialType) const;

    // this is a singleton class
    static MaterialFactoryMpm * instance();
    void destroy();

private:
    MaterialFactoryMpm() {}
    ~MaterialFactoryMpm() {}
    static MaterialFactoryMpm *instance_;

    // copy constructor and assignment operator
private:
    MaterialFactoryMpm(const MaterialFactoryMpm& mf);
    const MaterialFactoryMpm& operator = (const MaterialFactoryMpm& mf);

    CallBackMap_ callBacks_;
};

#include <cartesian/material/MaterialFactoryMpm.ipp>

#endif  // CARTESIAN_MATERIAL_MATERIALFACTORYMPM_H_
