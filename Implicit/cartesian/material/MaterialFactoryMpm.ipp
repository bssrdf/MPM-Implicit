// *********************************************************************************
// File: MaterialFactoryMpm.hpp
//
// Details: Implementation of container class for Materials
//
// Author: Fehmi Cirak, California Institute of Technology
//
// Version: 1.0 - 2005
// *********************************************************************************

#include <cartesian/material/MaterialFactoryMpm.hpp>
#include <utility>
#include <cassert>

using namespace mpm::material;

MaterialFactoryMpm* MaterialFactoryMpm::instance_ = NULL;


MaterialFactoryMpm* MaterialFactoryMpm::instance() {
    if (!instance_)
        instance_ = new MaterialFactoryMpm();

    return instance_;
}


void MaterialFactoryMpm::destroy() {
    callBacks_.clear();
    delete instance_;
}


bool MaterialFactoryMpm::registerMaterial(std::string materialType,
        CreateMaterialCallBack cb) {
    return callBacks_.insert(CallBackMap_::value_type(materialType,
                             cb)).second;
}


MaterialBaseMpm* MaterialFactoryMpm::createMaterial(std::string& materialType) const {
    CallBackMap_::const_iterator it = callBacks_.find(materialType);

    if (it == callBacks_.end())
        assert(!"Unknown material type");

    return (it->second());
}
