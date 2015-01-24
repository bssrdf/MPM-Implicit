// *********************************************************************************
// File: MaterialContainerMpm.hpp
//
// Details: Material Container Class
//
// Author: Fehmi Cirak, California Institute of Technology
//         Kosala Bandara, Computational Structural Mechanics Lab,
//         University of Cambridge
//
// Version: 1.0 - 2010
// *********************************************************************************

#ifndef CARTESIAN_MATERIAL_MATERIALCONTAINERMPM_H_
#define CARTESIAN_MATERIAL_MATERIALCONTAINERMPM_H_

#include <cartesian/material/MaterialBaseMpm.hpp>
#include <cartesian/material/MaterialFactoryMpm.hpp>

#include <corlib/misc.hpp>
#include <corlib/PropertiesParser.hpp>

#include <vector>
#include <cassert>
#include <iosfwd>

namespace mpm {
namespace material {
class MaterialContainerMpm;
class MaterialBaseMpm;
}
}

class mpm::material::MaterialContainerMpm {
public:
    // this is a singleton class
    static MaterialContainerMpm* instance();
    void destroy();

    void readMaterialStream(std::istream& is);
    void printMaterialStream(std::ostream& os);

    void addMaterial(MaterialBaseMpm* mat) {
        materials_.push_back(mat);
    }

    // acessor
    MaterialBaseMpm* getMaterial(const unsigned i) const {
        assert(i < materials_.size());
        return materials_[i];
    }

    std::vector<MaterialBaseMpm*> getMaterials() {
        return materials_;
    }

    // this is a singleton class
private:
    MaterialContainerMpm() {}
    ~MaterialContainerMpm() {}
    static MaterialContainerMpm *instance_;

    // copy constructor and assignment operator
private:
    MaterialContainerMpm(const MaterialContainerMpm& mc);
    const MaterialContainerMpm& operator = (const MaterialContainerMpm& mc);

private:
    typedef std::vector<MaterialBaseMpm*> MaterialCont_;
    typedef MaterialCont_::iterator MaterialIt_;

    MaterialCont_ materials_;
};
#include <cartesian/material/MaterialContainerMpm.ipp>

#endif  // CARTESIAN_MATERIAL_MATERIALCONTAINERMPM_H_

