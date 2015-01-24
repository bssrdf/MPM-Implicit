// *********************************************************************************
// File: MaterialBaseMpm.hpp
//
// Details: Implementation of base class for Materials
//
// Author: Fehmi Cirak, California Institute of Technology
//
// Version: 1.0 - 2005
// *********************************************************************************

#include <cassert>

#include <cartesian/material/MaterialBaseMpm.hpp>

using namespace mpm::material;


MaterialBaseMpm::MaterialBaseMpm() {
    parser_ = new corlib::PropertiesParser;
}


MaterialBaseMpm::~MaterialBaseMpm() {
    assert(parser_ != NULL);
    delete parser_;
}


void MaterialBaseMpm::readMaterialParams(std::istream& is) {
    char c;
    is >> c;
    if (c == '{') {
        (is >> std::ws).get(c);
        while (c != '}') {
            is.putback(c);

            corlib::skip_comment(is);
            (is >> std::ws).get(c);
            if (c== '}') return;

            is.putback(c);

            parser_->readVariable(is);
            (is >> std::ws).get(c);
        }
    }

    return;
}


void MaterialBaseMpm::printMaterialParams(std::ostream& os) {
    os << '{' << "\n";
    parser_->printValues(os, " ", " ");
    os << " } \n";

    return;
}
