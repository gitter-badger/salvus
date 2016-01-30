//
// Created by Michael Afanasiev on 2016-01-29.
//

#ifndef SALVUS_OPTIONS_H
#define SALVUS_OPTIONS_H

#include <iosfwd>
#include <string>
#include "petscsys.h"

class Options {

    // Integer options.
    PetscInt mPolynomialOrder;

    // String options.
    std::string mMeshType;
    std::string mExodusMeshFile;
    std::string mElementShape;
    std::string mPhysicsSystem;

public:

    void setOptions();

    // Integer geters
    inline PetscInt PolynomialOrder() { return mPolynomialOrder; }

    // String geters
    inline std::string PhysicsSystem() { return mPhysicsSystem; }
    inline std::string ExodusMeshFile() { return mExodusMeshFile; }
    inline std::string MeshType() { return mMeshType; }
    inline std::string ElementShape() { return mElementShape; }


};


#endif //SALVUS_OPTIONS_H_H
