//
// Created by Michael Afanasiev on 2016-01-27.
//

#ifndef SALVUS_MESH_H
#define SALVUS_MESH_H


#include <iosfwd>
#include <string>
#include <petscdmtypes.h>
#include <vector>
#include "Options.h"

class Mesh {

    PetscInt mNumberElementsLocal;
    PetscInt mNumberDimensions;

    std::string mExodusFileName;

    DM mDistributedMesh;

public:

    static Mesh *factory(Options options);

    void read(Options options);
    void setupGlobalDof(PetscInt number_dof_vertex, PetscInt number_dof_edge,
                        PetscInt number_dof_face, PetscInt number_dof_volume,
                        PetscInt number_dimensions);

    std::vector<PetscReal> ElementVertices(const PetscInt element_number);

    // Integer getattr.
    inline PetscInt NumberElementsLocal() { return mNumberElementsLocal; }

    // Distributed mesh getattr.
    inline DM &DistributedMesh() { return mDistributedMesh; }



};

class Exodus : public Mesh {


public:


};

class Simple : public Mesh {

};


#endif //SALVUS_MESH_H
