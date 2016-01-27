//
// Created by Michael Afanasiev on 2016-01-27.
//

#include "Mesh.h"

Mesh *Mesh::factory(std::string mesh_type) {
    if (mesh_type == "exodus") {
        return new Exodus;
    } else {
        return nullptr;
    }

}
