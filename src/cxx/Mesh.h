//
// Created by Michael Afanasiev on 2016-01-27.
//

#ifndef SALVUS_MESH_H
#define SALVUS_MESH_H


#include <iosfwd>
#include <string>

class Mesh {

public:

    static Mesh *factory(std::string mesh_type);

};

class Exodus : public Mesh {

};

class Simple : public Mesh {

};


#endif //SALVUS_MESH_H
