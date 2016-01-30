//
// Created by Michael Afanasiev on 2016-01-29.
//

#ifndef SALVUS_UTILITIES_H
#define SALVUS_UTILITIES_H

#include <iosfwd>
#include <string>
#include <openmpi/ompi/mpi/cxx/mpicxx.h>

namespace utilities {
    void print_from_root_mpi(const std::string msg);
}

#endif //SALVUS_UTILITIES_H
