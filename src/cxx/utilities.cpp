//
// Created by Michael Afanasiev on 2016-01-29.
//

#include "Utilities.h"

void utilities::print_from_root_mpi(const std::string msg) {
    if (MPI::COMM_WORLD.Get_rank() == 0) {
        std::cout << msg << std::endl;
    }
}