//
// Created by Michael Afanasiev on 2016-01-29.
//

#include "Utilities.h"

void utilities::print_from_root_mpi(const std::string msg) {
    if (MPI::COMM_WORLD.Get_rank() == 0) {
        std::cout << msg << std::endl;
    }
}

std::vector<double> utilities::broadcastStdVecFromRoot(std::vector<double> &send_buffer) {

    int root = 0;
    int int_size = 1;
    std::vector<double> receive_buffer;

    int length;
    if (!MPI::COMM_WORLD.Get_rank()) { length = send_buffer.size(); }
    MPI::COMM_WORLD.Bcast(&length, int_size, MPI::INT, root);

    receive_buffer.resize(length);
    if (!MPI::COMM_WORLD.Get_rank()) { receive_buffer = send_buffer; }
    MPI::COMM_WORLD.Bcast(receive_buffer.data(), length, MPI::DOUBLE, root);
    return receive_buffer;

}