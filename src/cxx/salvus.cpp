#include <iostream>
#include <mpi.h>
#include <petscsys.h>
#include "Solver.h"
#include "Mesh.h"

static constexpr char help[] = "Welcome to salvus.";

int main(int argc, char *argv[]) {

    PetscInitialize(&argc, &argv, NULL, help);

    Solver *solver = Solver::factory("time_domain");
    Mesh *mesh = Mesh::factory("exodus");

    std::cout << "Hello world." << std::endl;

    PetscFinalize();
}
