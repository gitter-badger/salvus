#include <iostream>
#include <mpi.h>
#include <petscsys.h>

static constexpr char help[] = "Welcome to salvus.";

int main(int argc, char *argv[]) {
    PetscInitialize(&argc, &argv, NULL, help);
    std::cout << "Hello world." << std::endl;
    PetscFinalize();
}