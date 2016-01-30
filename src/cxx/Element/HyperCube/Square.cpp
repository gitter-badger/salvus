//
// Created by Michael Afanasiev on 2016-01-30.
//

#include <iostream>
#include <petscdm.h>
#include <petscdmplex.h>
#include <assert.h>
#include <openmpi/ompi/mpi/cxx/mpicxx.h>
#include "Square.h"

void Square::attachVertexCoordinates() {

    Vec coordinates_local;
    PetscInt coordinate_buffer_size;
    PetscSection coordinate_section;
    PetscReal *coordinates_buffer = NULL;

    DMGetCoordinatesLocal(mDistributedMesh, &coordinates_local);
    DMGetCoordinateSection(mDistributedMesh, &coordinate_section);
    DMPlexVecGetClosure(mDistributedMesh, coordinate_section, coordinates_local, LocalElementNumber(),
                        &coordinate_buffer_size, &coordinates_buffer);
    std::vector<PetscReal> coordinates_element(coordinates_buffer, coordinates_buffer+coordinate_buffer_size);
    DMPlexVecRestoreClosure(mDistributedMesh, coordinate_section, coordinates_local, LocalElementNumber(),
                            &coordinate_buffer_size, &coordinates_buffer);

    // Reorder to desired vertex ordering.
    std::vector<double> vertex_coordinates_ordered ;
    std::vector<PetscInt> mapping_to_reference_element {6, 7, 0, 1, 4, 5, 2, 3};
    for (int i = 0; i < mNumberVertex; i++) {
        mVertexCoordinates(0,i) = coordinates_element[mapping_to_reference_element[mNumberDimensions*i+0]];
        mVertexCoordinates(1,i) = coordinates_element[mapping_to_reference_element[mNumberDimensions*i+1]];
    }

}

void Square::attachIntegrationPoints() {

    assert(NumberIntegrationPoints() == mNumberIntegrationPointsEps*mNumberIntegrationPointsEta);

    int point = 0;
    std::vector<PetscReal> ni(mVertexCoordinates.size());
    mIntegrationPoints.resize(mNumberIntegrationPoints * mNumberDimensions);
    for (int i = 0; i < mNumberIntegrationPointsEta; i++) {
        for (int j = 0; j < mNumberIntegrationPointsEps; j++) {

            PetscReal eps = mIntegrationCoordinatesEps[j];
            PetscReal eta = mIntegrationCoordinatesEta[i];

            ni[0] = n0(eps, eta);
            ni[1] = n1(eps, eta);
            ni[2] = n2(eps, eta);
            ni[3] = n3(eps, eta);

            for (int k = 0; k < mNumberVertex; k++) {
                mIntegrationPoints[point + 0] += ni[k] * mVertexCoordinates(0, k);
                mIntegrationPoints[point + 1] += ni[k] * mVertexCoordinates(1, k);
            }

            point += mNumberDimensions;
        }
    }
}

PetscReal Square::n0(const PetscReal &eps, const PetscReal &eta) {
    return 0.25 * (1.0 - eps) * (1.0 - eta);
}

PetscReal Square::n1(const PetscReal &eps, const PetscReal &eta) {
    return 0.25 * (1.0 + eps) * (1.0 - eta);
}

PetscReal Square::n2(const PetscReal &eps, const PetscReal &eta) {
    return 0.25 * (1.0 - eps) * (1.0 + eta);
}

PetscReal Square::n3(const PetscReal &eps, const PetscReal &eta) {
    return 0.25 * (1.0 + eps) * (1.0 + eta);
}

PetscReal Square::dn0deps(const PetscReal &eta) {
    return (-1) * (1 - eta) / 4.0;
}

PetscReal Square::dn1deps(const PetscReal &eta) {
    return (+1) * (1 - eta) / 4.0;
}

PetscReal Square::dn2deps(const PetscReal &eta) {
    return (-1) * (1 + eta) / 4.0;
}

PetscReal Square::dn3deps(const PetscReal &eta) {
    return (+1) * (1 + eta) / 4.0;
}

PetscReal Square::dn0deta(const PetscReal &eps) {
    return (1 - eps) * -1.0 / 4.0;
}

PetscReal Square::dn1deta(const PetscReal &eps) {
    return (1 + eps) * -1.0 / 4.0;
}

PetscReal Square::dn2deta(const PetscReal &eps) {
    return (1 - eps) * 1.0 / 4.0;
}

PetscReal Square::dn3deta(const PetscReal &eps) {
    return (1 + eps) * 1.0 / 4.0;

}

void Square::updateJacobian(PetscReal eps, PetscReal eta) {


    // Set local values.
    mJacobianMultiplier(0,0) = dn0deps(eta);
    mJacobianMultiplier(0,1) = dn1deps(eta);
    mJacobianMultiplier(0,2) = dn2deps(eta);
    mJacobianMultiplier(0,3) = dn3deps(eta);
    mJacobianMultiplier(1,0) = dn0deta(eps);
    mJacobianMultiplier(1,1) = dn1deta(eps);
    mJacobianMultiplier(1,2) = dn2deta(eps);
    mJacobianMultiplier(1,3) = dn3deta(eps);

    // Jacobian at (eps, eta).
    mJacobianBuffer = mJacobianMultiplier * mVertexCoordinates.transpose();

}
