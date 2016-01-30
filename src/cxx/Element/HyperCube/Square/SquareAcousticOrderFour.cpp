//
// Created by Michael Afanasiev on 2016-01-30.
//

#include <petscdm.h>
#include <iostream>
#include <openmpi/ompi/mpi/cxx/mpicxx.h>
#include "SquareAcousticOrderFour.h"
#include "../../../Options.h"

SquareAcousticOrderFour::SquareAcousticOrderFour(Options options) {

    // Basic properties.
    mElementShape = options.ElementShape();
    mPhysicsSystem = options.PhysicsSystem();
    mPolynomialOrder = options.PolynomialOrder();
    mNumberDimensions = 2;
    mNumberVertex = 4;

    // Gll points.
    mNumberDofVertex = 1;
    mNumberDofEdge = mPolynomialOrder - 1;
    mNumberDofFace = (mPolynomialOrder - 1) * (mPolynomialOrder - 1);
    mNumberDofVolume = 0;
    mNumberIntegrationPointsEps = 5;
    mNumberIntegrationPointsEta = 5;
    mNumberIntegrationPoints = 25;

    mIntegrationCoordinatesEps = {-1.0, -0.6546536707, 0.0, 0.6546536707, 1.0};
    mIntegrationCoordinatesEta = {-1.0, -0.6546536707, 0.0, 0.6546536707, 1.0};

}

void SquareAcousticOrderFour::registerFieldVectors() {

    PetscReal zero = 0;
    DMCreateLocalVector(mDistributedMesh, &mDisplacementLocal);
    DMCreateLocalVector(mDistributedMesh, &mAccelerationLocal);
    DMCreateLocalVector(mDistributedMesh, &mVelocityLocal);
    VecSet(mDisplacementLocal, zero);
    VecSet(mAccelerationLocal, zero);
    VecSet(mVelocityLocal, zero);

    DMCreateGlobalVector(mDistributedMesh, &mDisplacementGlobal);
    DMCreateGlobalVector(mDistributedMesh, &mAccelerationGlobal);
    DMCreateGlobalVector(mDistributedMesh, &mVelocityGlobal);
    VecSet(mDisplacementGlobal, zero);
    VecSet(mAccelerationGlobal, zero);
    VecSet(mVelocityGlobal, zero);

    PetscObjectSetName((PetscObject) mDisplacementGlobal, "displacement");
    PetscObjectSetName((PetscObject) mAccelerationGlobal, "acceleration");
    PetscObjectSetName((PetscObject) mVelocityGlobal, "velocity");

}

void SquareAcousticOrderFour::constructStiffnessMatrix() {

    for (int i = 0; i < mNumberIntegrationPointsEta; i++) {
        for (int j = 0; j < mNumberIntegrationPointsEps; j++) {

            PetscReal eps = mIntegrationCoordinatesEps[j];
            PetscReal eta = mIntegrationCoordinatesEta[i];
            updateJacobian(eps, eta);
        }
    }
}
