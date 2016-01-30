//
// Created by Michael Afanasiev on 2016-01-30.
//

#ifndef SALVUS_HYPERCUBE_H
#define SALVUS_HYPERCUBE_H

#include "../Element.h"
extern "C" {
#include "Square/Autogen/order4_square.h"
};

/*
 * Base class of an abstract four node square. The reference element is set up as below.
 *
 * (n2)______________(n3)
 * |                    |
 * |                    |
 * |                    |
 * |                    |
 * |                    |
 * |                    |
 * |                    |
 * |                    |
 * |                    |
 * (n0)______________(n1)
 *
 * (eta)
 *   ^
 *   |
 *   |
 *   |______> (eps)
*/


class Square : public Element {

    PetscReal n0(const PetscReal &eps, const PetscReal &eta);
    PetscReal n1(const PetscReal &eps, const PetscReal &eta);
    PetscReal n2(const PetscReal &eps, const PetscReal &eta);
    PetscReal n3(const PetscReal &eps, const PetscReal &eta);

    PetscReal dn0deps(const PetscReal &eta);
    PetscReal dn1deps(const PetscReal &eta);
    PetscReal dn2deps(const PetscReal &eta);
    PetscReal dn3deps(const PetscReal &eta);

    PetscReal dn0deta(const PetscReal &eps);
    PetscReal dn1deta(const PetscReal &eps);
    PetscReal dn2deta(const PetscReal &eps);
    PetscReal dn3deta(const PetscReal &eps);

    std::vector<PetscReal> mIntegrationPoints;

    // Fixed size eigen matrices for speed
    Eigen::Matrix<double,2,2> mJacobianBuffer;
    Eigen::Matrix<double,2,4> mVertexCoordinates;
    Eigen::Matrix<double,2,4> mJacobianMultiplier;

protected:

    PetscInt mNumberIntegrationPointsEps;
    PetscInt mNumberIntegrationPointsEta;

    std::vector<PetscReal> mIntegrationCoordinatesEps;
    std::vector<PetscReal> mIntegrationCoordinatesEta;

    // Local methods.
    void updateJacobian(PetscReal eps, PetscReal eta);

public:

    // Utility methods.
    virtual Element* clone() const = 0;

    // Superclass methods
    virtual void attachVertexCoordinates();
    virtual void attachIntegrationPoints();

    // Pure virtual methods.
    virtual void registerFieldVectors() = 0;
    virtual void constructStiffnessMatrix() = 0;

};


#endif //SALVUS_HYPERCUBE_H
