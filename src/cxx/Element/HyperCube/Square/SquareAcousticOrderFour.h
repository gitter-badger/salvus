//
// Created by Michael Afanasiev on 2016-01-30.
//

#ifndef SALVUS_SQUAREACOUSTICORDERFOUR_H
#define SALVUS_SQUAREACOUSTICORDERFOUR_H


#include <petscvec.h>
#include "../../Element.h"
#include "../../../Options.h"
#include "../Square.h"

class SquareAcousticOrderFour: public Square {

    Vec mDisplacementLocal;
    Vec mAccelerationLocal;
    Vec mVelocityLocal;

    Vec mDisplacementGlobal;
    Vec mAccelerationGlobal;
    Vec mVelocityGlobal;

public:

    SquareAcousticOrderFour(Options options);

    virtual Element *clone() const { return new SquareAcousticOrderFour(*this); }

    virtual void registerFieldVectors();
    virtual void constructStiffnessMatrix();
    virtual void interpolateMaterialProperties();

};


#endif //SALVUS_SQUAREACOUSTICORDERFOUR_H
