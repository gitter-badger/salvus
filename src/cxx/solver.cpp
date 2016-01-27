//
// Created by Michael Afanasiev on 2016-01-27.
//

#include <ostream>
#include <iostream>
#include "Solver.h"


Solver *Solver::factory(std::string solver_type) {
    if (solver_type == "time_domain") {
        return new TimeDomain;
    } else if (solver_type == "frequency_domain") {
        return new FrequencyDomain;
    } else {
        return nullptr;
    }
}

void TimeDomain::getType() {

    std::cout << "I AM A TIME DOMAIN SOLVER." << std::endl;

}

void FrequencyDomain::getType() {

    std::cout << "I AM A FREQUENCY DOMAIN SOLVER." << std::endl;

}