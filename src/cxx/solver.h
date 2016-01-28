//
// Created by Michael Afanasiev on 2016-01-27.
//

#ifndef SALVUS_SOLVER_H
#define SALVUS_SOLVER_H


#include <iosfwd>
#include <string>

class Solver {

public:

    static Solver *factory(std::string solver_type);

    virtual void getType() = 0;

};

class TimeDomain : public Solver {

public:

    void getType();

};

class FrequencyDomain : public Solver {

public:

    void getType();
};


#endif //SALVUS_SOLVE
