//
// Created by Jordan Dialpuri on 18/12/2023.
//

#ifndef PUCKER_VALIDATE_H
#define PUCKER_VALIDATE_H

#include "validate.h"


struct PuckerResult {
    float P;
    float theta_m;
};

//http://enzyme13.bt.a.u-tokyo.ac.jp/CP/altonas.html
enum PuckerConformation {
    _3T2, __3E, _3T4, __E4, _OT4, __OE, _OT1, __E1, _2T1, __2E,
    _2T3, __E3, _4T3, __4E, _4TO, __EO, _1TO, __1E, _1T2, __E2
};

struct PuckerType {
    PuckerConformation conformation;

    std::string to_string() const;
};

class PuckerValidate: public Validate {

};



#endif //PUCKER_VALIDATE_H
