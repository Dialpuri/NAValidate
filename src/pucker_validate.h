//
// Created by Jordan Dialpuri on 18/12/2023.
//

#ifndef PUCKER_VALIDATE_H
#define PUCKER_VALIDATE_H

#include "validate-util.h"

struct PuckerResult {
    PuckerResult(float i_P, float i_theta_m) {
        P = i_P;
        theta_m = i_theta_m;
    }

    float P;
    float theta_m;

    static PuckerResult null() {
        return PuckerResult(-100,-100);
    }

    static bool is_null(const PuckerResult& result) {
        return result.P == -100 && result.theta_m == -100;
    }

    bool is_null() const {
        return P == -100 && theta_m == -100;
    }
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

class PuckerValidate {
public:
    static PuckerResult calculate_pucker(const clipper::MMonomer& monomer);
    static PuckerType classify_pucker(const PuckerResult& result);
};



#endif //PUCKER_VALIDATE_H
