//
// Created by Jordan Dialpuri on 18/12/2023.
//

#ifndef BASEVALIDATION_H
#define BASEVALIDATION_H
#include <string>
#include "clipper/clipper-minimol.h"
#include "../util/matrix.h"
#include "../util/validate-util.h"

struct BaseConformationResult {
    float chi;
    std::string conformation;
};


class BaseConformationValidation {
public:
    static double calculate_chi(const clipper::MMonomer&monomer);
    static Matrix<float> calculate_plane_equation(const clipper::MMonomer& monomer);
};


#endif //BASEVALIDATION_H
