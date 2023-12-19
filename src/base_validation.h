//
// Created by Jordan Dialpuri on 18/12/2023.
//

#ifndef BASEVALIDATION_H
#define BASEVALIDATION_H
#include <set>
#include <string>
#include "clipper/clipper-minimol.h"
#include "matrix.h"


struct BaseConformationResult {
    float chi;
    std::string conformation;
};


class BaseConformationValidation {
public:
    static float calculate_chi(const clipper::MMonomer& monomer);
    static Matrix<float> calculate_plane_equation(const clipper::MMonomer& monomer);

private:
    std::set<std::string> base_atoms = {
        "C1",
        "C2",
        "C3",
        "C4",
        "C5",
        "C6",
        "C7",
        "C8",
        "N1",
        "N2",
        "N3",
        "N4",
        "N5",
        "N6",
        "N7",
        "N8",
        "N9",
        "O2",
        "O4",
        "O6"
    };
};


#endif //BASEVALIDATION_H
