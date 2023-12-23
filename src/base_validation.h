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
    static double calculate_chi(const clipper::MMonomer&monomer);
    static Matrix<float> calculate_plane_equation(const clipper::MMonomer& monomer);
    static clipper::Vec3<> calculate_plane(const clipper::MMonomer& monomer);
    static std::vector<clipper::Vec3<>> calculate_vector_calculation_point(const clipper::MMonomer& monomer);

    static float clip(float n, float lower, float upper) {
        return std::max(lower, std::min(n, upper));
    }
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
