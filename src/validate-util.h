//
// Created by Jordan Dialpuri on 18/12/2023.
//

#ifndef VALIDATE_UTIL_H
#define VALIDATE_UTIL_H

#include "clipper/clipper-minimol.h"
#include "fstream"

class ValidateUtil {
public:
    static float torsion(
        const clipper::Coord_orth& c1,
        const clipper::Coord_orth& c2,
        const clipper::Coord_orth& c3,
        const clipper::Coord_orth& c4) {

        return clipper::Util::rad2d(clipper::Coord_orth::torsion(c1,c2,c3,c4));
    }

    static std::string base_type(const std::string& base_type_id) {
        if (base_type_id == "C" || base_type_id == "T" || base_type_id == "U") {
            return "pyramidine";
        }
        if (base_type_id == "A" || base_type_id == "G") {
            return "purine";
        }
        return "UK";
    }

};

#endif //VALIDATE_UTIL_H
