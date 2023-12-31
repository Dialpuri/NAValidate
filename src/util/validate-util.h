//
// Created by Jordan Dialpuri on 18/12/2023.
//

#ifndef VALIDATE_UTIL_H
#define VALIDATE_UTIL_H

#include "clipper/clipper-minimol.h"
#include "fstream"
#include <map>
#include <set>
#include <string>
#include "../data/validate-data.h"


class ValidateUtil {
public:
    /** Math Functions */
    static float torsion(
        const clipper::Coord_orth&c1,
        const clipper::Coord_orth&c2,
        const clipper::Coord_orth&c3,
        const clipper::Coord_orth&c4) {
        return clipper::Util::rad2d(clipper::Coord_orth::torsion(c1, c2, c3, c4));
    }

    static float clip(float n, float lower, float upper) {
        return std::max(lower, std::min(n, upper));
    }

    /** Base Functiions*/
    static std::string base_type(const std::string&base_type_id);

    static HydrogenBondingAtoms get_h_atoms(const std::string&monomer_name);

    static clipper::MiniMol hydrogenate_model(clipper::MiniMol&mol);

    // Plane calculation
    static clipper::Vec3<> calculate_plane(const clipper::MMonomer&monomer);

    static std::vector<clipper::Vec3<>> calculate_vector_calculation_point(const clipper::MMonomer&monomer);

    static clipper::Mat33<> calculate_rodrigues_rotation_matrix(clipper::Vec3<> plane, float angle);

    // Hydrogenation
    static clipper::MMonomer hydrogenate_base(clipper::MMonomer&monomer);

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

#endif //VALIDATE_UTIL_H
