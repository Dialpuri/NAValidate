//
// Created by Jordan Dialpuri on 24/12/2023.
//

#include "validate-util.h"

std::string ValidateUtil::base_type(const std::string& base_type_id) {
    if (base_type_id == "C" || base_type_id == "T" || base_type_id == "U") {
        return "pyramidine";
    }
    if (base_type_id == "A" || base_type_id == "G") {
        return "purine";
    }
    return "UK";
}

HydrogenBondingAtoms ValidateUtil::get_h_atoms(const std::string& monomer_name) {
    HBondData data;
    const auto it = data.hydrogen_bonding_map.find(monomer_name);
    return it->second;
}

clipper::MiniMol ValidateUtil::hydrogenate_model(clipper::MiniMol& mol) {
    clipper::MiniMol h_mol = {mol.spacegroup(), mol.cell()};
    for (int p = 0; p < mol.size(); p++) {
        clipper::MPolymer mp;
        mp.set_id(mol[p].id());
        for (int m = 0; m < mol[p].size(); m++) {
            clipper::MMonomer h_mon = hydrogenate_base(mol[p][m]);
            mp.insert(h_mon);
        }
        h_mol.insert(mp);
    }
    return h_mol;
}


clipper::MMonomer ValidateUtil::hydrogenate_base(clipper::MMonomer& monomer) {
    HBondData data;
    clipper::MMonomer mon = monomer;

    const clipper::Vec3<> base_plane = calculate_plane(monomer);
    // const HydrogenBondingAtoms hydrogen_bonding_atoms = get_h_atoms(monomer.type());

    const auto h_bond_it = data.bonded_atoms.find(monomer.type().trim());
    if (h_bond_it == data.bonded_atoms.end()) {std::cout << "[Critical] Non recognised monomer! " << std::endl; return {};}

    for (const HydrogenPositionData& h_pos_data: h_bond_it->second.arr) {
        const std::string bonded_atom_name = h_pos_data.referance_atom;

        int no_h_atom = 1;

        for (const int& angle: h_pos_data.atom_angles) {
            const clipper::Mat33<> rot = calculate_rodrigues_rotation_matrix(base_plane.unit(), angle);
            const clipper::RTop_orth rotation = {rot,{0,0,0}};

            const int bonded_id = monomer.lookup(bonded_atom_name, clipper::MM::UNIQUE);
            if (bonded_id < 0) continue;

            const int current_id = monomer.lookup(h_pos_data.h_bond_atom, clipper::MM::UNIQUE);
            if (current_id < 0) continue;

            clipper::MAtom bonded_matom = monomer[bonded_id];
            const clipper::MAtom current_matom = monomer[current_id];

            const clipper::RTop_orth translate_to_origin = {clipper::Mat33<>::identity(), -current_matom.coord_orth()};
            const clipper::RTop_orth translate_from_origin = {clipper::Mat33<>::identity(), current_matom.coord_orth()};

            bonded_matom.transform(translate_to_origin);
            bonded_matom.transform(rotation);
            bonded_matom.transform(translate_from_origin);

            clipper::Vec3<> bond_vec = bonded_matom.coord_orth()-current_matom.coord_orth();
            clipper::Vec3<> bond_vec_unit = bond_vec.unit();
            clipper::Vec3<> hydrogen_vec = h_pos_data.h_bond_length*bond_vec_unit;
            bonded_matom.set_coord_orth(current_matom.coord_orth()+clipper::Coord_orth(hydrogen_vec));

            bonded_matom.set_id(h_pos_data.h_bond_atom);
            bonded_matom.set_element("H");

            std::string h_name = h_pos_data.h_bond_atom;
            h_name.replace(0,1, "H");
            h_name += std::to_string(no_h_atom);
            no_h_atom++;

            bonded_matom.set_name(h_name);
            mon.insert(bonded_matom);
        }
    }
    return mon;
}

clipper::Vec3<> ValidateUtil::calculate_plane(const clipper::MMonomer& monomer) {
    const std::vector<clipper::Vec3<>> points = calculate_vector_calculation_point(monomer);
    if (points.empty()) return {};

    const clipper::Vec3<> A = points[0];
    const clipper::Vec3<> B = points[1];
    const clipper::Vec3<> C = points[2];

    const clipper::Vec3<> BA = B-A;
    const clipper::Vec3<> CA = C-A;

    const clipper::Vec3<> n = clipper::Vec3<>::cross(BA, CA);
    return n;
}

std::vector<clipper::Vec3<>> ValidateUtil::calculate_vector_calculation_point(const clipper::MMonomer& monomer) {
        const std::set<std::string> purines = {"A", "G", "DG", "DA"};
        const std::set<std::string> pyrmidines = {"DT", "DC", "U", "C"};

        if (purines.find(monomer.type()) != purines.end()) {
            const int i_n9 = monomer.lookup("N9", clipper::MM::UNIQUE);
            const int i_n7 = monomer.lookup("N7", clipper::MM::UNIQUE);
            const int i_n1 = monomer.lookup("N1", clipper::MM::UNIQUE);

            if (i_n9 < 0 || i_n7 < 0 || i_n1 < 0) {
                return {};
            }

            return {
                monomer[i_n9].coord_orth(),
                monomer[i_n7].coord_orth(),
                monomer[i_n1].coord_orth()
            };
        }

        if (pyrmidines.find(monomer.type()) != pyrmidines.end()) {
            const int i_n1 = monomer.lookup("N1", clipper::MM::UNIQUE);
            const int i_c5 = monomer.lookup("C5", clipper::MM::UNIQUE);
            const int i_o2 = monomer.lookup("O2", clipper::MM::UNIQUE);

            if (i_n1 < 0 || i_c5 < 0 || i_o2 < 0) {
                return {};
            }

            return {
                monomer[i_n1].coord_orth(),
                monomer[i_c5].coord_orth(),
                monomer[i_o2].coord_orth()
            };
        }

        return {};

}

clipper::Mat33<> ValidateUtil::calculate_rodrigues_rotation_matrix(clipper::Vec3<> u, float theta) {
    const double r = theta * (M_PI / 180);
    const double x = u[0];
    const double y = u[1];
    const double z = u[2];
    const double c = cos(r);
    const double one_min_c = 1 - c;
    const double s = sin(r);

    const clipper::Mat33<> rot_mat = clipper::Mat33<>(c + (one_min_c * pow(x, 2)), (one_min_c * x * y) - (z * s),
                                                      (one_min_c * x * z) + (s * y),
                                                      (one_min_c * x * y) + (z * s), c + (one_min_c * pow(y, 2)),
                                                      (one_min_c * y * z) - (x * s),
                                                      (one_min_c * x * z) - (y * s), (one_min_c * y * z) + (x * s),
                                                      c + (one_min_c * pow(z, 2)));
    return rot_mat;
}
