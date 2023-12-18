//
// Created by Jordan Dialpuri on 17/12/2023.
//
#include "pucker_validate.h"

ValidateSugar::ValidateSugar(clipper::MiniMol& mol) {
    clipper::MiniMol na_only_model;

    for (int p = 0; p < mol.size(); p++) {
        clipper::MPolymer mp;
        mp.set_id(mol[p].id());

        for (int m = 0; m < mol[p].size(); m++) {
            auto it = m_NA_names.find(mol[p][m].type().trim());
            if (it != m_NA_names.end()) {
                mp.insert(mol[p][m]);
            }
        }

        na_only_model.insert(mp);
    }
    m_mol = na_only_model;
}

void ValidateSugar::validate() {
    std::vector<std::vector<std::string>> results;
    for (int p = 0; p < m_mol.size(); p++) {
        for (int m = 0; m < m_mol[p].size(); m++) {
            std::vector<std::string> result = {m_mol[p].id(), m_mol[p][m].id(), m_mol[p][m].type()};

            const PuckerResult pucker_result = calculate_pucker(m_mol[p][m]);
            result.emplace_back(std::to_string(pucker_result.P));
            result.emplace_back(std::to_string(pucker_result.theta_m));
            if (!pucker_result.is_null()) {
                const PuckerType pt = classify_pucker(pucker_result);
                result.push_back(pucker_type_to_string(pt));
            }
            result.push_back(type_classification(m_mol[p][m].type()));

            const BaseConformationResult conformation_result = calculate_base_conformation(m_mol[p][m]);
            result.emplace_back(std::to_string(pucker_result.P));
            result.emplace_back(std::to_string(pucker_result.theta_m));


            results.push_back(result);
        }
    }

    std::ofstream ofile;
    ofile.open("./results/1hr2.pdb");

    ofile << "ChainID,ResidueID,ResidueType,Pucker_P,Pucker_Theta,Pucker_Conformation,Base_Type,Base_Angle,Base_Conformation_Classification\n";
    for (const auto& r: results) {
        for (const auto& a: r) {
            ofile << a << ",";
        }
        ofile << "\n";
    }
    ofile.close();

}

PuckerResult ValidateSugar::calculate_pucker(clipper::MMonomer& monomer) {
// Torsions according to https://pubs.acs.org/doi/10.1021/acs.orglett.6b03626
    int i_c1 = monomer.lookup("C1'", clipper::MM::UNIQUE);
    int i_c2 = monomer.lookup("C2'", clipper::MM::UNIQUE);
    int i_c3 = monomer.lookup("C3'", clipper::MM::UNIQUE);
    int i_c4 = monomer.lookup("C4'", clipper::MM::UNIQUE);
    int i_c5 = monomer.lookup("C5'", clipper::MM::UNIQUE);
    int i_o4 = monomer.lookup("O4'", clipper::MM::UNIQUE);

    if (i_c1 < 0 || i_c2 < 0 || i_c3 < 0 || i_c4 < 0 || i_c5 < 0 || i_o4 < 0) {
        std::cout << "Invalid monomer configuration";
        return PuckerResult::null();
    }

    clipper::Coord_orth c1 = monomer[i_c1].coord_orth();
    clipper::Coord_orth c2 = monomer[i_c2].coord_orth();
    clipper::Coord_orth c3 = monomer[i_c3].coord_orth();
    clipper::Coord_orth c4 = monomer[i_c4].coord_orth();
    // clipper::Coord_orth c5 = monomer[i_c5].coord_orth();
    clipper::Coord_orth o4 = monomer[i_o4].coord_orth();

    float phi_0_d = Util::torsion(c1,c2,c3,c4);
    float phi_1_d = Util::torsion(c2,c3,c4,o4);
    float phi_2_d = Util::torsion(c3,c4,o4,c1);
    float phi_3_d = Util::torsion(c4,o4,c1,c2);
    float phi_4_d = Util::torsion(o4,c1,c2,c3);

    float numerator = (phi_2_d + phi_4_d) - (phi_1_d + phi_3_d);
    float denominator =  phi_0_d * 3.077;
    float eqn = numerator / denominator;

    float P = clipper::Util::rad2d(atan(eqn));
    float theta_m = phi_0_d / clipper::Util::rad2d(cos(clipper::Util::d2rad(P)));

    if (clipper::Util::is_nan(P)) {
        return PuckerResult::null();
    }

    if (P < 0) {
        P += 360;
    }

    PuckerResult pr;
    pr.P = P;
    pr.theta_m = theta_m;
    return pr;
}

PuckerType ValidateSugar::classify_pucker(PuckerResult pucker) {
    int index = round(pucker.P/18);
    if (index >= m_cyclic_pucker_types.size()) {
        index = index-m_cyclic_pucker_types.size();
    }
    return m_cyclic_pucker_types[index+1];
}

std::string PuckerType::to_string() const {

}
