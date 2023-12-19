//
// Created by Jordan Dialpuri on 17/12/2023.
//
#include "pucker_validate.h"

PuckerResult PuckerValidate::calculate_pucker(const clipper::MMonomer& monomer) {
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

    const clipper::Coord_orth c1 = monomer[i_c1].coord_orth();
    const clipper::Coord_orth c2 = monomer[i_c2].coord_orth();
    const clipper::Coord_orth c3 = monomer[i_c3].coord_orth();
    const clipper::Coord_orth c4 = monomer[i_c4].coord_orth();
    // clipper::Coord_orth c5 = monomer[i_c5].coord_orth();
    const clipper::Coord_orth o4 = monomer[i_o4].coord_orth();

    const float phi_0_d = ValidateUtil::torsion(c1,c2,c3,c4);
    const float phi_1_d = ValidateUtil::torsion(c2,c3,c4,o4);
    const float phi_2_d = ValidateUtil::torsion(c3,c4,o4,c1);
    const float phi_3_d = ValidateUtil::torsion(c4,o4,c1,c2);
    const float phi_4_d = ValidateUtil::torsion(o4,c1,c2,c3);

    const float numerator = (phi_2_d + phi_4_d) - (phi_1_d + phi_3_d);
    const float denominator =  phi_0_d * 3.077;
    const float eqn = numerator / denominator;

    float P = clipper::Util::rad2d(atan(eqn));
    const float theta_m = phi_0_d / clipper::Util::rad2d(cos(clipper::Util::d2rad(P)));

    if (clipper::Util::is_nan(P)) {
        return PuckerResult::null();
    }

    if (P < 0) {
        P += 360;
    }

    const PuckerResult pr = {P, theta_m};
    return pr;
}

PuckerType PuckerValidate::classify_pucker(const PuckerResult& pucker) {
    const std::vector<PuckerConformation> cyclic_pucker_types = {_3T2, __3E, _3T4, __E4, _OT4, __OE, _OT1, __E1, _2T1, __2E,
                                                        _2T3, __E3, _4T3, __4E, _4TO, __EO, _1TO, __1E, _1T2, __E2};
    int index = round(pucker.P/18);
    if (index >= cyclic_pucker_types.size()) {
        index = index-cyclic_pucker_types.size();
    }
    PuckerType pt;
    pt.conformation = cyclic_pucker_types[index+1];
    return pt;
}

std::string PuckerType::to_string() const {
    switch (conformation) {
        case _3T2: return "3T2";
        case __3E: return "3E";
        case _3T4: return "3T4";
        case __E4: return "E4";
        case _OT4: return "OT4";
        case __OE: return "OE";
        case _OT1: return "OT1";
        case __E1: return "E1";
        case _2T1: return "2T1";
        case __2E: return "2E";
        case _2T3: return "2T3";
        case __E3: return "E3";
        case _4T3: return "4T3";
        case __4E: return "4E";
        case _4TO: return "4TO";
        case __EO: return "EO";
        case _1TO: return "1TO";
        case __1E: return "1E";
        case _1T2: return "1T2";
        case __E2: return "E2";
        default: return "UK";
    }
}
