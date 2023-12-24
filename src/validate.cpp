//
// Created by Jordan Dialpuri on 18/12/2023.
//

#include "validate.h"


Validate::Validate(clipper::MiniMol& mol) {
    clipper::MiniMol na_only_model;

    for (int p = 0; p < mol.size(); p++) {
        clipper::MPolymer mp;
        mp.set_id(mol[p].id());

        for (int m = 0; m < mol[p].size(); m++) {
            auto it = m_na_names.find(mol[p][m].type().trim());
            if (it != m_na_names.end()) {
                mp.insert(mol[p][m]);
            }
        }

        na_only_model.insert(mp);
    }
    m_mol = na_only_model;
}

void Validate::validate() {


    BasePairValidation bpv = BasePairValidation(m_mol);
    bpv.validate();
    return;

    std::vector<std::vector<std::string>> results;
    for (int p = 0; p < m_mol.size(); p++) {
        for (int m = 0; m < m_mol[p].size(); m++) {
            std::vector<std::string> result = {m_mol[p].id(), m_mol[p][m].id(), m_mol[p][m].type()};

            const PuckerResult pucker_result = PuckerValidate::calculate_pucker(m_mol[p][m]);
            result.emplace_back(std::to_string(pucker_result.P));
            result.emplace_back(std::to_string(pucker_result.theta_m));
            if (!pucker_result.is_null()) {
                const PuckerType pt = PuckerValidate::classify_pucker(pucker_result);
                result.push_back(pt.to_string());
            }
            result.push_back(ValidateUtil::base_type(m_mol[p][m].type()));

            double chi = BaseConformationValidation::calculate_chi(m_mol[p][m]);
            result.push_back(std::to_string(chi));

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