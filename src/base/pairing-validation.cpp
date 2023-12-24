//
// Created by Jordan Dialpuri on 24/12/2023.
//

#include "pairing-validation.h"

#include "../util/validate-util.h"

BasePairValidation::BasePairValidation(clipper::MiniMol& mol) {

    // Hydrogenate the bases.
    m_mol = ValidateUtil::hydrogenate_model(mol);

    clipper::MMDBfile mf;
    mf.export_minimol(m_mol);
    mf.write_file("debug/1bna_hydrogenated.pdb");
}

void BasePairValidation::validate() {
    for (int p = 0; p < m_mol.size(); p++) {
        for (int m = 0; m < m_mol[p].size(); m++) {
            validate_basepair(m_mol[p][m]);
        }
    }
}

void BasePairValidation::validate_basepair(clipper::MMonomer& mon) {
    HydrogenBondingAtoms h_atoms = ValidateUtil::get_h_atoms(mon.type());



}

