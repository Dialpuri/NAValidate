//
// Created by Jordan Dialpuri on 18/12/2023.
//

#include "validate.h"

Validate::Validate(const clipper::MiniMol& mol) {
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