//
// Created by Jordan Dialpuri on 24/12/2023.
//

#include "pairing-validation.h"

#include "../util/validate-util.h"

BasePairValidation::BasePairValidation(clipper::MiniMol& mol) {

    // Hydrogenate the bases.
    m_mol = ValidateUtil::hydrogenate_model(mol);

    // clipper::MMDBfile mf;
    // mf.export_minimol(m_mol);
    // mf.write_file("debug/1bna_hydrogenated.pdb");
}

void BasePairValidation::validate() {
    clipper::MAtomNonBond ns = clipper::MAtomNonBond(m_mol, 5);

    std::map<int, HBondResult> found_bonds;

    clipper::MiniMol probe_points = {m_mol.spacegroup(), m_mol.cell()};

    for (int p = 0; p < m_mol.size(); p++) {
        clipper::MPolymer mp;
        mp.set_id(m_mol[p].id());
        for (int m = 0; m < m_mol[p].size(); m++) {

            clipper::MMonomer probe_monomer = add_debug_probes(m_mol[p][m]);
            mp.insert(probe_monomer);
            int monomer_id = std::stoi(m_mol[p][m].id().trim());
            HBondResult hbonds = calculate_h_bonds(m_mol[p][m], ns);
            found_bonds.insert({monomer_id, hbonds});

            // std::cout << monomer_id << " " << hbonds.paired_monomer_id << std::endl;
        }
        probe_points.insert(mp);
    }
    //
    // clipper::MMDBfile mfile;
    // mfile.export_minimol(probe_points);
    // mfile.write_file("debug/1bna_probe_points.pdb");

    std::map<int, int> found_pairs;
    std::set<int> checked;
    int at = 2;
    int cg = 3;

    for (auto const& bond: found_bonds) {
        if (checked.find(bond.first) != checked.end()) continue;

        int paired_id = bond.second.paired_monomer_id;
        std::string paired_type = bond.second.paired_monomer_type;

        if (found_bonds.find(paired_id) == found_bonds.end()) continue;
        auto it = found_bonds.find(paired_id);

        if (it->second.paired_monomer_id == bond.first) {
            int count = 0;
            count += bond.second.h_bonds;
            count += it->second.h_bonds;
            if (paired_type == "A" || paired_type == "T" ||
                paired_type == "DA" || paired_type == "DT" ||
                paired_type == "U" ) {
                if (count == at) {
                    found_pairs.insert({bond.first, bond.second.paired_monomer_id});
                }
            } else {
                if (count == cg) {
                    found_pairs.insert({bond.first, bond.second.paired_monomer_id});
                }
            }
        }
        checked.insert(bond.first);
        checked.insert(bond.second.paired_monomer_id);

    }

    for (auto const& pair: found_pairs) {
        std::cout << pair.first << "-" << pair.second << std::endl;
    }


}

std::vector<std::string> BasePairValidation::find_bonded_h_atoms(const clipper::MMonomer&mon,
                                                                 const std::string& donor_atom) {
    std::vector<std::string> bonded_h_atoms = {};
    std::string donor_h_atom = donor_atom;
    donor_h_atom.replace(0,1, "H");

    for (int i = 1; i <= 2; i++) {
        std::string potential_h_atom = donor_h_atom + std::to_string(i);
        int h_idx = mon.lookup(potential_h_atom, clipper::MM::UNIQUE);
        if (h_idx >= 0) {bonded_h_atoms.push_back(potential_h_atom);}
    }

    return bonded_h_atoms;
}

clipper::Vec3<> BasePairValidation::get_probe_positon(clipper::MMonomer& mon, std::string const& donor_atom,
                                                      const std::string& h_atom, double probe_length) noexcept {
    const int donor_atom_idx = mon.lookup(donor_atom, clipper::MM::UNIQUE);
    if (donor_atom_idx < 0) return clipper::Vec3<>::null();

    const int h_atom_idx = mon.lookup(h_atom, clipper::MM::UNIQUE);
    if (h_atom_idx < 0) return clipper::Vec3<>::null();

    const clipper::Vec3<> bond_vector = mon[h_atom_idx].coord_orth() - mon[donor_atom_idx].coord_orth();
    const clipper::Vec3<> probe_vector = bond_vector.unit() * probe_length;
    const clipper::Vec3<> probe_position = mon[donor_atom_idx].coord_orth() + probe_vector;
    return probe_position;
}

HBondResult BasePairValidation::calculate_h_bonds(clipper::MMonomer&mon, clipper::MAtomNonBond ns) {
    HydrogenBondingAtoms h_atoms = ValidateUtil::get_h_atoms(mon.type());

    int h_bond_count = 0;
    int paired_monomer_id = -1;
    std::string paired_monomer_type;

    for (auto const& donor_atom: h_atoms.donor_atoms) {
        std::cout << mon.type() << " " << mon.id().trim() << " " <<  donor_atom << std::endl;
        std::vector<std::string> bonded_h_atoms = find_bonded_h_atoms(mon, donor_atom);
        if (bonded_h_atoms.empty()) { std::cout << "No H atoms found on expected H bond donor, check hydrogenation\n";continue;}

        for (auto const& h_atom: bonded_h_atoms) {
            clipper::Vec3<> probe_position = get_probe_positon(mon, donor_atom, h_atom,2);
            if (probe_position.is_null()) continue;

            auto atoms_nearby = ns(clipper::Coord_orth(probe_position), 2);
            for (auto const& atom_near: atoms_nearby) {
                clipper::MMonomer residue_near = m_mol[atom_near.polymer()][atom_near.monomer()];
                std::string residue_type = residue_near.type();
                BasePairResult result = is_base_pair(mon.type(), residue_type);
                if (result != pair) {continue;}

                clipper::MAtom atom = residue_near[atom_near.atom()];

                if (std::find(h_atoms.acceptor_atoms.begin(), h_atoms.acceptor_atoms.end(),
                    atom.name().trim()) == h_atoms.acceptor_atoms.end()) {
                    continue;
                }

                h_bond_count++;
                paired_monomer_id = std::stoi(residue_near.id().trim());
                paired_monomer_type = residue_near.type().trim();
            }
        }
    }
    return {h_bond_count, paired_monomer_id, paired_monomer_type};
}

clipper::MMonomer BasePairValidation::add_debug_probes(clipper::MMonomer& mon) {
    HydrogenBondingAtoms h_atoms = ValidateUtil::get_h_atoms(mon.type());

    clipper::MMonomer mon_ = mon;

    for (auto const& donor_atom: h_atoms.donor_atoms) {
        std::vector<std::string> bonded_h_atoms = find_bonded_h_atoms(mon, donor_atom);
        if (bonded_h_atoms.empty()) { std::cout << "No H atoms found on expected H bond donor, check hydrogenation\n";continue;}

        for (auto const& h_atom: bonded_h_atoms) {
            clipper::Vec3<> probe_position = get_probe_positon(mon, donor_atom, h_atom,2);
            if (probe_position.is_null()) continue;

            clipper::MAtom ma;
            ma.set_coord_orth(clipper::Coord_orth(probe_position));
            ma.set_id("X");
            mon_.insert(ma);
        }
    }
    return mon_;
}
