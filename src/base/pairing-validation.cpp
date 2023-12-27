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
    clipper::MAtomNonBond ns = clipper::MAtomNonBond(m_mol, 5);

    for (int p = 0; p < m_mol.size(); p++) {
        for (int m = 0; m < m_mol[p].size(); m++) {
            validate_basepair(m_mol[p][m], ns);
        }
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

void BasePairValidation::validate_basepair(clipper::MMonomer& mon, clipper::MAtomNonBond ns) {
    HydrogenBondingAtoms h_atoms = ValidateUtil::get_h_atoms(mon.type());

    std::cout << mon.type() << std::endl;
    // Need to find the monomers which could possibly form a H bond to the H bond atoms.
    clipper::MPolymer probe_chain;
    probe_chain.set_id("P");
    clipper::MMonomer probes;
    probes.set_id(1);
    probes.set_seqnum(1);
    probes.set_type("X");

    int l = 0;
    for (auto const& donor_atom: h_atoms.donor_atoms) {
        std::cout << donor_atom << std::endl;
        std::vector<std::string> bonded_h_atoms = find_bonded_h_atoms(mon, donor_atom);
        if (bonded_h_atoms.empty()) { std::cout << "No H atoms found on expected H bond donor, check hydrogenation\n";continue;}

        for (auto const& h_atom: bonded_h_atoms) {
            clipper::Vec3<> probe_position = get_probe_positon(mon, donor_atom, h_atom,2);
            if (probe_position.is_null()) continue;
            clipper::MAtom m;
            m.set_id(l);
            m.set_name("X");
            m.set_element("X");
            m.set_coord_orth(clipper::Coord_orth(probe_position));
            m.set_u_iso(0);

            l += 1;
            probes.insert(m);

            auto atoms_nearby = ns(clipper::Coord_orth(probe_position), 1);
            for (auto const& atom_near: atoms_nearby) {
                std::string residue_type = m_mol[atom_near.polymer()][atom_near.monomer()].type();
                BasePairResult result = is_base_pair(mon.type(), residue_type);
                if (result != pair) {continue;}

                std::cout << mon.type() << "-" << residue_type << " link found" << std::endl;
            }

        }
    }

    probe_chain.insert(probes);
    clipper::MiniMol xm = m_mol;
    xm.insert(probe_chain);
    clipper::MMDBfile mfile;
    mfile.export_minimol(xm);
    mfile.write_file("debug/1bna_probe_points.pdb");


}

