//
// Created by Jordan Dialpuri on 27/12/2023.
//

#ifndef VALIDATE_DATA_H
#define VALIDATE_DATA_H


/**
 * \brief HydrogenBondingAtoms
 * contains a list of donor atom names and acceptor atom names
 * construtcted with two vectors of strings
 */
struct HydrogenBondingAtoms {
    HydrogenBondingAtoms(const std::vector<std::string>&donor,
                         const std::vector<std::string>&acceptor) {
        donor_atoms = donor, acceptor_atoms = acceptor;
    }

    std::vector<std::string> donor_atoms;
    std::vector<std::string> acceptor_atoms;
};



struct HydrogenPositionData {
    HydrogenPositionData(const std::string&atom,
                         const std::string&ref_atom,
                         const std::vector<int>&angles,
                         const double bond_length) {
        h_bond_atom = atom;
        reference_atom = ref_atom;
        atom_angles = angles;
        h_bond_length = bond_length;
    }

    std::string h_bond_atom;
    std::string reference_atom;
    std::vector<int> atom_angles;
    double h_bond_length;
};

struct HydrogenPositions {
    explicit HydrogenPositions(const std::vector<HydrogenPositionData>&hpos) { arr = hpos; }
    std::vector<HydrogenPositionData> arr;
};


struct HBondData {

    // Expected H Bond data
    HydrogenBondingAtoms adenine = HydrogenBondingAtoms({"N6"}, {"O4"});
    HydrogenBondingAtoms guanine = HydrogenBondingAtoms({"N2", "N1"}, {"O2", "N3"});
    HydrogenBondingAtoms cytosine = HydrogenBondingAtoms( {"N4"}, {"O6"});
    HydrogenBondingAtoms thymine_uracil = HydrogenBondingAtoms({"N3"}, {"N1"});

    std::unordered_map<std::string, HydrogenBondingAtoms> hydrogen_bonding_map = {
        {"A", adenine},
      {"DA", adenine},
      {"G", guanine},
      {"DG",guanine},

      {"C", cytosine },
      {"DC", cytosine},
      {"DT", thymine_uracil},
      {"U", thymine_uracil},
    };

    // HydrogenPositionData
    HydrogenPositions adenine_positions = HydrogenPositions({
            {"N6", "C6", {120, -120}, 1.017}
        });

    HydrogenPositions guanine_positions = HydrogenPositions({
        {"N1", "C2", {-120}, 1.017},
        {"N2", "C2", {-120, 120}, 1.017}
    });

    HydrogenPositions thymine_positions = HydrogenPositions({
        {"N3", "C4", {120}, 1.017}
    });

    HydrogenPositions cytosine_positions = HydrogenPositions({
        {"N4", "C4", {120, -120}, 1.017}
    });

    std::map<std::string, HydrogenPositions> bonded_atoms = {
        {"A", adenine_positions},
        {"DA", adenine_positions},
        {"DT", thymine_positions},
        {"C", cytosine_positions},
        {"DC", cytosine_positions},
        {"G", guanine_positions},
        {"DG", guanine_positions},
        {"U", thymine_positions}
    };



};
#endif //VALIDATE_DATA_H
