//
// Created by Jordan Dialpuri on 24/12/2023.

#ifndef BASE_PAIRING_VALIDATION_H
#define BASE_PAIRING_VALIDATION_H

#include "clipper/clipper-minimol.h"
#include <map>
#include <set>

enum BasePairResult {
    pair, non_pair, not_recognised
};

/**
 * \brief BasePairValidation
 * constructed with a MiniMol object
 * call validate() to validate the entire model
 */
class BasePairValidation {
public:
    explicit BasePairValidation(clipper::MiniMol& mol);

    void validate();

    static std::vector<std::string> find_bonded_h_atoms(const clipper::MMonomer&mon, const std::string&donor_atom);

    static clipper::Vec3<> get_probe_positon(clipper::MMonomer& mon, const std::string &donor_atom,
        const std::string& h_aton, double probe_length=3) noexcept;

    void validate_basepair(clipper::MMonomer& mon, clipper::MAtomNonBond ns);

private:
    clipper::MAtomNonBond m_ns;
    clipper::MiniMol m_mol;

    BasePairResult is_base_pair(const std::string& a, const std::string& b) {
        const auto a_it = base_pairing.find(a);
        if (a_it == base_pairing.end()) {return not_recognised;}

        const auto b_it = a_it->second.find(b);
        if (b_it == a_it->second.end()) {return non_pair;}

        return pair;

    }

    std::map<std::string, std::set<std::string>> base_pairing {
        {"A", {"U", "DT"}},
        {"DA", {"U", "DT"}},

        {"DT", {"A", "DT"}},

        {"C", {"G", "DG"}},
        {"DC", {"G", "DG"}},

        {"G", {"C", "DC"}},
        {"DG", {"C", "DC"}}
    };

};


#endif //BASE_PAIRING_VALIDATION_H
