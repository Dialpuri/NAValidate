//
// Created by Jordan Dialpuri on 18/12/2023.
//

#ifndef VALIDATE_H
#define VALIDATE_H

#include <set>
#include <clipper/clipper-minimol.h>
#include "validate-util.h"
#include "pucker_validate.h"

class Validate {
public:
    explicit Validate(clipper::MiniMol& mol);

    void validate();

protected:
    clipper::MiniMol m_mol;
    std::set<std::string> m_na_names = {"A", "C", "G", "U", "DA", "DC", "DG", "DT"};

};


#endif //VALIDATE_H
