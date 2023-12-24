//
// Created by Jordan Dialpuri on 24/12/2023.

#ifndef BASE_PAIRING_VALIDATION_H
#define BASE_PAIRING_VALIDATION_H
#include "clipper/clipper-minimol.h"


/**
 * \brief BasePairValidation
 * constructed with a MiniMol object
 * call validate() to validate the entire model
 */
class BasePairValidation {
public:
    explicit BasePairValidation(clipper::MiniMol& mol);

    void validate();

    static void validate_basepair(clipper::MMonomer&mon);

private:
    clipper::MAtomNonBond m_ns;
    clipper::MiniMol m_mol;

};


#endif //BASE_PAIRING_VALIDATION_H
