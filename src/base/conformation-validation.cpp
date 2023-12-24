//
// Created by Jordan Dialpuri on 18/12/2023.
//

#include "conformation-validation.h"

double BaseConformationValidation::calculate_chi(const clipper::MMonomer&monomer) {
    const clipper::Vec3<> plane = ValidateUtil::calculate_plane(monomer);
    const clipper::Vec3<> plane_u = plane.unit();

    const int i_c1 = monomer.lookup("C1'", clipper::MM::UNIQUE);
    const int i_o4 = monomer.lookup("O4'", clipper::MM::UNIQUE);
    if (i_c1 < 0 || i_o4 < 0) {
       return -0.0f;
    }

    const clipper::Coord_orth c1 = monomer[i_c1].coord_orth();
    const clipper::Coord_orth o4 = monomer[i_o4].coord_orth();

    const clipper::Coord_orth co = o4-c1;
    const clipper::Vec3<> co_vector = co;
    const clipper::Vec3<> co_vector_unit = co_vector.unit();

    const float dot = clipper::Vec3<>::dot(co_vector_unit, plane_u);

    const double angle = acos(ValidateUtil::clip(dot, -1, 1)) + M_PI/2;
    return clipper::Util::rad2d(angle);
}

Matrix<float> BaseConformationValidation::calculate_plane_equation(const clipper::MMonomer& monomer) {
    float xx_sum = 0.0;
    float xy_sum = 0.0;
    float x_sum = 0.0;
    float y_sum = 0.0;
    float yy_sum = 0.0;
    float xz_sum = 0.0;
    float yz_sum = 0.0;
    float z_sum = 0.0;

    for (int atom_index = 0; atom_index < monomer.size(); atom_index++) {

        clipper::Coord_orth atom_coord_orth = monomer[atom_index].coord_orth();

        xx_sum += pow(atom_coord_orth.x(),2);
        yy_sum += pow(atom_coord_orth.y(),2);
        xy_sum += atom_coord_orth.x() + atom_coord_orth.y();
        x_sum += atom_coord_orth.x();
        y_sum += atom_coord_orth.y();
        xz_sum += atom_coord_orth.x() + atom_coord_orth.z();
        yz_sum += atom_coord_orth.y() + atom_coord_orth.z();
        z_sum += atom_coord_orth.z();
    }

    Matrix<float> A(3, 3);

    A.set(0,0,xx_sum);
    A.set(0,1,xy_sum);
    A.set(0,2,x_sum);

    A.set(1,0, xy_sum);
    A.set(1,1, yy_sum);
    A.set(1,2, y_sum);

    A.set(2,0, x_sum);
    A.set(2,1, y_sum);
    A.set(2,2, monomer.size());

    Matrix<float> B(3, 1);

    B.set(0,0,xz_sum);
    B.set(0,1,yz_sum);
    B.set(0,2,z_sum);

    Matrix<float> A_inv = A.inverse();
    Matrix<float> X = A_inv.dot(B);

    // X.print();
    // X.print_eqn();

    return X;
}

