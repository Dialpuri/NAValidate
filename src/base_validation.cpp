//
// Created by Jordan Dialpuri on 18/12/2023.
//

#include "base_validation.h"

float BaseConformationValidation::calculate_chi(const clipper::MMonomer& monomer) {

    Matrix<float> plane = calculate_plane_equation(monomer);
    // plane.print_eqn();

    int i_c1 = monomer.lookup("C1'", clipper::MM::UNIQUE);
    int i_o4 = monomer.lookup("O4'", clipper::MM::UNIQUE);
    if (i_c1 < 0 || i_o4 < 0) {
        return -0.0f;
    }

    const clipper::Coord_orth c1 = monomer[i_c1].coord_orth();
    const clipper::Coord_orth o4 = monomer[i_o4].coord_orth();

    const clipper::Coord_orth co = o4-c1;
    clipper::Vec3<> co_vector = co;
    // clipper::Vec3<> co_unit_vector = co_vector.unit();
    clipper::Vec3<> normal = {plane(0,0), plane(0,1), plane(0,2)};
    // clipper::Vec3<> normal_unit = normal.unit();

    float co_length = sqrt(co.lengthsq());
    float n_length = sqrt(pow(normal[0],2)+pow(normal[1],2)+pow(normal[2],2));

    float angle = acos(clipper::Vec3<>::dot(co_vector, normal)/(co_length*n_length));
    std::cout << clipper::Util::rad2d(angle) << std::endl;

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