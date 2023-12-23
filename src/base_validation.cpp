//
// Created by Jordan Dialpuri on 18/12/2023.
//

#include "base_validation.h"

double BaseConformationValidation::calculate_chi(const clipper::MMonomer&monomer) {
    const clipper::Vec3<> plane = calculate_plane(monomer);
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

    const double angle = acos(clip(dot, -1, 1)) + M_PI/2;
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

clipper::Vec3<> BaseConformationValidation::calculate_plane(const clipper::MMonomer& monomer) {
    std::vector<clipper::Vec3<>> points = calculate_vector_calculation_point(monomer);
    if (points.empty()) return {};

    clipper::Vec3<> A = points[0];
    clipper::Vec3<> B = points[1];
    clipper::Vec3<> C = points[2];

    clipper::Vec3<> BA = B-A;
    clipper::Vec3<> CA = C-A;

    clipper::Vec3<> n = clipper::Vec3<>::cross(BA, CA);

    float d = clipper::Vec3<>::dot(n, B);
    // std::cout << monomer.type() << monomer.id() << n.format() << d << std::endl;
    // std::cout << "a, b, c, d = " << n[0] << "," << n[1] << "," << n[2] << "," << d << std::endl;
    return n;
}

std::vector<clipper::Vec3<>> BaseConformationValidation::calculate_vector_calculation_point(
    const clipper::MMonomer& monomer) {

    std::set<std::string> purines = {"A", "G", "DG", "DA"};
    std::set<std::string> pyrmidines = {"DT", "DC", "U", "C"};

    if (purines.find(monomer.type()) != purines.end()) {
        int i_n9 = monomer.lookup("N9", clipper::MM::UNIQUE);
        int i_n7 = monomer.lookup("N7", clipper::MM::UNIQUE);
        int i_n1 = monomer.lookup("N1", clipper::MM::UNIQUE);

        if (i_n9 < 0 || i_n7 < 0 || i_n1 < 0) {
            return {};
        }

        return {
            monomer[i_n9].coord_orth(),
            monomer[i_n7].coord_orth(),
            monomer[i_n1].coord_orth()
        };
    }

    if (pyrmidines.find(monomer.type()) != pyrmidines.end()) {
        int i_n1 = monomer.lookup("N1", clipper::MM::UNIQUE);
        int i_c5 = monomer.lookup("C5", clipper::MM::UNIQUE);
        int i_o2 = monomer.lookup("O2", clipper::MM::UNIQUE);

        if (i_n1 < 0 || i_c5 < 0 || i_o2 < 0) {
            return {};
        }

        return {
            monomer[i_n1].coord_orth(),
            monomer[i_c5].coord_orth(),
            monomer[i_o2].coord_orth()
        };
    }

}
