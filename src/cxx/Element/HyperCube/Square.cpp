//
// Created by Michael Afanasiev on 2016-01-30.
//

#include <iostream>
#include <petscdm.h>
#include <petscdmplex.h>
#include <openmpi/ompi/mpi/cxx/mpicxx.h>
#include "Square.h"
#include "Square/Autogen/order4_square.h"

/*
 * STATIC FUNCTIONS WHICH ARE ONLY ON THE REFERENCE ELEMENT.
 */
int Square::mNumberDofVertex;
int Square::mNumberDofEdge;
int Square::mNumberDofFace;
int Square::mNumberDofVolume;
int Square::mNumberIntegrationPointsEps;
int Square::mNumberIntegrationPointsEta;
int Square::mNumberIntegrationPoints;
int Square::mPolynomialOrder;

Eigen::VectorXi Square::mClosureMapping;
Eigen::MatrixXd Square::mGradientOperator;
Eigen::VectorXd Square::mIntegrationWeightsEps;
Eigen::VectorXd Square::mIntegrationWeightsEta;
Eigen::VectorXd Square::mIntegrationCoordinatesEps;
Eigen::VectorXd Square::mIntegrationCoordinatesEta;

double Square::n0(const double &eps, const double &eta) { return 0.25 * (1.0 - eps) * (1.0 - eta); }
double Square::n1(const double &eps, const double &eta) { return 0.25 * (1.0 + eps) * (1.0 - eta); }
double Square::n2(const double &eps, const double &eta) { return 0.25 * (1.0 - eps) * (1.0 + eta); }
double Square::n3(const double &eps, const double &eta) { return 0.25 * (1.0 + eps) * (1.0 + eta); }
double Square::dn0deps(const double &eta) { return (-1) * (1 - eta) / 4.0; }
double Square::dn1deps(const double &eta) { return (+1) * (1 - eta) / 4.0; }
double Square::dn2deps(const double &eta) { return (-1) * (1 + eta) / 4.0; }
double Square::dn3deps(const double &eta) { return (+1) * (1 + eta) / 4.0; }
double Square::dn0deta(const double &eps) { return (1 - eps) * -1.0 / 4.0; }
double Square::dn1deta(const double &eps) { return (1 + eps) * -1.0 / 4.0; }
double Square::dn2deta(const double &eps) { return (1 - eps) * 1.0 / 4.0; }
double Square::dn3deta(const double &eps) { return (1 + eps) * 1.0 / 4.0; }

Eigen::VectorXd Square::GllPointsForOrder(const int order) {
    Eigen::VectorXd gll_points(order+1);
    if (order == 4) {
        gll_points << -1.0, -0.6546536707, 0.0, 0.6546536707, 1.0;
    }
    return gll_points;
}

Eigen::VectorXd Square::GllIntegrationWeightForOrder(const int order) {
    Eigen::VectorXd integration_weights(order+1);
    if (order == 4) {
        integration_weights << 0.1, 0.5444444444, 0.7111111111, 0.5444444444, 0.1;
    }
    return integration_weights;
}

Eigen::VectorXi Square::ClosureMapping(const int order, const int dimension) {
    Eigen::VectorXi closure_mapping((order+1)*(order+1));
    if (dimension == 2) {
        if (order == 4) {
            closure_mapping << 6, 7, 8, 11, 12, 13, 16, 17, 18, 1, 2, 3,
                    9, 14, 19, 23, 22, 21, 15, 10, 5, 0, 4, 24, 20;
//            closure_mapping << 8, 13, 18, 7, 12, 17, 6, 11, 16, 9, 14, 19, 23, 22, 21, 15, 10, 5, 1, 2, 3, 4, 24, 20, 0;
        }
    }
    return closure_mapping;
}


Eigen::Map<const Eigen::VectorXd> Square::epsVectorStride(const Eigen::VectorXd &function,
                                                    const int &eta_index) {
    return Eigen::Map<const Eigen::VectorXd> (
            function.data() + eta_index * mNumberIntegrationPointsEta,
            mNumberIntegrationPointsEps);
}

Eigen::Map<const Eigen::VectorXd, 0, Eigen::InnerStride<>> Square::etaVectorStride(const Eigen::VectorXd &function,
                                                                                   const int &eta_index) {
    return Eigen::Map<const Eigen::VectorXd, 0, Eigen::InnerStride<>> (
            function.data() + eta_index, mNumberIntegrationPointsEta,
            Eigen::InnerStride<> (mNumberIntegrationPointsEps));
}


Eigen::Vector4d Square::interpolateShapeFunctions(const double &eps, const double &eta) {
    Eigen::Vector4d coefficients;
    coefficients(0) = n0(eps, eta);
    coefficients(1) = n1(eps, eta);
    coefficients(2) = n2(eps, eta);
    coefficients(3) = n3(eps, eta);
    return coefficients;
}

void Square::attachVertexCoordinates(DM &distributed_mesh) {

    Vec coordinates_local;
    PetscInt coordinate_buffer_size;
    PetscSection coordinate_section;
    PetscReal *coordinates_buffer = NULL;

    DMGetCoordinatesLocal(distributed_mesh, &coordinates_local);
    DMGetCoordinateSection(distributed_mesh, &coordinate_section);
    DMPlexVecGetClosure(distributed_mesh, coordinate_section, coordinates_local, mElementNumber,
                        &coordinate_buffer_size, &coordinates_buffer);
    std::vector<PetscReal> coordinates_element(coordinates_buffer, coordinates_buffer+coordinate_buffer_size);
    DMPlexVecRestoreClosure(distributed_mesh, coordinate_section, coordinates_local, mElementNumber,
                            &coordinate_buffer_size, &coordinates_buffer);

    // Reorder to desired vertex ordering.
    std::vector<double> vertex_coordinates_ordered ;
    std::vector<PetscInt> mapping_to_reference_element {6, 7, 0, 1, 4, 5, 2, 3};
    for (int i = 0; i < mNumberVertex; i++) {
        mVertexCoordinates(0,i) = coordinates_element[mapping_to_reference_element[mNumberDimensions*i+0]];
        mVertexCoordinates(1,i) = coordinates_element[mapping_to_reference_element[mNumberDimensions*i+1]];
    }

}

Eigen::Matrix<double,2,2> Square::jacobianAtPoint(PetscReal eps, PetscReal eta) {
    Eigen::Matrix<double,2,4> jacobian_multiplier;
    jacobian_multiplier(0,0) = dn0deps(eta);
    jacobian_multiplier(0,1) = dn1deps(eta);
    jacobian_multiplier(0,2) = dn2deps(eta);
    jacobian_multiplier(0,3) = dn3deps(eta);
    jacobian_multiplier(1,0) = dn0deta(eps);
    jacobian_multiplier(1,1) = dn1deta(eps);
    jacobian_multiplier(1,2) = dn2deta(eps);
    jacobian_multiplier(1,3) = dn3deta(eps);
    return jacobian_multiplier * mVertexCoordinates.transpose();
}

Eigen::Vector4d Square::__interpolateMaterialProperties(ExodusModel &model, std::string parameter_name) {

    Eigen::Vector4d material_at_vertices(mNumberVertex);

    // THIS IS STUPID JUST BECAUSE DENSITY IS NOT IN THE EXODUS FILE YET.
    if (parameter_name == "density") {
        material_at_vertices.setConstant(3.0);
        return material_at_vertices;
    }

    for (auto i = 0; i < mNumberVertex; i++) {
        material_at_vertices(i) = model.getMaterialParameterAtPoint({mVertexCoordinates(0, i),
                                                                     mVertexCoordinates(1, i)},
                                                                    parameter_name);
    }
    return material_at_vertices;

}

void Square::attachSource(std::vector<Source*> sources) {

    for (auto &source: sources) {
        if (mCheckHull(source->PhysicalLocationX(), source->PhysicalLocationZ())) {
            Eigen::Vector2d reference_location = inverseCoordinateTransform(source->PhysicalLocationX(),
                                                                            source->PhysicalLocationZ(),
                                                                            0.0, 0.0);
            source->setReferenceLocationEps(reference_location(0));
            source->setReferenceLocationEta(reference_location(1));
            mSources.push_back(source);
        }
    }

}

bool Square::mCheckHull(double x, double z) {
    int n_neg = 0;
    int n_pos = 0;
    std::vector<int> edge_mapping {0, 1, 3, 2, 0};
    Eigen::Vector2d test_point; test_point << x, z;
    for (auto i = 0; i < mNumberVertex; i++) {
        Eigen::Vector2d p0 = mVertexCoordinates.col(edge_mapping[i+0]);
        Eigen::Vector2d p1 = mVertexCoordinates.col(edge_mapping[i+1]);
        Eigen::Vector2d v_seg = p1 - p0;
        Eigen::Vector2d p_seg = test_point - p0;
        double x_0 = v_seg(0) * p_seg(1) - v_seg(1) * p_seg(0);
        if (x_0 <= 0) {
            n_neg++;
        } else {
            n_pos++;
        }
    }
    return n_neg == mNumberVertex || n_pos == mNumberVertex;
}

Eigen::Vector2d Square::inverseCoordinateTransform(const double &x_real, const double &z_real,
                                                   double eps, double eta) {

    double tol = 1e-6;
    Eigen::Vector2d solution {eps, eta};
    while (true) {

        eps = solution(0);
        eta = solution(1);

        Eigen::Matrix2d jacobian;
        Eigen::Vector4d shape_functions {n0(eps, eta), n1(eps, eta), n2(eps, eta), n3(eps, eta)};
        Eigen::Vector4d dNdEps {dn0deps(eta), dn1deps(eta), dn2deps(eta), dn3deps(eta)};
        Eigen::Vector4d dNdEta {dn0deta(eps), dn1deta(eps), dn2deta(eps), dn3deta(eps)};
        Eigen::Vector2d objective_function {x_real - shape_functions.dot(mVertexCoordinates.row(0)),
                                            z_real - shape_functions.dot(mVertexCoordinates.row(1))};

        jacobian << -1 * dNdEps.dot(mVertexCoordinates.row(0)),
                    -1 * dNdEta.dot(mVertexCoordinates.row(0)),
                    -1 * dNdEps.dot(mVertexCoordinates.row(1)),
                    -1 * dNdEta.dot(mVertexCoordinates.row(1));

        if ((objective_function.array().abs() < tol).all()) {
            return solution;
        } else {
            solution -= (jacobian.inverse() * objective_function);
        }

    }

}

void Square::readOperators() {

    int i = 0;
    double eta = mIntegrationCoordinatesEta[0];
    mGradientOperator.resize(mNumberIntegrationPointsEta, mNumberIntegrationPointsEps);
    Eigen::MatrixXd test(mNumberIntegrationPointsEta, mNumberIntegrationPointsEps);
    for (auto i=0; i < mNumberIntegrationPointsEps; i++) {
        double eps = mIntegrationCoordinatesEps[i];
        interpolate_eps_derivative_order4_square(eps, eta, test.data());
        mGradientOperator.row(i) = test.col(0);
    }

}

Square::Square(Options options) {

    // Basic properties.
    mPolynomialOrder = options.PolynomialOrder();

    // Gll points.
    mNumberDofVolume = 0;
    mNumberDofVertex = 1;
    mNumberDofEdge = mPolynomialOrder - 1;
    mNumberDofFace = (mPolynomialOrder - 1) * (mPolynomialOrder - 1);

    // Integration points.
    mIntegrationCoordinatesEps = Square::GllPointsForOrder(options.PolynomialOrder());
    mIntegrationCoordinatesEta = Square::GllPointsForOrder(options.PolynomialOrder());
    mIntegrationWeightsEps = Square::GllIntegrationWeightForOrder(options.PolynomialOrder());
    mIntegrationWeightsEta = Square::GllIntegrationWeightForOrder(options.PolynomialOrder());
    mClosureMapping = Square::ClosureMapping(options.PolynomialOrder(),mNumberDimensions);

    // Save number of integration points.
    mNumberIntegrationPointsEps = mIntegrationCoordinatesEps.size();
    mNumberIntegrationPointsEta = mIntegrationCoordinatesEta.size();
    mNumberIntegrationPoints = mNumberIntegrationPointsEps * mNumberIntegrationPointsEta;

}

void Square::scatterMassMatrix(Mesh *mesh) {

    mesh->setFieldOnElement((int) AcousticFields::mass_matrix, mElementNumber, mClosureMapping,
                            mMassMatrix);

}
