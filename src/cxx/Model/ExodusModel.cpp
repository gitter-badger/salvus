//
// Created by Michael Afanasiev on 2016-02-01.
//

#include <openmpi/ompi/mpi/cxx/mpicxx.h>
#include "exodusII.h"
#include "ExodusModel.h"
#include "../Utilities.h"

ExodusModel::ExodusModel(Options options) {

    mExodusFileName = options.ExodusModelFile();

}

void ExodusModel::initializeParallel() {

    if (!MPI::COMM_WORLD.Get_rank()) {
        getInitialization();
        readCoordinates();
    }

    mNodalX = utilities::broadcastStdVecFromRoot(mNodalX);
    mNodalY = utilities::broadcastStdVecFromRoot(mNodalY);
    if (mNumberDimension > 2) { mNodalZ = utilities::broadcastStdVecFromRoot(mNodalZ); }

    createKdTree();
//    readNodalVariables();

}

void ExodusModel::exodusError(const int retval, std::string func_name) {

    if (retval) {
        throw std::runtime_error("Error in exodus function: " + func_name);
    }

}



void ExodusModel::readCoordinates() {

    try {
        mNodalX.resize(mNumberVertices);
        mNodalY.resize(mNumberVertices);
        if (mNumberDimension > 2) { mNodalZ.resize(mNumberVertices); }
        if (mNumberDimension == 2) {
            exodusError(ex_get_coord(
                                mExodusId, mNodalX.data(), mNodalY.data(), NULL),
                        "ex_get_coord");
        } else {
            exodusError(ex_get_coord(
                                mExodusId, mNodalX.data(), mNodalY.data(), mNodalZ.data()),
                        "ex_get_coord");
        }
    } catch (std::exception &e) {
        utilities::print_from_root_mpi(e.what());
        MPI_Abort(PETSC_COMM_WORLD, -1);
    }

}

void ExodusModel::getInitialization() {

    int io_ws = 0;
    int comp_ws = 8;
    try {
        mExodusId = ex_open(mExodusFileName.c_str(), EX_READ, &comp_ws, &io_ws, &mExodusVersion);
        if (mExodusId < 0) { throw std::runtime_error("Error opening exodus model file."); }

        exodusError(ex_get_init(
                            mExodusId, mExodusTitle, &mNumberDimension, &mNumberVertices, &mNumberElements,
                            &mNumberElementBlocks, &mNumberNodeSets, &mNumberSideSets),
                    "ex_get_init");

    } catch (std::exception &e) {
        utilities::print_from_root_mpi(e.what());
        MPI_Abort(PETSC_COMM_WORLD, -1);
    }

}

void ExodusModel::createKdTree() {

    // Need to keep the index arrays allocated.
    mKdTreeData.resize(mNodalX.size());

    // Create.
    if (mNumberDimension == 2) {
        mKdTree = kd_create(2);
        for (auto i = 0; i < mNumberVertices; i++) {
            mKdTreeData[i] = i;
            kd_insert(mKdTree, std::vector<double> {mNodalX[i], mNodalY[i]}.data(), &mKdTreeData[i]);
        }
    } else {
        mKdTree = kd_create(3);
        for (auto i = 0; i < mNumberVertices; i++) {
            mKdTreeData[i] = i;
            kd_insert(mKdTree, std::vector<double> {mNodalX[i], mNodalY[i], mNodalZ[i]}.data(), &mKdTreeData[i]);
        }
    }

}

void ExodusModel::readNodalVariables() {

    exodusError(ex_get_var_param(
            mExodusId, "n", &mNumberNodalVariables),
                "ex_get_var_param");
    char *nm[mNumberNodalVariables];
    for (auto i = 0; i < mNumberNodalVariables; i++) {
        mNodalVariableNames.push_back(std::string((char *) calloc((MAX_STR_LENGTH+1), sizeof(char))));
    }

}
