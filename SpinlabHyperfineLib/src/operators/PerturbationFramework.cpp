#include "pch.h"
#include <aef/operators/PerturbationFramework.h>

using namespace aef;
using namespace aef::operators;

using basis_vec = j_basis_vec;
using System = HyperfineCalculator;

aef::operators::PerturbationFramework::PerturbationFramework(System *sys_):
    sys(sys_)
{
    opMap.clear();
    opMatMap.clear();
    pImpl = nullptr;
}

aef::operators::PerturbationFramework::~PerturbationFramework() {
    // deallocate operator matricies
    for (auto& [id, matrix] : opMatMap) {
        delete matrix;
    }
}

System* aef::operators::PerturbationFramework::get_system() {
    return sys;
}

void aef::operators::PerturbationFramework::set_basis_system(System* calc) {
    if (calc)
        sys = calc;
}

IKetOperator<basis_vec>* aef::operators::PerturbationFramework::getOperator(const std::string& id) {
    return opMap[id];
}

void aef::operators::PerturbationFramework::addOperator(const std::string& id, IKetOperator<basis_vec>* op) {
    opMap.emplace(id, op);
}

void aef::operators::PerturbationFramework::evaluate(void) {
    const Eigen::Index nBasisElts = sys->nBasisElts;
    for (auto &[id, op] : opMap) {
        Eigen::MatrixXcd* opMat;
        if (!opMatMap.contains(id)) {
            opMat = new Eigen::MatrixXcd(nBasisElts, nBasisElts);
            opMatMap.emplace(id, opMat);
        } else {
            opMat = opMatMap[id];
        }
        op->fillMatrix(*opMat);
    }
}

Eigen::MatrixXcd* aef::operators::PerturbationFramework::getOperatorMatrix(const std::string& id) {
    return opMatMap[id];
}

dcomplex aef::operators::PerturbationFramework::get_matrix_element(const std::string& id, int eidx1, int eidx2) {
    Eigen::VectorXcd state1 = sys->Vs.col(eidx1);
    Eigen::VectorXcd state2 = sys->Vs.col(eidx2);
    Eigen::MatrixXcd *op = getOperatorMatrix(id);

    dcomplex out;
    aef::matrix::matrix_element(out, state1, *op, state2);
    //return (state1.adjoint() * *op * state2);
    return out;
}

dcomplex aef::operators::PerturbationFramework::expectation_value(const std::string& id, int eidx1) {
    auto state1 = sys->Vs.col(eidx1);
    return (state1.adjoint() * *getOperatorMatrix(id) * state1);
}

aef::ResultCode aef::operators::PerturbationFramework::delta_E_lo(const std::string& id, Eigen::VectorXcd& output, Eigen::MatrixXcd *workspace) {
    bool internal_workspace = !workspace;
    if (internal_workspace) {
        workspace = new Eigen::MatrixXcd;
        workspace->resize(sys->nBasisElts, sys->nBasisElts);
    }
    Eigen::MatrixXcd* op = getOperatorMatrix(id);
    // get expectation values
    ResultCode rc = aef::matrix::group_action(*workspace, sys->Vs, *op);
    output = workspace->diagonal();
    // second order ???

    if (internal_workspace) {
        delete workspace;
    }
    return rc;
}
