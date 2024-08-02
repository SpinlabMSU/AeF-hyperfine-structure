#include "pch.h"
#include <aef/operators/PerturbationFramework.h>

using namespace aef;
using namespace aef::operators;

using basis_vec = j_basis_vec;
using system = HyperfineCalculator;

aef::operators::PerturbationFramework::PerturbationFramework(system *sys_):
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

system* aef::operators::PerturbationFramework::get_system() {
    return sys;
}

void aef::operators::PerturbationFramework::set_basis_system(system* calc) {
    if (calc)
        sys = calc;
}

IOperator<basis_vec>* aef::operators::PerturbationFramework::getOperator(std::string& id) {
    return &(opMap[id]);
}

void aef::operators::PerturbationFramework::addOperator(IOperator<basis_vec>* op, std::string& id) {
    opMap.emplace(id, op);
}

void aef::operators::PerturbationFramework::evaluate(void) {
    const Eigen::Index nBasisElts = sys->nBasisElts;
    for (auto &[id, op] : opMap) {
        if (!opMatMap.contains(id)) {
            Eigen::MatrixXcd* opMat = new Eigen::MatrixXcd(nBasisElts, nBasisElts);
        }
    }
}

Eigen::MatrixXcd* aef::operators::PerturbationFramework::getOperatorMatrix(std::string& id) {
    return nullptr;
}

dcomplex aef::operators::PerturbationFramework::get_matrix_element(std::string& id, int eidx1, int eidx2) {
    return dcomplex();
}

dcomplex aef::operators::PerturbationFramework::expectation_value(std::string& id, int eidx1) {
    return dcomplex();
}

aef::ResultCode aef::operators::PerturbationFramework::delta_E_lo(std::string& id, Eigen::VectorXcd& output, Eigen::MatrixXcd *workspace) {
    bool internal_workspace = !workspace;
    if (internal_workspace) {
        workspace = new Eigen::MatrixXcd;
        workspace->resize(sys->nBasisElts, sys->nBasisElts);
    }
    Eigen::MatrixXcd* op = getOperatorMatrix(id);
    // get expectation values
    aef::matrix::group_action(*workspace, sys->Vs, *op);
    output = workspace->diagonal();
    // second order ???

    if (internal_workspace) {
        delete workspace;
    }
}
