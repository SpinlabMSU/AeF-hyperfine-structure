#include "pch.h"
#include "aef/aef_run.h"

namespace fs = std::filesystem;


aef::aef_run::aef_run(char* path_) {
    runpath = fs::absolute(path_);

    //
    auto rp = runpath;
    bool fail = false;
    while (!fs::is_regular_file(rp / "out.log")) {
        //if (rp.)
    }
    if (!fail) {
        runpath = rp;
    }
}

aef::aef_run::aef_run(std::filesystem::path runpath_)
    :runpath(runpath_){
}

aef::aef_run::~aef_run() {}

aef::run_params aef::aef_run::parse_params(bool force_reparse) {
    return run_params();
}

bool aef::aef_run::is_valid() {
    if (!valid_checked) {

    }
    return valid;
}

bool aef::aef_run::has_matrix() {
    return false;
}

bool aef::aef_run::has_coeffs() {
    return false;
}

std::string aef::aef_run::get_run_name() {
    return std::string();
}

std::vector<double> aef::aef_run::get_stark_shifts() {
    return std::vector<double>();
}

fs::path aef::aef_run::get_run_path() {
    return runpath;
}

fs::path aef::aef_run::get_log_path() {
    return runpath / "out.log";
}

fs::path aef::aef_run::get_matrix_path() {
    return runpath / "matrix.dat";
}

fs::path aef::aef_run::get_coeff_dir_path() {
    return runpath / "state_coeffs";
}

