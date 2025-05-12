#include "pch.h"
#include "aef/aef_run.h"

namespace fs = std::filesystem;


aef::aef_run::aef_run(char* path_) {
    runpath = fs::absolute(path_);

    //
    fs::path rp = runpath;
    bool fail = false;
    while (!fs::is_regular_file(rp / "out.log")) {
        if (rp.parent_path() == rp) {
            fail = true;
        }
        rp = rp.parent_path();
    }
    if (!fail) {
        runpath = rp;
    } else {
        auto serr = fmt::format("Invalid aef run path {}", path_);
        throw std::runtime_error(serr);
    }
    stk_shifts.clear();
    durations.clear();
    run_name = (char*)(rp.filename().u8string().c_str()); 
}

aef::aef_run::aef_run(std::filesystem::path runpath_)
    :runpath(runpath_){
}

aef::aef_run::~aef_run() {}
namespace {
    std::string get_ssv_val(std::string line, std::string marker) {
        //
        auto val_idx = line.find(marker) + marker.length() + 1;
        auto next_space_idx = line.find(' ', val_idx);
        return line.substr(val_idx, next_space_idx);
    }

    std::string get_ssv_val(std::string line, std::string marker, int nExtraSpaces) {
        //
        auto val_idx = line.find(marker) + marker.length() + 1;
        auto next_space_idx = line.find(' ', val_idx);
        while (nExtraSpaces-- > 0) {
            next_space_idx = line.find(' ', val_idx);
        } 
        return line.substr(val_idx, next_space_idx);
    }
};
aef::run_params aef::aef_run::parse_params(bool force_reparse) {
    if (valid_params && !force_reparse) {
        // early out
        return params;
    }

    // parse parameter file
    std::ifstream in(get_log_path());

    // line-by-line parse
    std::string line;

    constexpr const char *n_search = "nmax is";
    constexpr const char *e_search = "E_z is";
    constexpr const char *k_search = "K is";
    constexpr const char *Emax_search = "Electric field strength is";
    constexpr const char *sta_search = "Start time is";
    constexpr const char *dur_search = "have taken";

    memset(&params, 0xff, sizeof(params));
    stk_shifts.clear();

    while(std::getline(in, line)) {
        // parse line
        if(line.contains(n_search) && line.contains(k_search)) {
            // this is the param line
            params.nmax = std::stoi(get_ssv_val(line, n_search));
            params.calc_E_z = std::stoi(get_ssv_val(line, e_search));
            params.K = std::stod(get_ssv_val(line, k_search));
        } else if (line.contains(Emax_search)) {
            // electric field strength
            double E_z = std::stod(get_ssv_val(line, Emax_search));
            params.max_E_z = std::max(E_z, params.max_E_z);
            stk_shifts.push_back(E_z);
        } else if (line.contains(sta_search)){
            // TODO implement
            std::string t = get_ssv_val(line, sta_search, 1);
        } else if (line.contains(dur_search)) {
            // time line
            // TODO implement
        }
    }

    if (params.nmax < 0 || !std::isfinite(params.calc_E_z) || !std::isfinite(params.K)) {
        // invalid
        valid_params = false;
        // what to do??
        constexpr const char * fmt_str = "log file {} valid but does not contain correct paramater line"; 
        std::string serr = std::format(fmt_str, get_log_path().generic_string());
        throw std::runtime_error(serr);
        aef::unreachable();
    }

    if (!std::isfinite(params.max_E_z)) {
        valid_params = false;
        constexpr const char *fmt_str = "log file {} valid and has correct parameter line but does not contain max E_z (no diagonalization output?)"; 
        std::string serr = std::format(fmt_str, get_log_path().generic_string());
        throw std::runtime_error(serr);
        aef::unreachable();
    }

    valid_params = true;
    return params;
}

bool aef::aef_run::is_valid(bool force_recheck) {
    if (valid_checked && !force_recheck) {
        return valid;
    }
    
    // check existence of the log, for right now this is the only check
    bool bLog = this->has_log();
    if(!bLog) {
        valid = false;
        goto end;
    }

    valid = true;
    end:
    valid_checked = true;
    return valid;
}

bool aef::aef_run::has_matrix() {
    return fs::is_regular_file(get_matrix_path());
}

bool aef::aef_run::has_log() {
    return fs::is_regular_file(get_log_path());
}

bool aef::aef_run::has_coeffs() {
    return fs::is_directory(get_matrix_path());
}

std::string aef::aef_run::get_run_name() {
    return run_name;
}

std::vector<double> aef::aef_run::get_stark_shifts() {
    return stk_shifts;
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

