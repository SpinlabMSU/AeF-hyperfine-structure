/*
    aef/aef_types.h -- contains the aef_run utility class

    This file is part of the AeF-hyperfine-structure program.

    AeF-hyperfine-structure is free software: you can redistribute it and/or
    modify it under the terms of the GNU General Public License as published
    by the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    AeF-hyperfine-structure is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
    or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for
    more details.

    You should have received a copy of the GNU General Public License along with
    AeF-hyperfine-structure. If not, see <https://www.gnu.org/licenses/>.
*/
#ifndef _AEF_AEF_RUN_H
#define _AEF_AEF_RUN_H 1
#pragma once
#include <aef/aef_types.h>
#include <filesystem>
#include <chrono>

namespace aef {
    namespace fs = std::filesystem;

    struct run_params {
        spin nmax;
        double calc_E_z;
        double max_E_z;
        double K;
    };

    class aef_run {
        fs::path runpath;
        bool valid_checked;
        bool valid;
        bool valid_params;
        run_params params;
        std::vector<double> stk_shifts;
        std::chrono::system_clock::time_point start_time;
        std::vector<std::chrono::nanoseconds> durations;
        std::string run_name;

    public:
        aef_run(char* path);
        aef_run(std::filesystem::path runpath);
        ~aef_run();

        run_params parse_params(bool force_reparse=false);
        
        bool is_valid(bool force_recheck=false);
        bool has_log();
        bool has_matrix();
        bool has_coeffs();

        std::string get_run_name();



        std::vector<double> get_stark_shifts();

        fs::path get_run_path();
        fs::path get_log_path();
        fs::path get_matrix_path();
        fs::path get_coeff_dir_path();

        

    };

};
#endif //_AEF_AEF_RUN_H