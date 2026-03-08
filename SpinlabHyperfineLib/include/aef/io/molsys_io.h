#pragma once
struct ExternalFieldParameters {
    double E_z;
    double K;
};

enum class molsys_save_version : uint16_t {
    // 
    invalid = 0,
    //
    initial = 1,
    //
    fix_hfs_stk_swap = 2,
    fix_payload_chunk = 3,
    _maxn,
    //
    max = _maxn - 1
};

constexpr molsys_save_version MINIMUM_LOAD_VERSION = molsys_save_version::initial;
constexpr molsys_save_version MAXIMUM_LOAD_VERSION = molsys_save_version::max;
constexpr molsys_save_version CURRENT_SAVE_VERSION = MAXIMUM_LOAD_VERSION;