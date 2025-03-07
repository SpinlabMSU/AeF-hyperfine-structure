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
    max
};

constexpr molsys_save_version MINIMUM_LOAD_VERSION = molsys_save_version::initial;
constexpr molsys_save_version MAX_SUPPORTED_VERSION = molsys_save_version::initial;