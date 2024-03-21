#include "cusolver_utils.h"
namespace cuda {
    // Returns cudaDataType value as defined in library_types.h for the string containing type name
    cudaDataType get_cuda_library_type(std::string type_string) {
        if (type_string.compare("CUDA_R_16F") == 0)
            return CUDA_R_16F;
        else if (type_string.compare("CUDA_C_16F") == 0)
            return CUDA_C_16F;
        else if (type_string.compare("CUDA_R_32F") == 0)
            return CUDA_R_32F;
        else if (type_string.compare("CUDA_C_32F") == 0)
            return CUDA_C_32F;
        else if (type_string.compare("CUDA_R_64F") == 0)
            return CUDA_R_64F;
        else if (type_string.compare("CUDA_C_64F") == 0)
            return CUDA_C_64F;
        else if (type_string.compare("CUDA_R_8I") == 0)
            return CUDA_R_8I;
        else if (type_string.compare("CUDA_C_8I") == 0)
            return CUDA_C_8I;
        else if (type_string.compare("CUDA_R_8U") == 0)
            return CUDA_R_8U;
        else if (type_string.compare("CUDA_C_8U") == 0)
            return CUDA_C_8U;
        else if (type_string.compare("CUDA_R_32I") == 0)
            return CUDA_R_32I;
        else if (type_string.compare("CUDA_C_32I") == 0)
            return CUDA_C_32I;
        else if (type_string.compare("CUDA_R_32U") == 0)
            return CUDA_R_32U;
        else if (type_string.compare("CUDA_C_32U") == 0)
            return CUDA_C_32U;
        else
            throw std::runtime_error("Unknown CUDA datatype");
    }

    cusolverIRSRefinement_t get_cusolver_refinement_solver(std::string solver_string) {
        if (solver_string.compare("CUSOLVER_IRS_REFINE_NONE") == 0)
            return CUSOLVER_IRS_REFINE_NONE;
        else if (solver_string.compare("CUSOLVER_IRS_REFINE_CLASSICAL") == 0)
            return CUSOLVER_IRS_REFINE_CLASSICAL;
        else if (solver_string.compare("CUSOLVER_IRS_REFINE_GMRES") == 0)
            return CUSOLVER_IRS_REFINE_GMRES;
        else if (solver_string.compare("CUSOLVER_IRS_REFINE_CLASSICAL_GMRES") == 0)
            return CUSOLVER_IRS_REFINE_CLASSICAL_GMRES;
        else if (solver_string.compare("CUSOLVER_IRS_REFINE_GMRES_GMRES") == 0)
            return CUSOLVER_IRS_REFINE_GMRES_GMRES;
        else
            printf("Unknown solver parameter: \"%s\"\n", solver_string.c_str());

        return CUSOLVER_IRS_REFINE_NOT_SET;
    }

};