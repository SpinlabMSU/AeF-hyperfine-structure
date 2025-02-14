#include <pch.h>
#include <aef/root_matrix_io.h>

TTree *aef::write_self_adj_op_tree(Eigen::MatrixXcd& op, const char* name, const char* title, double mag_thresh, bool use_rel_thresh) {
    TTree *matrix_tree = new TTree(name, title);

    // threshold calculation
    const Eigen::Index nBasisElts = op.rows();

    double threshold = mag_thresh;

    if (use_rel_thresh) {
        double min_mag = std::numeric_limits<double>::infinity();
        double max_mag = -1;
        for (Eigen::Index jdx = 0; jdx < nBasisElts; jdx++) {
            for (Eigen::Index idx = 0; idx < nBasisElts; idx++) {
                aef::dcomplex val = op(idx, jdx);
                double mag = std::norm(val);
                min_mag = std::min(min_mag, mag);
                max_mag = std::max(max_mag, mag);
            }
        }
        threshold = mag_thresh * max_mag;
    }


    uint32_t sidx;
    Eigen::Index idx, jdx;
    aef::dcomplex val;
    double re, im, mag;

    size_t nFills = 0; // number of matrix elements large enough to be "filled" and written out
    size_t nZeros = 0; // number of matrix elements that are exactly zero
    size_t nFails = 0; // number of matrix elements that are 
    size_t nElts = nBasisElts * (nBasisElts + 1) / 2; // only need to write out half of the matrix

    matrix_tree->Branch("sidx", &sidx, "sidx/I");
    matrix_tree->Branch("re", &re, "re/D");
    matrix_tree->Branch("re", &re, "re/D");

    for (Eigen::Index jdx = 0; jdx < nBasisElts; jdx++) {
        for (Eigen::Index idx = 0; idx <= jdx; idx++) {
            val = op(idx, jdx);
            re = std::real(val);
            im = std::imag(val);
            mag = std::abs(val);
            if (mag > threshold) {
                matrix_tree->Fill();
            } else if (mag > 0) {
                nFails++;
            } else {
                nZeros++;
            }
        }
    }
    std::cout << fmt::format("Filled {} entries of {} ({} %), {} zeros, {} fails", nFills, nElts, (100 * nFills) / (double)nElts, nZeros, nFails) << std::endl;
    matrix_tree->Write();
    return matrix_tree;
}

TTree* aef::write_matrix_tree(Eigen::MatrixXcd& op, const char* name, const char* title, double mag_thresh, bool use_rel_thres, bool cartesian) {
    TTree* matrix_tree = new TTree(name, title);

    const Eigen::Index nBasisElts = op.rows();

    double threshold = mag_thresh;

    if (use_rel_thres) {
        double min_mag = std::numeric_limits<double>::infinity();
        double max_mag = -1;
        for (Eigen::Index jdx = 0; jdx < nBasisElts; jdx++) {
            for (Eigen::Index idx = 0; idx < nBasisElts; idx++) {
                aef::dcomplex val = op(idx, jdx);
                double mag = std::norm(val);
                min_mag = std::min(min_mag, mag);
                max_mag = std::max(max_mag, mag);
            }
        }
        threshold = mag_thresh * max_mag;
    }

    Eigen::Index idx, jdx;
    double mag, phase;
    double real, imag;
    aef::dcomplex val;

    matrix_tree->Branch("idx", &idx, "idx/I");
    matrix_tree->Branch("jdx", &jdx, "jdx/I");
    matrix_tree->Branch("mag", &mag, "mag/D");
    //matrix_tree->Branch("val", &val, "val/D");
    matrix_tree->Branch("phase", &phase, "phase/D");

    size_t nFills = 0;
    size_t nZeros = 0;
    size_t nFails = 0;
    size_t nElts = nBasisElts * nBasisElts;
    for (jdx = 0; jdx < nBasisElts; jdx++) {
        for (idx = 0; idx < nBasisElts; idx++) {
            val = op(idx, jdx);
            mag = std::abs(val); // note: norm --> mag squared
            phase = std::arg(val);
            // only fill
            if (mag > threshold) {
                matrix_tree->Fill();
                nFills++;
            } else if (mag > 0) {
                nFails++;
            } else {
                nZeros++;
            }
        }
    }
    std::cout << fmt::format("Filled {} entries of {} ({} %), {} zeros, {} fails", nFills, nElts, (100 * nFills) / (double)nElts, nZeros, nFails) << std::endl;
    matrix_tree->Write();
    return matrix_tree;
}
