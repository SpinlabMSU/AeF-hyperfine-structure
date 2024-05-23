#include <string>
#include <TDirectory.h>
#include <TTree.h>
#include <TH2.h>

#define eprintf(...) fprintf(stderr, __VA_ARGS__)

constexpr double MHz_D_per_V_m = 0.005034;
/// <summary>
/// Electric field conversion factor from V/cm to MHz/D
/// </summary>
constexpr double MHz_D_per_V_cm = 100 * MHz_D_per_V_m;
constexpr double V_cm_per_MHz_D = 1.0 / MHz_D_per_V_cm;

// Root C++ file
struct run_params {
    double nmax; // units: hbar
    double K; // units:
    double E_z; // units: V/cm
    bool enableDev; // bool
    bool starkOnly;
    std::string *run; // run name string
};


run_params get_run_params(TFile &f) {
    TDirectory* dir = f.GetDirectory("parameters");
    if (!dir) {
        printf("Error directory parameters was null\n");
        exit(11111);
    }
    struct run_params params = {};
    params.nmax = dir->Get<TParameter<double>>("nmax")->GetVal();
    params.K = dir->Get<TParameter<double>>("K")->GetVal();
    params.E_z = dir->Get<TParameter<double>>("E_z")->GetVal() * V_cm_per_MHz_D;
    params.enableDev = dir->Get<TParameter<bool>>("enableDev")->GetVal();
    params.starkOnly = dir->Get<TParameter<bool>>("stark_only")->GetVal();
    //params.run = new std::string("Run <placeholder>");
    params.run = dir->Get<std::string>("run");
    printf("nmax %f, K %f, E_z %f, enableDev %d, run %s\n, starkOnly %d", params.nmax, params.K, params.E_z, params.enableDev, params.run->c_str(), params.starkOnly);
    return params;
}

void plot_Hamiltonian() {
    printf("aaaaa\n");
    TString fname("matrix.root");
    TFile file(fname, "READ");
    file.cd();

    auto H_tot = (TTree*)file.Get("H_tot");
    struct run_params params = get_run_params(file);

    TFile *f2 = new TFile("out2.root", "RECREATE");
    f2->cd();
    gStyle->SetOptStat(0);
    // alignment = 10 * horiz + vert
    // horiz: 1=left adjusted, 2=centered, 3=right adjusted
    // vert: 1=bottom adjusted, 2=centered, 3=top adjusted
    gStyle->SetTitleAlign(23);
    printf("progress #\n");
    // plot ham
    auto *canv = new TCanvas("c1", "ROOT TCanvas", 1920, 1080);
    canv->cd();
    H_tot->Draw("jdx:idx:mag>>hHamiltonian", "", "col");
    printf("progress # 1\n");
    auto hHamiltonian = (TH2F*)gDirectory->Get("hHamiltonian");
    hHamiltonian->SetTitle("Absolute values of the matrix elements of ^{138}BaF X^{2}#Sigma^{+} Rotohyperfine Hamiltonian as expressed in the j-basis");
    hHamiltonian->SetXTitle("index of the j-basis basis-ket");
    hHamiltonian->SetYTitle("index of the j-basis basis-bra");
    printf("progress # 2\n");
    hHamiltonian->Draw("colz");
    printf("progress # 3\n");
    TLatex l;
    l.SetTextSize(0.025);
    const char* endis = "";// params.enableDev ? "in-medium" : "in-vacuum";
    if (params.enableDev) {
        endis = "in-medium";
    } else if (params.starkOnly) {
        endis = "stark-only";
    } else {
        endis = "in-vacuum";
    }
    TString str = TString::Format("n_{max}=%f, %s, E_z=%f #frac{kV}{cm}, run %s", params.nmax, endis, params.E_z/1000, params.run->c_str());
    printf("MF nmax %f, K %f, E_z %f, enableDev %d, run %s\n", params.nmax, params.K, params.E_z, params.enableDev, params.run->c_str());
    printf("String data is %s\n", str.Data());
    printf("progress # 3.5\n");
    l.DrawLatex(980, 280, str.Data());
    printf("progress # 4\n");

    canv->Update();
    printf("progress # 4.5\n");
    canv->SaveAs("Hamiltonian-plot.pdf");
    canv->SaveAs("Hamiltonian-plot.png");
    eprintf("progress # 5\n");
}