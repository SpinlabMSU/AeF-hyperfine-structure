#include <string>
#include <TDirectory.h>
#include <TTree.h>
#include <TH2.h>

#define eprintf(...) fprintf(stderr, __VA_ARGS__)

// Root C++ file
struct run_params {
    double nmax;
    double K;
    double E_z;
    bool enableDev;
    std::string *run;
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
    params.E_z = dir->Get<TParameter<double>>("E_z")->GetVal();
    params.K = dir->Get<TParameter<bool>>("enableDev")->GetVal();
    //params.run = new std::string("Run <placeholder>");
    params.run = dir->Get<std::string>("run");
    return params;
}

void plot_Hamiltonian() {
    printf("aaaaa\n");
    TString fname("matrix.root");
    TFile file(fname, "READ");
    file.cd();

    auto H_tot = (TTree*)file.Get("H_tot");
    struct run_params params = get_run_params(file);

    TFile *f2 = new TFile("out.root", "RECREATE");
    f2->cd();
    gStyle->SetOptStat(0);
    printf("progress #\n");
    // plot ham
    auto *canv = new TCanvas("c1", "ROOT TCanvas", 1920, 1080);
    canv->cd();
    H_tot->Draw("jdx:idx:mag>>hHamiltonian", "", "col");
    printf("progress # 1\n");
    auto hHamiltonian = (TH2F*)gDirectory->Get("hHamiltonian");
    hHamiltonian->SetTitle("^{138}BaF X^{2}#Sigma^{+} Rotohyperfine Hamiltonian in matrix form as expressed in the j-basis");
    hHamiltonian->SetXTitle("index of the j-basis basis-ket");
    hHamiltonian->SetYTitle("index of the j-basis basis-bra");
    printf("progress # 2\n");
    hHamiltonian->Draw("colz");
    printf("progress # 3\n");
    TLatex l;
    l.SetTextSize(0.025);
    const char* endis = params.enableDev ? "in-medium" : "in-vacuum";
    TString str = TString::Format("n_{max}=%d, %s, E_z=%d #frac{kV}{cm}, run %s", params.nmax, endis, params.E_z, params.run->c_str());
    printf("progress # 3.5\n");
    l.DrawLatex(980, 280, str.Data());
    printf("progress # 4\n");

    canv->Update();
    printf("progress # 4.5\n");
    canv->SaveAs("Hamiltonian-plot.pdf");
    canv->SaveAs("Hamiltonian-plot.png");
    eprintf("progress # 5\n");
}