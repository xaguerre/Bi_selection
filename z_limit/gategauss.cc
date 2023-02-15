#include "TStyle.h"
#include "TMath.h"
#include "TCanvas.h"
#include "TPaveText.h"
#include <iostream>
#include <memory>
using namespace RooFit;
#define _USE_MATH_DEFINES
#include <fstream>
#include <sstream>
#include <vector>
#include <TH1D.h>
#include <TH2D.h>
#include <TH3D.h>
#include <TF1.h>
#include <TStyle.h>
#include "TKey.h"
#include "TFile.h"
#include "TTree.h"
#include "TLine.h"
#include "TROOT.h"
#include <TText.h>
#include <TLatex.h>
#include <TRandom3.h>
#include <TLegend.h>
#include <TSystem.h>
#include <TGraphErrors.h>
#include <TGraphAsymmErrors.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <TParameter.h>
using namespace std;

double gategaus_func (double *x, double *par)
{
  const double GATE_WIDTH = par[0];
  const double GATE_MEAN = par[1];
  const double TOT_INT = par[2];
  const double GAUS_WIDTH = par[3];

  //Fit parameters:
  //par[0]= Sigma of the gate
  //par[1]= Mean of the gate
  //par[2]=Total area (integral -inf to inf, normalization constant)
  //par[3]=Width (sigma) of convoluted Gaussian function

  // Numeric constants
  double invsq2pi = 0.3989422804014;   // (2 pi)^(-1/2)

  // Control constants
  double np = 100.0;      // number of convolution steps
  double sc = 5.0;      // convolution extends to +-sc Gaussian sigmas

  // Variables
  double xx;
  double fgate;
  double sum = 0.0;
  double xlow,xupp;
  double step;
  double i;

  // Range of convolution integral
  xlow = x[0] - sc * GAUS_WIDTH;
  xupp = x[0] + sc * GAUS_WIDTH;

  step = (xupp-xlow) / np;
  // Convolution integral of Landau and Gaussian by sum
  for(i=1.0; i<=np/2; i++) {
    xx = xlow + (i-.5) * step;
    // fland = TMath::Landau(mpc - xx,0,LANDAU_WIDTH) / LANDAU_WIDTH ;
    fgate = (TMath::Abs(xx - GATE_MEAN) < GATE_WIDTH) ? 1:0;
    sum += fgate * TMath::Gaus(x[0],xx,GAUS_WIDTH);

    xx = xupp - (i-.5) * step;
    // fland = TMath::Landau(mpc - xx,0,LANDAU_WIDTH) / LANDAU_WIDTH ;
    fgate = (TMath::Abs(xx - GATE_MEAN) < GATE_WIDTH) ? 1:0;
    sum += fgate * TMath::Gaus(x[0],xx,GAUS_WIDTH);
  }

  return (TOT_INT * step * sum * invsq2pi / GAUS_WIDTH);
}

void gategauss() {
  gStyle->SetOptStat(1111);
  gStyle->SetOptFit(111);
  gStyle->SetLabelSize(0.03,"x");
  gStyle->SetLabelSize(0.03,"y");

  int om_number, row_number;
  double mean, mean_error, Khi2, sigma, sigma_error;
  TTree Result_tree("Result_tree","");
  // Result_tree.Branch("run_number", &run_number);
  Result_tree.Branch("mean", &mean);
  Result_tree.Branch("mean_error", &mean_error);
  Result_tree.Branch("sigma", &sigma);
  Result_tree.Branch("sigma_error", &sigma_error);
  Result_tree.Branch("Khi2", &Khi2);
  Result_tree.Branch("row_number", &row_number);

  std::vector<int> *vec_om_number = new std::vector<int>;
  std::vector<double> *z_last_gg = new std::vector<double>;
  int nelec, ncalo_tot;

  TFile tree_file("../test2.root", "READ");
  TTree* tree = (TTree*)tree_file.Get("Result_tree");
  tree->SetBranchStatus("*",0);
  tree->SetBranchStatus("om_number",1);
  tree->SetBranchAddress("om_number", &vec_om_number);
  tree->SetBranchStatus("z_last_gg",1);
  tree->SetBranchAddress("z_last_gg", &z_last_gg);
  tree->SetBranchStatus("nelec",1);
  tree->SetBranchAddress("nelec", &nelec);
  tree->SetBranchStatus("ncalo_tot",1);
  tree->SetBranchAddress("ncalo_tot", &ncalo_tot);
  gROOT->cd();

  for (int i = 2; i < 12; i++) {
    row_number = i;
    TH1D *z_distrib = new TH1D ("z_distrib", "", 200, -1, 1);
    tree->Project("z_distrib", "z_last_gg", Form("z_last_gg > -1.5 && om_number%%13 == %d && nelec == 1 && ncalo_tot == 1", row_number));
    // tree->Project("z_distrib", "charge_tree", Form("charge_tree > 2500 && om_number == %d && time < %f", i, t + 9999999*(i+1)));
    // z_distrib->Draw();
    // return;

    // z_distrib->Scale(1/z_distrib->GetEntries());
    TCanvas* canvas = new TCanvas;
    TF1 *f_gategauss = new TF1 ("f_gategauss", gategaus_func, -1, 1, 4);
    // f_gategauss->SetParameters(0.1, 0.01, 200, 0.03);
    f_gategauss->SetParameters(0.1, z_distrib->GetMean(), z_distrib->GetMaximum(), 0.03);
    f_gategauss->SetParNames("gate sigma", "gate mean", "normalization", "gaus sigma");
    f_gategauss->SetNpx(1000);
    // z_distrib->Draw();
    // f_gategauss->Draw("same");
    // return;


    z_distrib->Fit("f_gategauss", "R0");
    // z_distrib->Fit("f_gategauss", "R");
    // f_gategauss->SetRange(f_gategauss->GetParameter(1)-4*f_gategauss->GetParameter(3), f_gategauss->GetParameter(1)+4*f_gategauss->GetParameter(3));
    // z_distrib->Fit("f_gategauss", "R");
    // f_gategauss->SetRange(f_gategauss->GetParameter(1)-4*f_gategauss->GetParameter(3), f_gategauss->GetParameter(1)+4*f_gategauss->GetParameter(3));

    z_distrib->Draw();
    f_gategauss->Draw("lsame");
    return;
    // canvas->SetLogy();
    // z_distrib->GetYaxis()->SetRangeUser(0.1,10000);

    sigma = f_gategauss->GetParameter(0);
    sigma_error = f_gategauss->GetParError(0);
    mean = f_gategauss->GetParameter(1);
    mean_error = f_gategauss->GetParError(1);
    Khi2 = f_gategauss->GetChisquare()/f_gategauss->GetNDF();
    std::cout << Khi2 << i << '\n';
    Result_tree.Fill();

    canvas->SaveAs(Form("png/z_distribution_row_%d.png", i));
    delete canvas;
    delete f_gategauss;

    delete z_distrib;
    // return;
  }

  TFile file("z_distrib.root", "RECREATE");
  file.cd();
  Result_tree.Write();
  file.Close();

}


int main(int argc, char const *argv[]) {

  return 0;
}
