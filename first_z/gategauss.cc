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


void gauss_fitter() {
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

  std::vector<double> *z_last_gg = new std::vector<double>;
  int nelec, ncalo_tot;

  TFile tree_file("../test.root", "READ");
  TTree* tree = (TTree*)tree_file.Get("Result_tree");
  tree->SetBranchStatus("*",0);
  tree->SetBranchStatus("z_first_gg",1);
  tree->SetBranchAddress("z_first_gg", &z_last_gg);
  tree->SetBranchStatus("first_column",1);
  tree->SetBranchAddress("first_column", &first_column);
  tree->SetBranchStatus("nelec",1);
  tree->SetBranchAddress("nelec", &nelec);
  tree->SetBranchStatus("ncalo_tot",1);
  tree->SetBranchAddress("ncalo_tot", &ncalo_tot);
  gROOT->cd();

  for(int i = 0; i < 6; i++){                   //loop on source column
    for(int i = 0; i < 7; i++){                 //loop on source row
      TH1D *z_distrib = new TH1D ("z_distrib", "", 200, -1, 1);
      tree->Project("z_distrib", "z_first_gg", Form("z_first_gg > -10 && ncalo_tot == 1 && first_column < 68 && first_column > 62", row_number));
      // tree->Project("z_distrib", "z_last_gg", Form("z_last_gg > -1.5 && om_number == %d && om_number%%13 == 6 && nelec == 1 && ncalo_tot == 1", 6 + 13*row_number));

      TCanvas* canvas = new TCanvas;
      TF1 *f_gategauss = new TF1 ("f_gategauss", gategaus_func, -1, 1, 4);
      // f_gategauss->SetParameters(0.1, 0.01, 200, 0.03);
      f_gategauss->SetParameters(0.1, z_distrib->GetMean(), z_distrib->GetMaximum(), 0.03);
      f_gategauss->SetParNames("gate sigma", "gate mean", "normalization", "gaus sigma");
      f_gategauss->SetNpx(1000);
      // z_distrib->Draw();
      // f_gategauss->Draw("same");
      // return;
      for (size_t t = 0; t < 400; t++) {
        z_distrib->Fit("f_gategauss", "R0");
        f_gategauss->SetRange(f_gategauss->GetParameter(1) - 2*f_gategauss->GetParameter(0), f_gategauss->GetParameter(1) + 2*f_gategauss->GetParameter(0));
        f_gategauss->SetParameters(f_gategauss->GetParameter(0), f_gategauss->GetParameter(1), f_gategauss->GetParameter(2), f_gategauss->GetParameter(3));
      }

      z_distrib->Draw();
      f_gategauss->Draw("lsame");
      // return;
      // canvas->SetLogy();
      // z_distrib->GetYaxis()->SetRangeUser(0.1,10000);

      sigma = f_gategauss->GetParameter(3);
      sigma_error = f_gategauss->GetParError(3);
      mean = f_gategauss->GetParameter(1);
      mean_error = f_gategauss->GetParError(1);
      Khi2 = f_gategauss->GetChisquare()/f_gategauss->GetNDF();
      std::cout << Khi2 << i << '\n';
      Result_tree.Fill();
      // canvas->SaveAs(Form("png_om/z_distribution_om_%d.png", i));
      canvas->SaveAs(Form("png/z_distribution_row_%d.png", i));
      delete canvas;
      delete f_gategauss;

      delete z_distrib;
      // return;
    }
  }

  TFile file("z_distrib.root", "RECREATE");
  file.cd();
  Result_tree.Write();
  file.Close();

}

int main(int argc, char const *argv[]) {

  return 0;
}
