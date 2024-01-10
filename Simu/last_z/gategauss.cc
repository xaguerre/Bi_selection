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

TF1 *f_gategauss, *f_gategauss_2, *f_gategauss_3, *f_gategauss_4, *f_gategauss_5, *f_gategauss_6, *f_gategauss_7, *f_gategauss_8, *f_gategauss_9, *f_gategauss_10, *f_gategauss_11;

TF1 *f_gategauss_error, *f_gategauss_error2;

double finter(double *x, double*par) {
  if (f_gategauss->EvalPar(x,par) > 0 && f_gategauss_2->EvalPar(x,par) > 0) {
    return TMath::Abs(f_gategauss->EvalPar(x,par) - f_gategauss_2->EvalPar(x,par));
  }
  else return 100000;
}

double finter_error(double *x, double*par) {
  if (f_gategauss_error->EvalPar(x,par) > 0 && f_gategauss_error2->EvalPar(x,par) > 0) {
    return TMath::Abs(f_gategauss_error->EvalPar(x,par) - f_gategauss_error2->EvalPar(x,par));
  }
  else return 100000;
}

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
  std::vector<int> *flag = new std::vector<int>;

  TFile tree_file("../cut_bi.root", "READ");
  TTree* tree = (TTree*)tree_file.Get("Result_tree");
  tree->SetBranchStatus("*",0);
  tree->SetBranchStatus("om_num",1);
  tree->SetBranchAddress("om_num", &vec_om_number);
  tree->SetBranchStatus("z_last_gg",1);
  tree->SetBranchAddress("z_last_gg", &z_last_gg);
  tree->SetBranchStatus("flag",1);
  tree->SetBranchAddress("flag", &flag);
  gROOT->cd();

  for (int i = 1; i < 12; i++) {
    row_number = i;
    TH1D *z_distrib = new TH1D ("z_distrib", "", 200, -1, 1);
    tree->Project("z_distrib", "z_last_gg", Form("z_last_gg > -1.5 && om_num%%13 == %d && (flag == 1111110 ||flag == 11111111)", row_number),"",tree->GetEntries()/10);
    // z_distrib->Draw();
    // return;
    TCanvas* canvas = new TCanvas;
    TF1 *f_gategauss = new TF1 ("f_gategauss", gategaus_func, -1, 1, 4);
    // f_gategauss->SetParameters(0.1, 0.01, 200, 0.03);
    f_gategauss->SetParameters(0.2, z_distrib->GetMean(), z_distrib->GetMaximum(), 0.03);
    f_gategauss->SetParNames("gate sigma", "gate mean", "normalization", "gaus sigma");
    f_gategauss->SetNpx(1000);
    // z_distrib->Draw();
    // f_gategauss->Draw("same");
    // return;
    for (size_t t = 0; t < 400; t++) {
      z_distrib->Fit("f_gategauss", "R0");
      f_gategauss->SetRange(f_gategauss->GetParameter(1) - 1.5*f_gategauss->GetParameter(0), f_gategauss->GetParameter(1) + 1.5*f_gategauss->GetParameter(0));
      f_gategauss->SetParameters(f_gategauss->GetParameter(0), f_gategauss->GetParameter(1), f_gategauss->GetParameter(2), f_gategauss->GetParameter(3));
    }
    if (i == 5) {
      for (size_t t = 0; t < 600; t++) {
        z_distrib->Fit("f_gategauss", "R0");
        f_gategauss->SetRange(f_gategauss->GetParameter(1) - 2*f_gategauss->GetParameter(0), f_gategauss->GetParameter(1) + 2*f_gategauss->GetParameter(0));
        f_gategauss->SetParameters(f_gategauss->GetParameter(0), f_gategauss->GetParameter(1), f_gategauss->GetParameter(2), f_gategauss->GetParameter(3));
      }
    }
    //   z_distrib->Fit("f_gategauss", "R");
    // f_gategauss->SetParameters(f_gategauss->GetParameter(0), f_gategauss->GetParameter(1), f_gategauss->GetParameter(2), f_gategauss->GetParameter(3));
    // z_distrib->Fit("f_gategauss", "R");
    // f_gategauss->SetParameters(f_gategauss->GetParameter(0), f_gategauss->GetParameter(1), f_gategauss->GetParameter(2), f_gategauss->GetParameter(3));

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

  TFile file("z_distrib.root", "RECREATE");
  file.cd();
  Result_tree.Write();
  file.Close();

}

double error_calculation(double par1, double par2, double par3, double par4, double par5, double par6, double par7, double par8){

  f_gategauss_error = new TF1 ("f_gategauss_error", gategaus_func, -1, 1, 4);
  f_gategauss_error->SetParameters(par1, par2, par3, par4);
  f_gategauss_error->SetNpx(1000);
  f_gategauss_error2 = new TF1 ("f_gategauss_error2", gategaus_func, -1, 1, 4);
  f_gategauss_error2->SetParameters(par5, par6, par7, par8);
  f_gategauss_error2->SetNpx(1000);

  TCanvas c;
  f_gategauss_error->Draw();
  f_gategauss_error2->Draw("same");

  TF1 *fint = new TF1("fint",finter_error,-1,1,0);
  double cross_error = fint->GetMinimumX();
  return cross_error;
}

void cross_gategauss() {
  gStyle->SetOptStat(1111);
  gStyle->SetOptFit(111);
  gStyle->SetLabelSize(0.03,"x");
  gStyle->SetLabelSize(0.03,"y");

  int om_number, row_number;
  double mean, mean_error, Khi2, sigma, sigma_error, cross, cross_error_plus, cross_error_moins;
  TTree Result_tree("Result_tree","");
  // Result_tree.Branch("run_number", &run_number);
  Result_tree.Branch("mean", &mean);
  Result_tree.Branch("mean_error", &mean_error);
  Result_tree.Branch("sigma", &sigma);
  Result_tree.Branch("sigma_error", &sigma_error);
  Result_tree.Branch("Khi2", &Khi2);
  Result_tree.Branch("row_number", &row_number);
  Result_tree.Branch("cross", &cross);
  Result_tree.Branch("cross_error_plus", &cross_error_plus);
  Result_tree.Branch("cross_error_moins", &cross_error_moins);

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

  float cross_tab[11] = {0};
  float cross_tab_plus[11] = {0};
  float cross_tab_moins[11] = {0};

  for (int i = 1; i < 11; i++) {
    TCanvas* canvas = new TCanvas;
    row_number = i;
    TH1D *z_distrib = new TH1D ("z_distrib", "", 200, -1, 1);
    tree->Project("z_distrib", "z_last_gg", Form("z_last_gg > -1.5 && om_number%%13 == %d && nelec == 1 && ncalo_tot == 1", row_number));
    // tree->Project("z_distrib", "z_last_gg", Form("z_last_gg > -1.5 && om_number == %d && om_number%%13 == 6 && nelec == 1 && ncalo_tot == 1", 6 + 13*row_number));

    f_gategauss = new TF1 ("f_gategauss", gategaus_func, -1, 1, 4);
    f_gategauss->SetParameters(0.1, z_distrib->GetMean(), z_distrib->GetMaximum(), 0.03);
    f_gategauss->SetParNames("gate sigma", "gate mean", "normalization", "gaus sigma");
    f_gategauss->SetNpx(1000);

    for (size_t t = 0; t < 400; t++) {
      z_distrib->Fit("f_gategauss", "R0");
      f_gategauss->SetParameters(f_gategauss->GetParameter(0), f_gategauss->GetParameter(1), f_gategauss->GetParameter(2), f_gategauss->GetParameter(3));
    }
    if (i == 5) {
      for (size_t t = 0; t < 600; t++) {
        z_distrib->Fit("f_gategauss", "R0");
        f_gategauss->SetParameters(f_gategauss->GetParameter(0), f_gategauss->GetParameter(1), f_gategauss->GetParameter(2), f_gategauss->GetParameter(3));
      }
    }
    TH1D *z_distrib_2 = new TH1D ("z_distrib_2", "", 200, -1, 1);
    tree->Project("z_distrib_2", "z_last_gg", Form("z_last_gg > -1.5 && om_number%%13 == %d && nelec == 1 && ncalo_tot == 1", row_number+1));
    // tree->Project("z_distrib", "z_last_gg", Form("z_last_gg > -1.5 && om_number == %d && om_number%%13 == 6 && nelec == 1 && ncalo_tot == 1", 6 + 13*row_number));
    f_gategauss_2 = new TF1 ("f_gategauss_2", gategaus_func, -1, 1, 4);
    f_gategauss_2->SetParameters(0.1, z_distrib_2->GetMean(), z_distrib_2->GetMaximum(), 0.03);
    f_gategauss_2->SetParNames("gate sigma", "gate mean", "normalization", "gaus sigma");
    f_gategauss_2->SetNpx(1000);

    for (size_t t = 0; t < 400; t++) {
      z_distrib_2->Fit("f_gategauss_2", "R0");
      f_gategauss_2->SetParameters(f_gategauss_2->GetParameter(0), f_gategauss_2->GetParameter(1), f_gategauss_2->GetParameter(2), f_gategauss_2->GetParameter(3));
    }
    if (i == 5) {
      for (size_t t = 0; t < 600; t++) {
        z_distrib_2->Fit("f_gategauss_2", "R0");
        f_gategauss_2->SetParameters(f_gategauss_2->GetParameter(0), f_gategauss_2->GetParameter(1), f_gategauss_2->GetParameter(2), f_gategauss_2->GetParameter(3));
      }
    }
    f_gategauss->SetParameter(2, f_gategauss->GetParameter(2)/f_gategauss->GetMaximum());
    f_gategauss_2->SetParameter(2, f_gategauss_2->GetParameter(2)/f_gategauss_2->GetMaximum());

    // z_distrib->Draw();
    // f_gategauss->Draw();
    // z_distrib_2->Draw("lsame");
    // f_gategauss_2->Draw("lsame");

    TF1 *fint = new TF1("fint",finter,-1,1,0);
    cross = fint->GetMinimumX();
    cross_tab[i] = cross;

    cross_error_plus = abs(cross_tab[i] - error_calculation(f_gategauss->GetParameter(0) + f_gategauss->GetParError(0), f_gategauss->GetParameter(1), f_gategauss->GetParameter(2), f_gategauss->GetParameter(3), f_gategauss_2->GetParameter(0) + f_gategauss_2->GetParError(0), f_gategauss_2->GetParameter(1), f_gategauss_2->GetParameter(2), f_gategauss_2->GetParameter(3)));
    cross_tab_plus[i] = cross_error_plus;
    cross_error_moins = abs(cross_tab[i] - error_calculation(f_gategauss->GetParameter(0) - f_gategauss->GetParError(0), f_gategauss->GetParameter(1), f_gategauss->GetParameter(2), f_gategauss->GetParameter(3), f_gategauss_2->GetParameter(0) - f_gategauss_2->GetParError(0), f_gategauss_2->GetParameter(1), f_gategauss_2->GetParameter(2), f_gategauss_2->GetParameter(3)));
    cross_tab_moins[i] = cross_error_moins;
    // canvas->SetLogy();
    // z_distrib->GetYaxis()->SetRangeUser(0.1,10000);


    sigma = f_gategauss->GetParameter(3);
    sigma_error = f_gategauss->GetParError(3);
    mean = f_gategauss->GetParameter(1);
    mean_error = f_gategauss->GetParError(1);
    Khi2 = f_gategauss->GetChisquare()/f_gategauss->GetNDF();
    Result_tree.Fill();
    TLatex l;
    l.SetTextFont(40);
    l.DrawLatex(0, 0.01, Form("Cross = %.2f", cross));
    canvas->SaveAs(Form("png_cross/z_distribution_row_%d.png", i));
    delete canvas;
    delete f_gategauss;
    delete f_gategauss_2;
    delete z_distrib_2;
    delete z_distrib;

  }


    TCanvas* canvas = new TCanvas;
    TH1D *z_distrib = new TH1D ("z_distrib", "", 200, -1, 1);
    tree->Project("z_distrib", "z_last_gg", Form("z_last_gg > -1.5 && om_number%%13 == %d && nelec == 1 && ncalo_tot == 1", 1));
    f_gategauss = new TF1 ("f_gategauss", gategaus_func, -1, 1, 4);
    f_gategauss->SetParameters(0.1, z_distrib->GetMean(), z_distrib->GetMaximum(), 0.03);
    f_gategauss->SetParNames("gate sigma", "gate mean", "normalization", "gaus sigma");
    f_gategauss->SetNpx(1000);

    for (size_t i = 0; i < 400; i++) {
      z_distrib->Fit("f_gategauss", "R0");
      f_gategauss->SetParameters(f_gategauss->GetParameter(0), f_gategauss->GetParameter(1), f_gategauss->GetParameter(2), f_gategauss->GetParameter(3));
    }

    TH1D *z_distrib_2 = new TH1D ("z_distrib_2", "", 200, -1, 1);
    tree->Project("z_distrib_2", "z_last_gg", Form("z_last_gg > -1.5 && om_number%%13 == %d && nelec == 1 && ncalo_tot == 1", 2));
    f_gategauss_2 = new TF1 ("f_gategauss_2", gategaus_func, -1, 1, 4);
    f_gategauss_2->SetParameters(0.1, z_distrib_2->GetMean(), z_distrib_2->GetMaximum(), 0.03);
    f_gategauss_2->SetParNames("gate sigma", "gate mean", "normalization", "gaus sigma");
    f_gategauss_2->SetNpx(1000);

    for (size_t i = 0; i < 400; i++) {
      z_distrib_2->Fit("f_gategauss_2", "R0");
      f_gategauss_2->SetParameters(f_gategauss_2->GetParameter(0), f_gategauss_2->GetParameter(1), f_gategauss_2->GetParameter(2), f_gategauss_2->GetParameter(3));
    }

    TH1D *z_distrib_3 = new TH1D ("z_distrib_3", "", 200, -1, 1);
    tree->Project("z_distrib_3", "z_last_gg", Form("z_last_gg > -1.5 && om_number%%13 == %d && nelec == 1 && ncalo_tot == 1", 3));
    f_gategauss_3 = new TF1 ("f_gategauss_3", gategaus_func, -1, 1, 4);
    f_gategauss_3->SetParameters(0.1, z_distrib_3->GetMean(), z_distrib_3->GetMaximum(), 0.03);
    f_gategauss_3->SetParNames("gate sigma", "gate mean", "normalization", "gaus sigma");
    f_gategauss_3->SetNpx(1000);

    for (size_t i = 0; i < 400; i++) {
      z_distrib_3->Fit("f_gategauss_3", "R0");
      f_gategauss_3->SetParameters(f_gategauss_3->GetParameter(0), f_gategauss_3->GetParameter(1), f_gategauss_3->GetParameter(2), f_gategauss_3->GetParameter(3));
    }

    TH1D *z_distrib_4 = new TH1D ("z_distrib_4", "", 200, -1, 1);
    tree->Project("z_distrib_4", "z_last_gg", Form("z_last_gg > -1.5 && om_number%%13 == %d && nelec == 1 && ncalo_tot == 1", 4));
    f_gategauss_4 = new TF1 ("f_gategauss_4", gategaus_func, -1, 1, 4);
    f_gategauss_4->SetParameters(0.1, z_distrib_4->GetMean(), z_distrib_4->GetMaximum(), 0.03);
    f_gategauss_4->SetParNames("gate sigma", "gate mean", "normalization", "gaus sigma");
    f_gategauss_4->SetNpx(1000);

    for (size_t i = 0; i < 400; i++) {
      z_distrib_4->Fit("f_gategauss_4", "R0");
      f_gategauss_4->SetParameters(f_gategauss_4->GetParameter(0), f_gategauss_4->GetParameter(1), f_gategauss_4->GetParameter(2), f_gategauss_4->GetParameter(3));
    }

    TH1D *z_distrib_5 = new TH1D ("z_distrib_5", "", 200, -1, 1);
    tree->Project("z_distrib_5", "z_last_gg", Form("z_last_gg > -1.5 && om_number%%13 == %d && nelec == 1 && ncalo_tot == 1", 5));
    f_gategauss_5 = new TF1 ("f_gategauss_5", gategaus_func, -1, 1, 4);
    f_gategauss_5->SetParameters(0.1, z_distrib_5->GetMean(), z_distrib_5->GetMaximum(), 0.03);
    f_gategauss_5->SetParNames("gate sigma", "gate mean", "normalization", "gaus sigma");
    f_gategauss_5->SetNpx(1000);

    for (size_t i = 0; i < 1000; i++) {
      z_distrib_5->Fit("f_gategauss_5", "R0");
      f_gategauss_5->SetParameters(f_gategauss_5->GetParameter(0), f_gategauss_5->GetParameter(1), f_gategauss_5->GetParameter(2), f_gategauss_5->GetParameter(3));
    }

    TH1D *z_distrib_6 = new TH1D ("z_distrib_6", "", 200, -1, 1);
    tree->Project("z_distrib_6", "z_last_gg", Form("z_last_gg > -1.5 && om_number%%13 == %d && nelec == 1 && ncalo_tot == 1", 6));
    f_gategauss_6 = new TF1 ("f_gategauss_6", gategaus_func, -1, 1, 4);
    f_gategauss_6->SetParameters(0.1, z_distrib_6->GetMean(), z_distrib_6->GetMaximum(), 0.03);
    f_gategauss_6->SetParNames("gate sigma", "gate mean", "normalization", "gaus sigma");
    f_gategauss_6->SetNpx(1000);

    for (size_t i = 0; i < 400; i++) {
      z_distrib_6->Fit("f_gategauss_6", "R0");
      f_gategauss_6->SetParameters(f_gategauss_6->GetParameter(0), f_gategauss_6->GetParameter(1), f_gategauss_6->GetParameter(2), f_gategauss_6->GetParameter(3));
    }

    TH1D *z_distrib_7 = new TH1D ("z_distrib_7", "", 200, -1, 1);
    tree->Project("z_distrib_7", "z_last_gg", Form("z_last_gg > -1.5 && om_number%%13 == %d && nelec == 1 && ncalo_tot == 1", 7));
    f_gategauss_7 = new TF1 ("f_gategauss_7", gategaus_func, -1, 1, 4);
    f_gategauss_7->SetParameters(0.1, z_distrib_7->GetMean(), z_distrib_7->GetMaximum(), 0.03);
    f_gategauss_7->SetParNames("gate sigma", "gate mean", "normalization", "gaus sigma");
    f_gategauss_7->SetNpx(1000);

    for (size_t i = 0; i < 400; i++) {
      z_distrib_7->Fit("f_gategauss_7", "R0");
      f_gategauss_7->SetParameters(f_gategauss_7->GetParameter(0), f_gategauss_7->GetParameter(1), f_gategauss_7->GetParameter(2), f_gategauss_7->GetParameter(3));
    }

    TH1D *z_distrib_8 = new TH1D ("z_distrib_8", "", 200, -1, 1);
    tree->Project("z_distrib_8", "z_last_gg", Form("z_last_gg > -1.5 && om_number%%13 == %d && nelec == 1 && ncalo_tot == 1", 8));
    f_gategauss_8 = new TF1 ("f_gategauss_8", gategaus_func, -1, 1, 4);
    f_gategauss_8->SetParameters(0.1, z_distrib_8->GetMean(), z_distrib_8->GetMaximum(), 0.03);
    f_gategauss_8->SetParNames("gate sigma", "gate mean", "normalization", "gaus sigma");
    f_gategauss_8->SetNpx(1000);

    for (size_t i = 0; i < 400; i++) {
      z_distrib_8->Fit("f_gategauss_8", "R0");
      f_gategauss_8->SetParameters(f_gategauss_8->GetParameter(0), f_gategauss_8->GetParameter(1), f_gategauss_8->GetParameter(2), f_gategauss_8->GetParameter(3));
    }

    TH1D *z_distrib_9 = new TH1D ("z_distrib_9", "", 200, -1, 1);
    tree->Project("z_distrib_9", "z_last_gg", Form("z_last_gg > -1.5 && om_number%%13 == %d && nelec == 1 && ncalo_tot == 1", 9));
    f_gategauss_9 = new TF1 ("f_gategauss_9", gategaus_func, -1, 1, 4);
    f_gategauss_9->SetParameters(0.1, z_distrib_9->GetMean(), z_distrib_9->GetMaximum(), 0.03);
    f_gategauss_9->SetParNames("gate sigma", "gate mean", "normalization", "gaus sigma");
    f_gategauss_9->SetNpx(1000);

    for (size_t i = 0; i < 400; i++) {
      z_distrib_9->Fit("f_gategauss_9", "R0");
      f_gategauss_9->SetParameters(f_gategauss_9->GetParameter(0), f_gategauss_9->GetParameter(1), f_gategauss_9->GetParameter(2), f_gategauss_9->GetParameter(3));
    }

    TH1D *z_distrib_10 = new TH1D ("z_distrib_10", "", 200, -1, 1);
    tree->Project("z_distrib_10", "z_last_gg", Form("z_last_gg > -1.5 && om_number%%13 == %d && nelec == 1 && ncalo_tot == 1", 10));
    f_gategauss_10 = new TF1 ("f_gategauss_10", gategaus_func, -1, 1, 4);
    f_gategauss_10->SetParameters(0.1, z_distrib_10->GetMean(), z_distrib_10->GetMaximum(), 0.03);
    f_gategauss_10->SetParNames("gate sigma", "gate mean", "normalization", "gaus sigma");
    f_gategauss_10->SetNpx(1000);

    for (size_t i = 0; i < 400; i++) {
      z_distrib_10->Fit("f_gategauss_10", "R0");
      f_gategauss_10->SetParameters(f_gategauss_10->GetParameter(0), f_gategauss_10->GetParameter(1), f_gategauss_10->GetParameter(2), f_gategauss_10->GetParameter(3));
    }

    TH1D *z_distrib_11 = new TH1D ("z_distrib_11", "", 200, -1, 1);
    tree->Project("z_distrib_11", "z_last_gg", Form("z_last_gg > -1.5 && om_number%%13 == %d && nelec == 1 && ncalo_tot == 1", 11));
    f_gategauss_11 = new TF1 ("f_gategauss_11", gategaus_func, -1, 1, 4);
    f_gategauss_11->SetParameters(0.1, z_distrib_11->GetMean(), z_distrib_11->GetMaximum(), 0.03);
    f_gategauss_11->SetParNames("gate sigma", "gate mean", "normalization", "gaus sigma");
    f_gategauss_11->SetNpx(1000);

    for (size_t i = 0; i < 400; i++) {
      z_distrib_11->Fit("f_gategauss_11", "R0");
      f_gategauss_11->SetParameters(f_gategauss_11->GetParameter(0), f_gategauss_11->GetParameter(1), f_gategauss_11->GetParameter(2), f_gategauss_11->GetParameter(3));
    }


    f_gategauss->SetParameter(2, f_gategauss->GetParameter(2)/f_gategauss->GetMaximum());
    f_gategauss_2->SetParameter(2, f_gategauss_2->GetParameter(2)/f_gategauss_2->GetMaximum());
    f_gategauss_3->SetParameter(2, f_gategauss_3->GetParameter(2)/f_gategauss_3->GetMaximum());
    f_gategauss_4->SetParameter(2, f_gategauss_4->GetParameter(2)/f_gategauss_4->GetMaximum());
    f_gategauss_5->SetParameter(2, f_gategauss_5->GetParameter(2)/f_gategauss_5->GetMaximum());
    f_gategauss_6->SetParameter(2, f_gategauss_6->GetParameter(2)/f_gategauss_6->GetMaximum());
    f_gategauss_7->SetParameter(2, f_gategauss_7->GetParameter(2)/f_gategauss_7->GetMaximum());
    f_gategauss_8->SetParameter(2, f_gategauss_8->GetParameter(2)/f_gategauss_8->GetMaximum());
    f_gategauss_9->SetParameter(2, f_gategauss_9->GetParameter(2)/f_gategauss_9->GetMaximum());
    f_gategauss_10->SetParameter(2, f_gategauss_10->GetParameter(2)/f_gategauss_10->GetMaximum());
    f_gategauss_11->SetParameter(2, f_gategauss_11->GetParameter(2)/f_gategauss_11->GetMaximum());

    f_gategauss->Draw();
    f_gategauss_2->Draw("lsame");
    f_gategauss_3->Draw("lsame");
    f_gategauss_4->Draw("lsame");
    f_gategauss_5->Draw("lsame");
    f_gategauss_6->Draw("lsame");
    f_gategauss_7->Draw("lsame");
    f_gategauss_8->Draw("lsame");
    f_gategauss_9->Draw("lsame");
    f_gategauss_10->Draw("lsame");
    f_gategauss_11->Draw("lsame");

    TH2D DeltaZ("DeltaZ", "DeltaZ", 9,1,10,100,20,30);

    // TLatex l1;
    // l1.SetTextFont(40);
    // l1.SetTextAngle(90);
    // l1.DrawLatex(((cross_tab[0] + cross_tab[1])/2), 0.5, Form("Cross = %.2f cm", (cross_tab[1] - cross_tab[0])*2.9/0.02));

    TLatex l2;
    l2.SetTextFont(40);
    l2.SetTextAngle(90);
    l2.DrawLatex(((cross_tab[1] + cross_tab[2])/2), 0.5, Form("Cross = %.2f cm", (cross_tab[2] - cross_tab[1])*2.9/0.02));

    TLatex l3;
    l3.SetTextFont(40);
    l3.SetTextAngle(90);
    l3.DrawLatex(((cross_tab[2] + cross_tab[3])/2), 0.5, Form("Cross = %.2f cm", (cross_tab[3] - cross_tab[2])*2.9/0.02));

    TLatex l4;
    l4.SetTextFont(40);
    l4.SetTextAngle(90);
    l4.DrawLatex(((cross_tab[3] + cross_tab[4])/2), 0.5, Form("Cross = %.2f cm", (cross_tab[4] - cross_tab[3])*2.9/0.02));

    TLatex l5;
    l5.SetTextFont(40);
    l5.SetTextAngle(90);
    l5.DrawLatex(((cross_tab[4] + cross_tab[5])/2), 0.5, Form("Cross = %.2f cm", (cross_tab[5] - cross_tab[4])*2.9/0.02));

    TLatex l6;
    l6.SetTextFont(40);
    l6.SetTextAngle(90);
    l6.DrawLatex(((cross_tab[5] + cross_tab[6])/2), 0.5, Form("Cross = %.2f cm", (cross_tab[6] - cross_tab[5])*2.9/0.02));

    TLatex l7;
    l7.SetTextFont(40);
    l7.SetTextAngle(90);
    l7.DrawLatex(((cross_tab[6] + cross_tab[7])/2), 0.5, Form("Cross = %.2f cm", (cross_tab[7] - cross_tab[6])*2.9/0.02));

    TLatex l8;
    l8.SetTextFont(40);
    l8.SetTextAngle(90);
    l8.DrawLatex(((cross_tab[7] + cross_tab[8])/2), 0.5, Form("Cross = %.2f cm", (cross_tab[8] - cross_tab[7])*2.9/0.02));

    TLatex l9;
    l9.SetTextFont(40);
    l9.SetTextAngle(90);
    l9.DrawLatex(((cross_tab[8] + cross_tab[9])/2), 0.5, Form("Cross = %.2f cm", (cross_tab[9] - cross_tab[8])*2.9/0.02));

    // TLatex l10;
    // l10.SetTextFont(40);
    // l10.SetTextAngle(90);
    // l10.DrawLatex(((cross_tab[9] + cross_tab[10])/2), 0.5, Form("Cross = %.2f cm", (cross_tab[10] - cross_tab[9])*2.9/0.02));

    // TLatex l11;
    // l11.SetTextFont(40);
    // l11.SetTextAngle(90);
    // l11.DrawLatex(((cross_tab[10] + cross_tab[11])), 0.5, Form("Cross = %.2f cm", cross_tab[11] - cross_tab[10]));


  for (int i = 1; i < 10; i++) {
    cout << "for the line " << i << " z ∈ [" << cross_tab[i] << " - " << cross_tab[i+1] << ", so delta_z = " << (cross_tab[i+1] - cross_tab[i])*2.9/0.02 << "cm"<< endl;
    DeltaZ.Fill(i,(cross_tab[i+1] - cross_tab[i])*2.9/0.02);
  }


  TFile file("cross_z.root", "RECREATE");
  file.cd();
  DeltaZ.Write();
  Result_tree.Write();
  file.Close();

}

void moyenneur(){
  TFile *file = new TFile("cross_z.root", "READ");

  double cross, cross_error_plus, cross_error_moins;

  TTree* tree = (TTree*)file->Get("Result_tree");
  tree->SetBranchStatus("*",0);
  tree->SetBranchStatus("cross",1);
  tree->SetBranchAddress("cross", &cross);
  tree->SetBranchStatus("cross_error_plus",1);
  tree->SetBranchAddress("cross_error_plus", &cross_error_plus);
  tree->SetBranchStatus("cross_error_moins",1);
  tree->SetBranchAddress("cross_error_moins", &cross_error_moins);

  double cross_tab[tree->GetEntries()];
  double cross_tab_plus[tree->GetEntries()];
  double cross_tab_moins[tree->GetEntries()];

  for (int i = 0; i < tree->GetEntries(); i++) {
    tree->GetEntry(i);
    cross_tab[i] = cross*2.9/2;
    cross_tab_plus[i] = cross_error_plus*2.9/2;
    cross_tab_moins[i] = cross_error_moins*2.9/2;
  }

  double moyenne = 0;
  double moyenne_moins = 0;
  double moyenne_plus = 0;
  for (int i = 0; i < 8; i++) {
    moyenne += cross_tab[i+1] - cross_tab[i];
    moyenne_plus += abs(cross_tab_plus[i+1]- cross_tab_plus[i]);
    moyenne_moins += abs(cross_tab_moins[i+1]- cross_tab_moins[i]);
    // cout << "cross = " << cross_tab[i+1] << " + " << cross_tab_plus[i+1] << " - " << cross_tab_moins[i+1] << endl;
  }

  cout << "delta z moyen = " << moyenne/8 << " + " << moyenne_plus/8 << " - " << moyenne_moins/8 << endl;

}

void TGrapher() {
  TFile *file = new TFile("cross_z.root", "READ");

  double cross, cross_error_plus, cross_error_moins;

  TTree* tree = (TTree*)file->Get("Result_tree");
  tree->SetBranchStatus("*",0);
  tree->SetBranchStatus("cross",1);
  tree->SetBranchAddress("cross", &cross);
  tree->SetBranchStatus("cross_error_plus",1);
  tree->SetBranchAddress("cross_error_plus", &cross_error_plus);
  tree->SetBranchStatus("cross_error_moins",1);
  tree->SetBranchAddress("cross_error_moins", &cross_error_moins);

  double yaxis[tree->GetEntries()];
  double yaxis_error_moins[tree->GetEntries()];
  double yaxis_error_plus[tree->GetEntries()];
  double xaxis[tree->GetEntries()];
  double xaxis_error_moins[tree->GetEntries()];
  double xaxis_error_plus[tree->GetEntries()];

  for (int i = 0; i < tree->GetEntries(); i++) {
    tree->GetEntry(i);
    yaxis[i] = cross;
    yaxis_error_plus[i] = cross_error_plus;
    yaxis_error_moins[i] = cross_error_moins;
    xaxis[i] = i+1;
    xaxis_error_moins[i] = 0;
    xaxis_error_plus[i] = 0;
  }

  TGraphAsymmErrors *gain_graph = new TGraphAsymmErrors(tree->GetEntries(), xaxis, yaxis, xaxis_error_moins, xaxis_error_plus, yaxis_error_moins, yaxis_error_plus);
  gain_graph->Draw();



}

int main(int argc, char const *argv[]) {

  return 0;
}
