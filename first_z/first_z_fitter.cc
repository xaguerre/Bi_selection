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


double z_bas[7] = {-1,-0.75,-0.45, -0.15, 0.18, 0.5, 0.8};
double z_haut[7] = {-0.8,-0.5,-0.15, 0.15, 0.4, 0.7, 1};
int colonne_gauche[6] = {4, 22, 42, 60, 78, 97};
int colonne_droite[6] = {14, 33, 52, 72, 91, 110};


void gauss_fitter() {
  gStyle->SetOptStat(1111);
  gStyle->SetOptFit(111);
  gStyle->SetLabelSize(0.03,"x");
  gStyle->SetLabelSize(0.03,"y");

  int row_number, column_number;
  double mean, mean_error, Khi2, sigma, sigma_error;
  TTree Result_tree("Result_tree","");
  Result_tree.Branch("mean", &mean);
  Result_tree.Branch("mean_error", &mean_error);
  Result_tree.Branch("sigma", &sigma);
  Result_tree.Branch("sigma_error", &sigma_error);
  Result_tree.Branch("Khi2", &Khi2);
  Result_tree.Branch("row_number", &row_number);
  Result_tree.Branch("column_number", &column_number);

  std::vector<double> *z_first_gg = new std::vector<double>;
  std::vector<double> *first_column = new std::vector<double>;
  int nelec, ncalo_tot;
  std::vector<int> *flag_e_event = new std::vector<int>;
  std::vector<int> *flag_associated_nohit = new std::vector<int>;


  TFile tree_file("../cut_974.root", "READ");
  TTree* tree = (TTree*)tree_file.Get("Result_tree");
  tree->SetBranchStatus("*",0);
  tree->SetBranchStatus("z_first_gg",1);
  tree->SetBranchAddress("z_first_gg", &z_first_gg);
  tree->SetBranchStatus("first_column",1);
  tree->SetBranchAddress("first_column", &first_column);
  tree->SetBranchStatus("flag_e_event",1);
  tree->SetBranchAddress("flag_e_event", &flag_e_event);
  tree->SetBranchStatus("flag_associated_nohit",1);
  tree->SetBranchAddress("flag_associated_nohit", &flag_associated_nohit);
  gROOT->cd();

  for(int i = 0; i < 6; i++){                   //loop on source column
    for(int j = 0; j < 7; j++){                 //loop on source row
      column_number = i;
      row_number = j;
      TH1D *z_distrib = new TH1D ("z_distrib", "", 100, -1, 1);
      tree->Project("z_distrib", "z_first_gg", Form("z_first_gg > %f && z_first_gg < %f && first_column > %d && first_column < %d && flag_e_event == 1 && flag_associated_nohit == 1", z_bas[j], z_haut[j], colonne_gauche[i], colonne_droite[i]));

      TCanvas* canvas = new TCanvas;
      TF1 *f_gauss = new TF1 ("f_gauss", "gaus(0)", z_bas[j], z_haut[j]);
      // f_gauss->SetParameters(0.1, 0.01, 200, 0.03);
      f_gauss->SetParameters(z_distrib->GetMaximum(), z_distrib->GetMean(), 0.05); //);
      f_gauss->SetParNames("constant", "mean", "sigma");
      f_gauss->SetNpx(1000);
      // z_distrib->Draw();
      // f_gauss->Draw("same");
      // return;
      for (size_t t = 0; t < 100; t++) {
        z_distrib->Fit("f_gauss", "R0");
        f_gauss->SetRange(f_gauss->GetParameter(1) - 2*f_gauss->GetParameter(2), f_gauss->GetParameter(1) + 2*f_gauss->GetParameter(2));
        f_gauss->SetParameters(f_gauss->GetParameter(0), f_gauss->GetParameter(1), f_gauss->GetParameter(2));
      }

      z_distrib->Draw();
      f_gauss->Draw("lsame");
      // return;
      // return;
      // canvas->SetLogy();
      // z_distrib->GetYaxis()->SetRangeUser(0.1,10000);

      sigma = f_gauss->GetParameter(2);
      sigma_error = f_gauss->GetParError(2);
      mean = f_gauss->GetParameter(1);
      mean_error = f_gauss->GetParError(1);
      Khi2 = f_gauss->GetChisquare()/f_gauss->GetNDF();
      std::cout << Khi2 << i << '\n';
      Result_tree.Fill();
      // canvas->SaveAs(Form("png_om/z_distribution_om_%d.png", i));
      canvas->SaveAs(Form("png/z_distribution_column_%d_row_%d.png", i, j));
      delete canvas;
      delete f_gauss;

      delete z_distrib;
      // return;
    }
  }

  TFile file("z_distrib.root", "RECREATE");
  file.cd();
  Result_tree.Write();
  file.Close();

}

void moyenneur(){
  TFile *file = new TFile("z_distrib.root", "READ");

  double mean, mean_error, sigma, sigma_error, Khi2;
  int row_number, column_number;

  TTree* tree = (TTree*)file->Get("Result_tree");
  tree->SetBranchStatus("*",0);
  tree->SetBranchStatus("mean",1);
  tree->SetBranchAddress("mean", &mean);
  tree->SetBranchStatus("mean_error",1);
  tree->SetBranchAddress("mean_error", &mean_error);
  tree->SetBranchStatus("sigma",1);
  tree->SetBranchAddress("sigma", &sigma);
  tree->SetBranchStatus("sigma_error",1);
  tree->SetBranchAddress("sigma_error", &sigma_error);
  tree->SetBranchStatus("Khi2",1);
  tree->SetBranchAddress("Khi2", &Khi2);
  tree->SetBranchStatus("row_number",1);
  tree->SetBranchAddress("row_number", &row_number);
  tree->SetBranchStatus("column_number",1);
  tree->SetBranchAddress("column_number", &column_number);

  double delta_z_tab[tree->GetEntries()];
  double error_tab[tree->GetEntries()];


  int entry = tree->GetEntries();

  for (int i = 0; i < entry; i++) {
    tree->GetEntry(i);
    delta_z_tab[i] =  2*3*sigma;
    error_tab[i] = 6*sigma_error;
    std::cout << "sigma = " << sigma << "+-" << sigma_error  << '\n';
  }

  double moyenne = 0;
  double moyenne_plus_moins = 0;

  for (int i = 0; i < entry ; i++) {
    moyenne += delta_z_tab[i];
    moyenne_plus_moins += error_tab[i];
  }
  cout << moyenne << "     " << entry << endl;
  cout << "delta z moyen = " << moyenne/entry << " +- " << moyenne_plus_moins/entry << endl;

}

void TGrapher() {
  TFile *file = new TFile("z_distrib.root", "READ");

  int row_number;
  double mean, mean_error;

  TTree* tree = (TTree*)file->Get("Result_tree");
  tree->SetBranchStatus("*",0);
  tree->SetBranchStatus("mean",1);
  tree->SetBranchAddress("mean", &mean);
  tree->SetBranchStatus("mean_error",1);
  tree->SetBranchAddress("mean_error", &mean_error);
  tree->SetBranchStatus("row_number",1);
  tree->SetBranchAddress("row_number", &row_number);

  double yaxis[tree->GetEntries()];
  double yaxis_error[tree->GetEntries()];
  double xaxis[tree->GetEntries()];
  double xaxis_error[tree->GetEntries()];

  for (int i = 0; i < tree->GetEntries(); i++) {
    tree->GetEntry(i);
    yaxis[i] = mean*2.9/2;
    yaxis_error[i] = mean_error*2.9/2;
    xaxis[i] = row_number;
    xaxis_error[i] = 0;
  }

  TGraphErrors *gain_graph = new TGraphErrors(tree->GetEntries(), xaxis, yaxis, xaxis_error, yaxis_error);
  gain_graph->GetXaxis()->SetTitle("Row");
  gain_graph->GetYaxis()->SetTitle("Source position (cm)");
  gain_graph->Draw();
}


int main(int argc, char const *argv[]) {

  return 0;
}
