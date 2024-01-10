#include "TStyle.h"
#include "TMath.h"
#include "TCanvas.h"
#include "TPaveText.h"
#include <iostream>
#include <memory>
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
#include <TPaletteAxis.h>
using namespace std;
// #include "sndisplay.cc"


void MC_Simu(){

  TFile *file = new TFile("data/simu_Bi_207.root", "READ");

  std::vector<int> *om_num = new std::vector<int>;
  std::vector<int> *calo_nohits = new std::vector<int>;
  std::vector<double> *energy_ubc = new std::vector<double>;
  std::vector<double> *energy_u = new std::vector<double>;
  std::vector<double> *energy_bc = new std::vector<double>;
  std::vector<double> *energy = new std::vector<double>;
  std::vector<double> *calo_timestamp = new std::vector<double>;
  std::vector<int> *tracker_nohits = new std::vector<int>;
  std::vector<double> *tracker_side = new std::vector<double>;
  std::vector<double> *tracker_layer = new std::vector<double>;
  std::vector<double> *tracker_column = new std::vector<double>;
  std::vector<double> *x_gg = new std::vector<double>;
  std::vector<double> *y_gg = new std::vector<double>;
  std::vector<double> *z_gg = new std::vector<double>;

  TRandom3 rando;

  TTree* tree = (TTree*)file->Get("Result_tree");
  tree->SetBranchStatus("*",0);
  tree->SetBranchStatus("calo_nohits",1);
  tree->SetBranchAddress("calo_nohits", &calo_nohits);
  tree->SetBranchStatus("om_id",1);
  tree->SetBranchAddress("om_id", &om_num);
  tree->SetBranchStatus("energy_ubc",1);
  tree->SetBranchAddress("energy_ubc", &energy_ubc);
  tree->SetBranchStatus("energy_u",1);
  tree->SetBranchAddress("energy_u", &energy_u);
  tree->SetBranchStatus("energy_bc",1);
  tree->SetBranchAddress("energy_bc", &energy_bc);
  tree->SetBranchStatus("energy",1);
  tree->SetBranchAddress("energy", &energy);
  tree->SetBranchStatus("calo_timestamp",1);
  tree->SetBranchAddress("calo_timestamp", &calo_timestamp);
  tree->SetBranchStatus("tracker_nohits",1);
  tree->SetBranchAddress("tracker_nohits", &tracker_nohits);
  tree->SetBranchStatus("tracker_side",1);
  tree->SetBranchAddress("tracker_side", &tracker_side);
  tree->SetBranchStatus("tracker_layer",1);
  tree->SetBranchAddress("tracker_layer", &tracker_layer);
  tree->SetBranchStatus("tracker_column",1);
  tree->SetBranchAddress("tracker_column", &tracker_column);
  tree->SetBranchStatus("x_gg",1);
  tree->SetBranchAddress("x_gg", &x_gg);
  tree->SetBranchStatus("y_gg",1);
  tree->SetBranchAddress("y_gg", &y_gg);
  tree->SetBranchStatus("z_gg",1);
  tree->SetBranchAddress("z_gg", &z_gg);

  std::vector<double> energy_new;
  std::vector<double> energy_u_new;
  std::vector<double> energy_bc_new;
  std::vector<double> energy_ubc_new;

  TFile *newfile = new TFile("data/simu_Bi_207_eres.root","RECREATE");
  TTree Result_tree("Result_tree","");
  Result_tree.Branch("calo_nohits", &calo_nohits);
  Result_tree.Branch("om_number", &om_num);
  Result_tree.Branch("energy", &energy_new);
  Result_tree.Branch("energy_u", &energy_u_new);
  Result_tree.Branch("energy_bc", &energy_bc_new);
  Result_tree.Branch("energy_ubc", &energy_ubc_new);
  Result_tree.Branch("calo_timestamp", &calo_timestamp);
  Result_tree.Branch("tracker_nohits", &tracker_nohits);
  Result_tree.Branch("tracker_side", &tracker_side);
  Result_tree.Branch("tracker_layer", &tracker_layer);
  Result_tree.Branch("tracker_column", &tracker_column);
  Result_tree.Branch("x_gg", &x_gg);
  Result_tree.Branch("y_gg", &y_gg);
  Result_tree.Branch("z_gg", &z_gg);

  for (int i = 0; i < tree->GetEntries(); i++) {
    if (i%1000000 == 0) {
      std::cout << "/* entry = " << i << '\n';
    }
    double E_kolmo =0;
    tree->GetEntry(i);

    for (int j = 0; j < energy_ubc->size(); j++) {
      float Evis_ubc = rando.Gaus(energy_ubc->at(j), (8/235.482)*sqrt(energy_ubc->at(j)));
      float Evis_u = rando.Gaus(energy_u->at(j), (8/235.482)*sqrt(energy_u->at(j)));
      float Evis_bc = rando.Gaus(energy_bc->at(j), (8/235.482)*sqrt(energy_bc->at(j)));
      float Evis = rando.Gaus(energy->at(j), (8/235.482)*sqrt(energy->at(j)));

      energy_new.push_back(Evis);
      energy_ubc_new.push_back(Evis_ubc);
      energy_bc_new.push_back(Evis_bc);
      energy_u_new.push_back(Evis_u);


    }
    Result_tree.Fill();
    energy_new.clear();
    energy_ubc_new.clear();
    energy_u_new.clear();
    energy_bc_new.clear();
  }


  newfile->cd();
  Result_tree.Write();
  newfile->Close();

  return;
}
