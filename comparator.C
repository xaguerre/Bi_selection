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


void biplotter() {


  int om_number;
  double mean, sigma, Chi2_ndf;
  double fw_simu[520];

  TFile tree_file("fitted_bi.root", "READ");
  TTree* tree = (TTree*)tree_file.Get("Result_tree");
  tree->SetBranchStatus("*",0);
  tree->SetBranchStatus("om_number",1);
  tree->SetBranchAddress("om_number", &om_number);
  tree->SetBranchStatus("mean",1);
  tree->SetBranchAddress("mean", &mean);
  tree->SetBranchStatus("sigma", 1);
  tree->SetBranchAddress("sigma", &sigma);
  tree->SetBranchStatus("Chi2_ndf",1);
  tree->SetBranchAddress("Chi2_ndf", &Chi2_ndf);

  for (int i = 0; i < tree->GetEntries(); i++) {
    tree->GetEntry(i);
    if (Chi2_ndf > 0 && Chi2_ndf < 3) fw_simu[om_number] = 236*sigma/mean;
    else fw_simu[om_number] = -10;
    }



  TFile tree_file2("../../Bi_fit/fitted_bi.root", "READ");
  TTree* tree2 = (TTree*)tree_file2.Get("Result_tree");
  tree2->SetBranchStatus("*",0);
  tree2->SetBranchStatus("om_number",1);
  tree2->SetBranchAddress("om_number", &om_number);
  tree2->SetBranchStatus("mean",1);
  tree2->SetBranchAddress("mean", &mean);
  tree2->SetBranchStatus("sigma", 1);
  tree2->SetBranchAddress("sigma", &sigma);
  tree2->SetBranchStatus("Chi2_ndf",1);
  tree2->SetBranchAddress("Chi2_ndf", &Chi2_ndf);

  TH2D* MC_spectre_intern = new TH2D("biplot", "biplot", 520, 0, 100, 520, 0, 30);

  for (int i = 0; i < tree->GetEntries(); i++) {
    tree2->GetEntry(i);
    if (Chi2_ndf > 0 && Chi2_ndf < 3) MC_spectre_intern->Fill(fw_simu[i], 236*sigma);
    else MC_spectre_intern->Fill(-10, -10);
  }

  TFile newfile("biplot.root", "RECREATE");
  newfile.cd();
  MC_spectre_intern->Write();
  newfile.Close();

}

void comparator() {
  // SetOptStat(0);

  int ncalo_tot;
  std::vector<int> *flag_e_event = new std::vector<int>;
  std::vector<int> *om_number = new std::vector<int>;
  std::vector<double> *energyvis_ubc = new std::vector<double>;
  std::vector<double> *charge = new std::vector<double>;
  std::vector<int> *flag_charge = new std::vector<int>;
  std::vector<int> *flag_associated_nohit = new std::vector<int>;
  std::vector<int> *flag_MW = new std::vector<int>;
  std::vector<int> *flag_calo_square = new std::vector<int>;
  std::vector<int> *flag_source_square = new std::vector<int>;
  std::vector<int> *flag_last_column = new std::vector<int>;

  TFile *file = new TFile("Simu/cut_bi.root", "READ");
  TTree* tree = (TTree*)file->Get("Result_tree");
  tree->SetBranchStatus("*",0);
  tree->SetBranchStatus("ncalo_tot",1);
  tree->SetBranchAddress("ncalo_tot", &ncalo_tot);
  tree->SetBranchStatus("om_num",1);
  tree->SetBranchAddress("om_num", &om_number);
  tree->SetBranchStatus("flag_e_event",1);
  tree->SetBranchAddress("flag_e_event", &flag_e_event);
  tree->SetBranchStatus("flag_charge",1);
  tree->SetBranchAddress("flag_charge", &flag_charge);
  tree->SetBranchStatus("flag_associated_nohit",1);
  tree->SetBranchAddress("flag_associated_nohit", &flag_associated_nohit);
  tree->SetBranchStatus("flag_MW",1);
  tree->SetBranchAddress("flag_MW", &flag_MW);
  tree->SetBranchStatus("flag_calo_square",1);
  tree->SetBranchAddress("flag_calo_square", &flag_calo_square);
  tree->SetBranchStatus("flag_source_square",1);
  tree->SetBranchAddress("flag_source_square", &flag_source_square);
  tree->SetBranchStatus("flag_last_column",1);
  tree->SetBranchAddress("flag_last_column", &flag_last_column);
  tree->SetBranchStatus("energyvis_ubc",1);
  tree->SetBranchAddress("energyvis_ubc", &energyvis_ubc);


  TH1D *spectre_simu = new TH1D ("spectre_simu", "", 520, 0, 520);
  // tree->Project("spectre_simu", "om_num", "flag_last_column == 1 && ncalo_tot < 2 && energyvis_ubc > 0.6", "", tree->GetEntries()/10);
  // tree->Project("spectre_simu", "om_num", "flag_charge == 1 && flag_associated_nohit == 1 && flag_MW == 1 && flag_calo_square == 1 && flag_source_square == 1 && flag_last_column == 1 && ncalo_tot < 2 && energyvis_ubc > 0.6", "", tree->GetEntries()/10);
  tree->Project("spectre_simu", "om_num", "flag_charge == 1 && (ncalo_tot < 2) && energyvis_ubc > 0.6", "", tree->GetEntries()/10);
  cout << "simu :  " << spectre_simu->Integral() << endl;


  TFile *file2 = new TFile("cut_974.root", "READ");
  TTree* tree2 = (TTree*)file2->Get("Result_tree");
  tree2->SetBranchStatus("*",0);
  tree2->SetBranchStatus("ncalo_tot",1);
  tree2->SetBranchAddress("ncalo_tot", &ncalo_tot);
  tree2->SetBranchStatus("om_number",1);
  tree2->SetBranchAddress("om_number", &om_number);
  tree2->SetBranchStatus("flag_e_event",1);
  tree2->SetBranchAddress("flag_e_event", &flag_e_event);
  tree2->SetBranchStatus("flag_charge",1);
  tree2->SetBranchAddress("flag_charge", &flag_charge);
  tree2->SetBranchStatus("flag_associated_nohit",1);
  tree2->SetBranchAddress("flag_associated_nohit", &flag_associated_nohit);
  tree2->SetBranchStatus("flag_MW",1);
  tree2->SetBranchAddress("flag_MW", &flag_MW);
  tree2->SetBranchStatus("flag_calo_square",1);
  tree2->SetBranchAddress("flag_calo_square", &flag_calo_square);
  tree2->SetBranchStatus("flag_source_square",1);
  tree2->SetBranchAddress("flag_source_square", &flag_source_square);
  tree2->SetBranchStatus("flag_last_column",1);
  tree2->SetBranchAddress("flag_last_column", &flag_last_column);
  tree2->SetBranchStatus("charge",1);
  tree2->SetBranchAddress("charge", &charge);

  TH1D *spectre_data = new TH1D ("spectre_data", "", 520, 0, 520);
  // tree2->Project("spectre_data", "om_number", "flag_last_column == 1 && ncalo_tot < 2 && charge/22000. > 0.6", "", tree2->GetEntries()/10);
  // tree2->Project("spectre_data", "om_number", "flag_charge == 1 && flag_associated_nohit == 1 && flag_MW == 1 && flag_calo_square == 1 && flag_source_square == 1 && flag_last_column == 1 && ncalo_tot < 2 && charge/22000. > 0.6", "", tree2->GetEntries()/10);
  tree2->Project("spectre_data", "om_number", "flag_charge == 1 && ncalo_tot < 2 && charge/22000. > 0.6", "", tree2->GetEntries()/10);
  gROOT->cd();
  cout << "data :  " << spectre_data->Integral() << endl;

  TCanvas* can = new TCanvas("can");
  can->SetBottomMargin(0.3);

  for (int i = 0; i < 520; i++) {
      spectre_data->SetBinError(i, 0);
      spectre_simu->SetBinError(i, 0);
  }

  TLegend *legend = new TLegend(0.3,0.8,0.45,0.9);
  legend->AddEntry(spectre_data, "data", "F");
  legend->AddEntry(spectre_simu, "simu", "F");

  spectre_simu->Draw();
  spectre_simu->SetFillColorAlpha(kBlue, 0.35);
  spectre_data->Scale(1/1.78);
  spectre_data->Draw("same");
  spectre_data->SetLineColor(kRed);
  spectre_data->SetFillColorAlpha(kRed, 0.35);

  legend->Draw("same");

  TPad* pad = new TPad("pad", "pad", 0., 0., 1., 1.);
  pad->SetTopMargin(0.7);
  pad->Draw();
  pad->SetFillStyle(0);
  pad->cd();

  TH1D* residual = new TH1D ("residu", "", 520, 0, 520);
  for(Int_t i = 1; i < 520; i++){
    residual->SetBinContent(i, spectre_simu->GetBinContent(i) / (spectre_data->GetBinContent(i)+1));
  } // i
residual->Draw("ape");

  // can->SaveAs("test.png");

}
//
// void TH2_comparator ()
// {
//   int calo_nohits;
//   std::vector<int> *calo_type = new std::vector<int>;
//   std::vector<int> *calo_side = new std::vector<int>;
//   std::vector<int> *calo_wall = new std::vector<int>;
//   std::vector<int> *calo_column = new std::vector<int>;
//   std::vector<int> *calo_row = new std::vector<int>;
//   std::vector<int> *calo_charge = new std::vector<int>;
//
//   TFile *file = new TFile("data/snemo_run-974_udd.root", "READ");
//   TTree* tree = (TTree*)file->Get("SimData");
//   tree->SetBranchStatus("*",0);
//   tree->SetBranchStatus("digicalo.type",1);
//   tree->SetBranchAddress("digicalo.type", &calo_type);
//   tree->SetBranchStatus("digicalo.side",1);
//   tree->SetBranchAddress("digicalo.side", &calo_side);
//   tree->SetBranchStatus("digicalo.wall",1);
//   tree->SetBranchAddress("digicalo.wall", &calo_wall);
//   tree->SetBranchStatus("digicalo.column",1);
//   tree->SetBranchAddress("digicalo.column", &calo_column);
//   tree->SetBranchStatus("digicalo.row",1);
//   tree->SetBranchAddress("digicalo.row", &calo_row);
//   tree->SetBranchStatus("digicalo.charge",1);
//   tree->SetBranchAddress("digicalo.charge", &calo_charge);
//   tree->SetBranchStatus("digicalo.nohits",1);
//   tree->SetBranchAddress("digicalo.nohits", &calo_nohits);
//
//   double rate[712];
//   memset (rate, 0, 712*sizeof(double));
//   for (int i = 0; i < tree->GetEntries()/10; i++) {
//     tree->GetEntry(i);
//     for (int j = 0; j < calo_type->size(); j++) {
//       if (calo_type->at(j) == 0 && -calo_charge->at(j)/22000. > 0.6 && calo_nohits < 3) {
//         rate[calo_side->at(j)*260 + calo_column->at(j)*13 + calo_row->at(j)]++;
//         // if (calo_side->at(j)*260 + calo_column->at(j)*13 + calo_row->at(j) == 103) cout << "rate 103 = " << rate[103] << endl;
//         // cout << "om : " << calo_side->at(j)*260 + calo_column->at(j)*13 + calo_row->at(j) << " -> charge = " << calo_charge->at(j) << endl;
//       }
//     }
//   }
//
//   // for (size_t i = 0; i < 520; i++) {
//   //   cout << "om " << i << " : " << rate[i] << endl;
//   // }
//
//   std::vector<int> *om_number = new std::vector<int>;
//   std::vector<double> *energyvis_ubc = new std::vector<double>;
//
//   TFile *file2 = new TFile("Simu/data/simu_Bi_207.root", "READ");
//   TTree* tree2 = (TTree*)file2->Get("Result_tree");
//   tree2->SetBranchStatus("*",0);
//   tree2->SetBranchStatus("om_id",1);
//   tree2->SetBranchAddress("om_id", &om_number);
//   tree2->SetBranchStatus("energyvis_ubc",1);
//   tree2->SetBranchAddress("energyvis_ubc", &energyvis_ubc);
//
//
//   double rate_simu[712];
//   memset (rate_simu, 0, 712*sizeof(double));
//   for (int i = 0; i < tree2->GetEntries()/10; i++) {
//     tree2->GetEntry(i);
//     // cout << "size = " << om_number->size() << endl;
//     for (int j = 0; j < om_number->size(); j++) {
//       if (energyvis_ubc->at(j) > 0.6 && om_number->size() < 2) {
//         rate_simu[om_number->at(j)]++;
//
//       }
//     }
//   }
//
//   TH2D *map_it = new TH2D("map_it","map_it",20,0,20,13,0,13);
//   TH2D *map_fr = new TH2D("map_fr","map_fr",20,0,20,13,0,13);
//
//
//   for (int i = 0; i < 520; i++) {
//     cout << "om : " << i << " -> " << rate_simu[i] << " : " << rate[i] << endl;
//   }
//
//   for (int i = 0; i < 260; i++) {
//     map_it->SetBinContent(i/13 +1, i%13 +1, rate_simu[i]/(rate[i]/1.78));
//     map_fr->SetBinContent(i/13 +1, i%13 +1, rate_simu[i+260]/(rate[i+260]/1.78));
//   }
//
//   TFile *newfile = new TFile("test_comparator.root", "RECREATE");
//   newfile->cd();
//   map_it->Write();
//   map_fr->Write();
//   newfile->Close();
// }
//
//
// void TH2_comparator_flag ()
// {
//   gStyle->SetOptStat(0);
//
//   int ncalo_tot;
//   std::vector<int> *flag_e_event = new std::vector<int>;
//   std::vector<int> *om_number = new std::vector<int>;
//   std::vector<double> *energyvis_ubc = new std::vector<double>;
//   std::vector<double> *charge = new std::vector<double>;
//   std::vector<int> *flag_charge = new std::vector<int>;
//   std::vector<int> *flag_associated_nohit = new std::vector<int>;
//   std::vector<int> *flag_MW = new std::vector<int>;
//   std::vector<int> *flag_calo_square = new std::vector<int>;
//   std::vector<int> *flag_source_square = new std::vector<int>;
//   std::vector<int> *flag_last_column = new std::vector<int>;
//   std::vector<int> *flag_last_z = new std::vector<int>;
//
//   TFile *file = new TFile("Simu/cut_bi.root", "READ");
//   TTree* tree = (TTree*)file->Get("Result_tree");
//   tree->SetBranchStatus("*",0);
//   tree->SetBranchStatus("ncalo_tot",1);
//   tree->SetBranchAddress("ncalo_tot", &ncalo_tot);
//   tree->SetBranchStatus("om_num",1);
//   tree->SetBranchAddress("om_num", &om_number);
//   tree->SetBranchStatus("flag_e_event",1);
//   tree->SetBranchAddress("flag_e_event", &flag_e_event);
//   tree->SetBranchStatus("flag_charge",1);
//   tree->SetBranchAddress("flag_charge", &flag_charge);
//   tree->SetBranchStatus("flag_associated_nohit",1);
//   tree->SetBranchAddress("flag_associated_nohit", &flag_associated_nohit);
//   tree->SetBranchStatus("flag_MW",1);
//   tree->SetBranchAddress("flag_MW", &flag_MW);
//   tree->SetBranchStatus("flag_calo_square",1);
//   tree->SetBranchAddress("flag_calo_square", &flag_calo_square);
//   tree->SetBranchStatus("flag_source_square",1);
//   tree->SetBranchAddress("flag_source_square", &flag_source_square);
//   tree->SetBranchStatus("flag_last_column",1);
//   tree->SetBranchAddress("flag_last_column", &flag_last_column);
//   tree->SetBranchStatus("flag_last_z",1);
//   tree->SetBranchAddress("flag_last_z", &flag_last_z);
//   tree->SetBranchStatus("energyvis_ubc",1);
//   tree->SetBranchAddress("energyvis_ubc", &energyvis_ubc);
//
//
//   double rate_simu[128][712];
//   memset (rate_simu, 0, 120*712*sizeof(double));
//
//   for (int i = 0; i < tree->GetEntries()/10; i++) {
//     tree->GetEntry(i);
//     for (int j = 0; j < ncalo_tot; j++) {
//       for (size_t k = 0; k < 2; k++) {
//         for (size_t l = 0; l < 2; l++) {
//           for (size_t m = 0; m < 2; m++) {
//             for (size_t n = 0; n < 2; n++) {
//               for (size_t o = 0; o < 2; o++) {
//                 for (size_t p = 0; p < 2; p++) {
//                   for (size_t q = 0; q < 2; q++) {
//                     if (ncalo_tot < 2 && flag_charge->at(j) == p && flag_MW->at(j) == q && flag_associated_nohit->at(j) == k && flag_calo_square->at(j) == l && flag_source_square->at(j) == m && flag_last_column->at(j) == n && flag_last_z->at(j) == o) {
//                       rate_simu[64*k+32*l+16*m+8*n+4*o+2*p+q][om_number->at(j)]++;
//                     }
//                     // cout << "numero :" << 64*k+32*l+16*m+8*n+4*o+2*p+q << endl;
//                   }
//                 }
//               }
//             }
//           }
//         }
//       }
//     }
//   }
//
//   TFile *file2 = new TFile("cut_974.root", "READ");
//   TTree* tree2 = (TTree*)file2->Get("Result_tree");
//   tree2->SetBranchStatus("*",0);
//   tree2->SetBranchStatus("ncalo_tot",1);
//   tree2->SetBranchAddress("ncalo_tot", &ncalo_tot);
//   tree2->SetBranchStatus("om_number",1);
//   tree2->SetBranchAddress("om_number", &om_number);
//   tree2->SetBranchStatus("flag_e_event",1);
//   tree2->SetBranchAddress("flag_e_event", &flag_e_event);
//   tree2->SetBranchStatus("flag_charge",1);
//   tree2->SetBranchAddress("flag_charge", &flag_charge);
//   tree2->SetBranchStatus("flag_associated_nohit",1);
//   tree2->SetBranchAddress("flag_associated_nohit", &flag_associated_nohit);
//   tree2->SetBranchStatus("flag_MW",1);
//   tree2->SetBranchAddress("flag_MW", &flag_MW);
//   tree2->SetBranchStatus("flag_calo_square",1);
//   tree2->SetBranchAddress("flag_calo_square", &flag_calo_square);
//   tree2->SetBranchStatus("flag_source_square",1);
//   tree2->SetBranchAddress("flag_source_square", &flag_source_square);
//   tree2->SetBranchStatus("flag_last_column",1);
//   tree2->SetBranchAddress("flag_last_column", &flag_last_column);
//   tree2->SetBranchStatus("flag_last_z",1);
//   tree2->SetBranchAddress("flag_last_z", &flag_last_z);
//   tree2->SetBranchStatus("charge",1);
//   tree2->SetBranchAddress("charge", &charge);
//
//   double rate[128][712];
//   memset (rate, 0, 120*712*sizeof(double));
//
//   for (int i = 0; i < tree->GetEntries()/10; i++) {
//     tree->GetEntry(i);
//     for (int j = 0; j < ncalo_tot; j++) {
//       for (size_t k = 0; k < 2; k++) {
//         for (size_t l = 0; l < 2; l++) {
//           for (size_t m = 0; m < 2; m++) {
//             for (size_t n = 0; n < 2; n++) {
//               for (size_t o = 0; o < 2; o++) {
//                 for (size_t p = 0; p < 2; p++) {
//                   for (size_t q = 0; q < 2; q++) {
//                     if (ncalo_tot < 2 && flag_charge->at(j) == p && flag_MW->at(j) == q && flag_associated_nohit->at(j) == k && flag_calo_square->at(j) == l && flag_source_square->at(j) == m && flag_last_column->at(j) == n && flag_last_z->at(j) == o) {
//                       rate[64*k+32*l+16*m+8*n+4*o+2*p+q][om_number->at(j)]++;
//                     }
//                     // cout << "numero :" << 64*k+32*l+16*m+8*n+4*o+2*p+q << endl;
//                   }
//                 }
//               }
//             }
//           }
//         }
//       }
//     }
//   }
//
//
//   TFile *newfile = new TFile("flag_comparator.root", "RECREATE");
//   newfile->cd();
//
//   TH3D map_it3("map_it", "map_it",20,0,20,13,0,13,128,0,128);
//   TH3D map_fr3("map_fr", "map_fr",20,0,20,13,0,13,128,0,128);
//
//   for (int j = 0; j < 128; j++) {
//     TH2D map_it(Form("map_it_%d", j), Form("map_it_%d", j),20,0,20,13,0,13);
//     TH2D map_fr(Form("map_fr_%d", j), Form("map_fr_%d", j),20,0,20,13,0,13);
//     for (int i = 0; i < 260; i++) {
//       cout << "flag : " << j  << "  om : " << i << "  |rate simu " <<  rate_simu[j][i] << " and rate " << rate[j][i] << endl;
//       map_it.SetBinContent(i/13 +1, i%13 +1, rate_simu[j][i]/(rate[j][i]/1.78));
//       map_fr.SetBinContent(i/13 +1, i%13 +1, rate_simu[j][i+260]/(rate[j][i+260]/1.78));
//       map_it3.SetBinContent(i/13 +1, i%13 +1,j, rate_simu[j][i]/(rate[j][i]/1.78));
//       map_fr3.SetBinContent(i/13 +1, i%13 +1,j, rate_simu[j][i+260]/(rate[j][i+260]/1.78));
//     }
//     map_it.Write();
//     map_fr.Write();
//
//
//   }
//     map_it3.Write();
//     map_fr3.Write();
//
//
//   // for (int i = 0; i < 520; i++) {
//   //   cout << "om : " << i << " -> " << rate_simu[i] << " : " << rate[i] << endl;
//   // }
//
//
//
//
//   newfile->Close();
// }
