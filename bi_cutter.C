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

int calo_track(vector<int> track_column, vector<int> track_layer, int calo_column){        /// comparison between the calo column and the tracker layer
  int compteur = 0;
  if (calo_column == 0 && compteur == 0)
    for (int i = 0; i < track_column.size(); i++)
      if(track_column.at(i) >= 0 && track_column.at(i) <= 5 && track_layer.at(i) >= 6 && track_layer.at(i) < 8) compteur++;
  if (calo_column == 1 && compteur == 0)
    for (int i = 0; i < track_column.size(); i++)
      if(track_column.at(i) >= 4 && track_column.at(i) <= 10 && track_layer.at(i) >= 6 && track_layer.at(i) < 8) compteur++;
  if (calo_column == 2 && compteur == 0)
    for (int i = 0; i < track_column.size(); i++)
      if(track_column.at(i) >= 8 && track_column.at(i) <= 16 && track_layer.at(i) >= 6 && track_layer.at(i) < 8) compteur++;
  if (calo_column == 3 && compteur == 0)
    for (int i = 0; i < track_column.size(); i++)
      if(track_column.at(i) >= 14 && track_column.at(i) <= 22 && track_layer.at(i) >= 6 && track_layer.at(i) < 8) compteur++;
  if (calo_column == 4 && compteur == 0)
    for (int i = 0; i < track_column.size(); i++)
      if(track_column.at(i) >= 21 && track_column.at(i) <= 28 && track_layer.at(i) >= 6 && track_layer.at(i) < 8) compteur++;
  if (calo_column == 5 && compteur == 0)
    for (int i = 0; i < track_column.size(); i++)
      if(track_column.at(i) >= 26 && track_column.at(i) <= 33 && track_layer.at(i) >= 6 && track_layer.at(i) < 8) compteur++;
  if (calo_column == 6 && compteur == 0)
    for (int i = 0; i < track_column.size(); i++)
      if(track_column.at(i) >= 31 && track_column.at(i) <= 39 && track_layer.at(i) >= 6 && track_layer.at(i) < 8) compteur++;
  if (calo_column == 7 && compteur == 0)
    for (int i = 0; i < track_column.size(); i++)
      if(track_column.at(i) >= 38 && track_column.at(i) <= 45 && track_layer.at(i) >= 6 && track_layer.at(i) < 8) compteur++;
  if (calo_column == 8 && compteur == 0)
    for (int i = 0; i < track_column.size(); i++)
      if(track_column.at(i) >= 44 && track_column.at(i) <= 51 && track_layer.at(i) >= 6 && track_layer.at(i) < 8) compteur++;
  if (calo_column == 9 && compteur == 0)
    for (int i = 0; i < track_column.size(); i++)
      if(track_column.at(i) >= 49 && track_column.at(i) <= 57 && track_layer.at(i) >= 6 && track_layer.at(i) < 8) compteur++;
  if (calo_column == 10 && compteur == 0)
    for (int i = 0; i < track_column.size(); i++)
      if(track_column.at(i) >= 56 && track_column.at(i) <= 63 && track_layer.at(i) >= 6 && track_layer.at(i) < 8) compteur++;
  if (calo_column == 11 && compteur == 0)
    for (int i = 0; i < track_column.size(); i++)
      if(track_column.at(i) >= 62 && track_column.at(i) <= 67 && track_layer.at(i) >= 6 && track_layer.at(i) < 8) compteur++;
  if (calo_column == 12 && compteur == 0)
    for (int i = 0; i < track_column.size(); i++)
      if(track_column.at(i) >= 67 && track_column.at(i) <= 73 && track_layer.at(i) >= 6 && track_layer.at(i) < 8) compteur++;
  if (calo_column == 13 && compteur == 0)
    for (int i = 0; i < track_column.size(); i++)
      if(track_column.at(i) >= 72 && track_column.at(i) <= 80 && track_layer.at(i) >= 6 && track_layer.at(i) < 8) compteur++;
  if (calo_column == 14 && compteur == 0)
    for (int i = 0; i < track_column.size(); i++)
      if(track_column.at(i) >= 79 && track_column.at(i) <= 86 && track_layer.at(i) >= 6 && track_layer.at(i) < 8) compteur++;
  if (calo_column == 15 && compteur == 0)
    for (int i = 0; i < track_column.size(); i++)
      if(track_column.at(i) >= 85 && track_column.at(i) <= 92 && track_layer.at(i) >= 6 && track_layer.at(i) < 8) compteur++;
  if (calo_column == 16 && compteur == 0)
    for (int i = 0; i < track_column.size(); i++)
      if(track_column.at(i) >= 91 && track_column.at(i) <= 98 && track_layer.at(i) >= 6 && track_layer.at(i) < 8) compteur++;
  if (calo_column == 17 && compteur == 0)
    for (int i = 0; i < track_column.size(); i++)
      if(track_column.at(i) >= 97 && track_column.at(i) <= 104 && track_layer.at(i) >= 6 && track_layer.at(i) < 8) compteur++;
  if (calo_column == 18 && compteur == 0)
    for (int i = 0; i < track_column.size(); i++)
      if(track_column.at(i) >= 102 && track_column.at(i) <= 109 && track_layer.at(i) >= 6 && track_layer.at(i) < 8) compteur++;
  if (calo_column == 19 && compteur == 0)
    for (int i = 0; i < track_column.size(); i++)
      if(track_column.at(i) >= 108 && track_column.at(i) <= 112 && track_layer.at(i) >= 6 && track_layer.at(i) < 8) compteur++;

  if (compteur > 3) return 1;
  else return 0;
}

int Bi_source_checker(int track_layer){        /// comparison between the calo column and the tracker layer
  if (track_layer >= 6 && track_layer <= 11) return 1;
  if (track_layer >= 25 && track_layer <= 30) return 1;
  if (track_layer >= 44 && track_layer <= 49) return 1;
  if (track_layer >= 63 && track_layer <= 68) return 1;
  if (track_layer >= 82 && track_layer <= 87) return 1;
  if (track_layer >= 101 && track_layer <= 106) return 1;
  else return 0;
}

int z_comparator(double R0, double R5, double R6, int om_number) {
  double z_calo = (om_number % 13) * 0.256 - 0.164 + 0.128;            // om_number % 13 = row number, 0.256 = OM height. to start at the same height than the 3m tracker, substract 0.164 ((0.256*13-3)/2) + 0.128 to be on the middle of the OM
  double z_track;
  double L_cell = 2.9;
  double t_top = (R6 - R0)*12.5E-3;                                 // Put it in µs
  double t_bottom = (R5 - R0)*12.5E-3;
  double prop_rate = L_cell/(65/*anode_tdc_r1*12.5E-3 + anode_tdc_r2*12.5E-3*/);

  if (t_bottom > 0 && t_bottom < 1.1E+17 && t_top > 0 && t_top < 1.1E+17) z_track =  ((t_bottom - t_top)/(t_top + t_bottom))*L_cell/2;
  else if ((t_bottom <= 0 || t_bottom > 1.1E+17) && (t_top > 0 && t_top < 1.1E+17)) {z_track = 3 - t_top*prop_rate;}
  else if ((t_top <= 0 || t_top > 1.1E+17) && (t_bottom > 0 && t_bottom < 1.1E+17)) {z_track = t_bottom*prop_rate;}
  else if ((t_bottom <= 0 && t_top <= 0) || (t_bottom > 1.1E+17 && t_top > 1.1E+17) || (t_bottom <= 0 && t_top > 1.1E+17) || (t_bottom  > 1.1E+17 && t_top <= 0)) {z_track = 1000;}


  if (z_track < -L_cell/2) z_track = -L_cell/2;
  if (z_track > L_cell/2 && z_track != 1000) z_track = L_cell/2;


  if (z_track <= z_calo + 0.256 && z_track >= z_calo - 256) return 1;
  else return 0;
}

void calo_track_selectionner(vector<int> track_side, vector<int> track_layer, vector<int> track_column, int compteur_calo[][20], vector<int> &last_column, vector<int> &last_column_side, vector<vector<long>> R0, double timestamp, vector<int> &cut_track_side, vector<int> &cut_track_layer , vector<int> &cut_track_column){        /// comparison between the calo column and the tracker layer
  for (int i = 0; i < track_side.size(); i++) {
    if ((2*R0.at(i).at(0) - timestamp)*6.25e-9 > (-0.2e-6) && (2*R0.at(i).at(0) - timestamp)*6.25e-9 < (5e-6)) {
      cut_track_side.push_back(track_side.at(i));
      cut_track_layer.push_back(track_layer.at(i));
      cut_track_column.push_back(track_column.at(i));
      if (track_layer.at(i) > 4) {
        if (track_column.at(i) >= 0 && track_column.at(i) <= 5) compteur_calo[track_side.at(i)][0]++;
        if (track_column.at(i) >= 4 && track_column.at(i) <= 10) compteur_calo[track_side.at(i)][1]++;
        if (track_column.at(i) >= 8 && track_column.at(i) <= 16) compteur_calo[track_side.at(i)][2]++;
        if (track_column.at(i) >= 14 && track_column.at(i) <= 22) compteur_calo[track_side.at(i)][3]++;
        if (track_column.at(i) >= 21 && track_column.at(i) <= 28) compteur_calo[track_side.at(i)][4]++;
        if (track_column.at(i) >= 26 && track_column.at(i) <= 33) compteur_calo[track_side.at(i)][5]++;
        if (track_column.at(i) >= 31 && track_column.at(i) <= 39) compteur_calo[track_side.at(i)][6]++;
        if (track_column.at(i) >= 38 && track_column.at(i) <= 45) compteur_calo[track_side.at(i)][7]++;
        if (track_column.at(i) >= 44 && track_column.at(i) <= 51) compteur_calo[track_side.at(i)][8]++;
        if (track_column.at(i) >= 49 && track_column.at(i) <= 57) compteur_calo[track_side.at(i)][9]++;
        if (track_column.at(i) >= 56 && track_column.at(i) <= 63) compteur_calo[track_side.at(i)][10]++;
        if (track_column.at(i) >= 62 && track_column.at(i) <= 67) compteur_calo[track_side.at(i)][11]++;
        if (track_column.at(i) >= 67 && track_column.at(i) <= 73) compteur_calo[track_side.at(i)][12]++;
        if (track_column.at(i) >= 72 && track_column.at(i) <= 80) compteur_calo[track_side.at(i)][13]++;
        if (track_column.at(i) >= 79 && track_column.at(i) <= 86) compteur_calo[track_side.at(i)][14]++;
        if (track_column.at(i) >= 85 && track_column.at(i) <= 92) compteur_calo[track_side.at(i)][15]++;
        if (track_column.at(i) >= 91 && track_column.at(i) <= 98) compteur_calo[track_side.at(i)][16]++;
        if (track_column.at(i) >= 97 && track_column.at(i) <= 104) compteur_calo[track_side.at(i)][17]++;
        if (track_column.at(i) >= 102 && track_column.at(i) <= 109) compteur_calo[track_side.at(i)][18]++;
        if (track_column.at(i) >= 108 && track_column.at(i) <= 112) compteur_calo[track_side.at(i)][19]++;
        if (track_layer.at(i) == 8) {
          last_column.push_back(track_column.at(i));
          last_column_side.push_back(track_side.at(i));
        }
      }
    }
  }
}

void source_track_selectionner(vector<int> track_side, vector<int> track_layer, vector<int> track_column, int compteur_source[][6], vector<vector<long>> R0, double timestamp) {
  for (int i = 0; i < track_layer.size(); i++) {
    if ((2*R0.at(i).at(0) - timestamp)*6.25e-9 > (-0.2e-6) && (2*R0.at(i).at(0) - timestamp)*6.25e-9 < (5e-6)) {
      if (track_layer.at(i) < 4) {
        for (int j = 0; j < 6; j++) {
          if (abs(track_column.at(i) - (8.5 + j * 19)) < 3) {
            compteur_source[track_side.at(i)][j]++;
          }
        }
      }
    }
  }
}

int calo_source_track(int compteur_source[][6], int compteur_calo[][20], int calo_column, int calo_side){        /// source localiser
  if (compteur_calo[calo_side][calo_column] < 3) return 0;
  for (int i = 0; i < 20; i++) {
    if (compteur_source[calo_side][i] > 3) return 1;
  }
  return 0;
}

// void backtoback(){
//
//
//
//
// }

void track_cutter() {
  TFile *file = new TFile("snemo_run-902_udd.root", "READ");
  std::vector<vector<short>> *wave = new std::vector<vector<short>>;
  std::vector<vector<long>> *R0 = new std::vector<vector<long>>;
  std::vector<vector<long>> *R5 = new std::vector<vector<long>>;
  std::vector<vector<long>> *R6 = new std::vector<vector<long>>;
  std::vector<long> *timestamp = new std::vector<long>;
  std::vector<int> *calo_wall = new std::vector<int>;
  std::vector<int> *calo_side = new std::vector<int>;
  std::vector<int> *calo_column = new std::vector<int>;
  std::vector<int> *calo_row = new std::vector<int>;
  std::vector<int> *calo_charge = new std::vector<int>;
  std::vector<int> *calo_ampl = new std::vector<int>;
  std::vector<int> *tracker_side = new std::vector<int>;
  std::vector<int> *tracker_column = new std::vector<int>;
  std::vector<int> *tracker_layer = new std::vector<int>;
  std::vector<int> *last_column = new std::vector<int>;
  std::vector<int> *last_column_side = new std::vector<int>;
  std::vector<int> *cut_track_side = new std::vector<int>;
  std::vector<int> *cut_track_column = new std::vector<int>;
  std::vector<int> *cut_track_layer = new std::vector<int>;

  int eventnumber, calo_nohits, tracker_nohits;

  TTree* tree = (TTree*)file->Get("SimData");
  tree->SetBranchStatus("*",0);
  tree->SetBranchStatus("header.eventnumber",1);
  tree->SetBranchAddress("header.eventnumber", &eventnumber);
  tree->SetBranchStatus("digicalo.nohits",1);
  tree->SetBranchAddress("digicalo.nohits", &calo_nohits);
  tree->SetBranchStatus("digicalo.timestamp",1);
  tree->SetBranchAddress("digicalo.timestamp", &timestamp);
  tree->SetBranchStatus("digicalo.waveform",1);
  tree->SetBranchAddress("digicalo.waveform", &wave);
  tree->SetBranchStatus("digicalo.wall",1);
  tree->SetBranchAddress("digicalo.wall", &calo_wall);
  tree->SetBranchStatus("digicalo.side",1);
  tree->SetBranchAddress("digicalo.side", &calo_side);
  tree->SetBranchStatus("digicalo.column",1);
  tree->SetBranchAddress("digicalo.column", &calo_column);
  tree->SetBranchStatus("digicalo.row",1);
  tree->SetBranchAddress("digicalo.row", &calo_row);
  tree->SetBranchStatus("digicalo.charge",1);
  tree->SetBranchAddress("digicalo.charge", &calo_charge);
  tree->SetBranchStatus("digicalo.peakamplitude",1);
  tree->SetBranchAddress("digicalo.peakamplitude", &calo_ampl);
  tree->SetBranchStatus("digitracker.nohits",1);
  tree->SetBranchAddress("digitracker.nohits", &tracker_nohits);
  tree->SetBranchStatus("digitracker.side",1);
  tree->SetBranchAddress("digitracker.side", &tracker_side);
  tree->SetBranchStatus("digitracker.layer",1);
  tree->SetBranchAddress("digitracker.layer", &tracker_layer);
  tree->SetBranchStatus("digitracker.column",1);
  tree->SetBranchAddress("digitracker.column", &tracker_column);
  tree->SetBranchStatus("digitracker.anodetimestampR0",1);
  tree->SetBranchAddress("digitracker.anodetimestampR0", &R0);
  tree->SetBranchStatus("digitracker.bottomcathodetimestamp",1);
  tree->SetBranchAddress("digitracker.bottomcathodetimestamp", &R5);
  tree->SetBranchStatus("digitracker.topcathodetimestamp",1);
  tree->SetBranchAddress("digitracker.topcathodetimestamp", &R6);

  TH2D spectre_ampl_full("spectre_ampl_full","spectre_ampl_full", 520, 0, 520, 1000, 0, 20000);
  TH2D spectre_charge_full("spectre_charge_full","spectre_charge_full", 520, 0, 520, 1000, 0, 60000);
  TH2D spectre_ampl("spectre_ampl","spectre_ampl", 520, 0, 520, 1000, 0, 20000);
  TH2D spectre_charge("spectre_charge","spectre_charge", 520, 0, 520, 1000, 0, 60000);
  TH2D spectre_ampl_single("spectre_ampl_single","spectre_ampl_single", 520, 0, 520, 1000, 0, 20000);
  TH2D spectre_charge_single("spectre_charge_single","spectre_charge_single", 520, 0, 520, 1000, 0, 60000);
  int compteur = 0;
  double tracker_hit_side;
  double tracker_hit_layer;
  double tracker_hit_column;

  int compteur_source[2][6] = {0};
  int compteur_calo[2][20] = {0};
  int om_num;
  TFile *newfile = new TFile("test.root", "RECREATE");

  int calo_compteur_nohits = 0;
  int tracker_compteur_nohits = 0;
  int event, calo_hit;
  TTree Result_tree("Result_tree","");
  Result_tree.Branch("eventnumber", &event);
  Result_tree.Branch("om_num", &om_num);
  Result_tree.Branch("calo_nohits", &calo_compteur_nohits);
  Result_tree.Branch("calo_hit", &calo_hit);
  Result_tree.Branch("tracker_side", &cut_track_side);
  Result_tree.Branch("tracker_layer", &cut_track_layer);
  Result_tree.Branch("tracker_column", &cut_track_column);
  Result_tree.Branch("tracker_nohits", &tracker_compteur_nohits);
  Result_tree.Branch("last_column", &last_column);
  Result_tree.Branch("last_column_side", &last_column_side);

  for (int i = 0; i < tree->GetEntries(); i++) {      //loop on event number
  // for (int i = 62; i < 63; i++) {      //loop on event number
    tree->GetEntry(i);
    event = i;
    // if (calo_nohits < 3) {

      for (int k = 0; k < calo_nohits; k++) {      //k : loop on calo hit number
        om_num = calo_side->at(k)*260 + calo_column->at(k)*13 + calo_row->at(k);
        memset(compteur_source, 0, sizeof(int) * 2 * 6);
        memset(compteur_calo, 0, sizeof(int) * 2 * 20);
        source_track_selectionner(*tracker_side, *tracker_layer, *tracker_column, compteur_source, *R0, timestamp->at(k));
        calo_track_selectionner(*tracker_side, *tracker_layer, *tracker_column, compteur_calo, *last_column, *last_column_side, *R0, timestamp->at(k), *cut_track_side, *cut_track_layer, *cut_track_column);
        calo_hit = k;

        // if (-calo_ampl->at(k) > 200 && om_num < 520 && om_num % 13 != 0 && om_num % 13 != 12 && tracker_nohits > 5 && tracker_nohits < 16) {      // condition to cut small charge and keep only MW OM
        if (-calo_ampl->at(k) > 200 && om_num < 520 && om_num % 13 != 0 && om_num % 13 != 12) {      // condition to cut small charge and keep only MW OM
          calo_compteur_nohits++;
          spectre_ampl_full.Fill(om_num, -calo_ampl->at(k));            //TH2 on amplitude and charge spectra without cut
          spectre_charge_full.Fill(om_num, -calo_charge->at(k));
          // last_column_selecteur(*last_column, calo_column->at(k));

          if (calo_source_track(compteur_source, compteur_calo, calo_column->at(k), calo_side->at(k)) == 1 && last_column->size() > 0){
            for (int l = 0; l < tracker_side->size(); l++) {
              if (tracker_side->at(l) == calo_side->at(k)) {
                tracker_compteur_nohits++;
              }
              else {
                tracker_side->at(l) = -1;
                tracker_layer->at(l) = -1;
                tracker_column->at(l) = -1;
              }
            }
            spectre_ampl.Fill(om_num, -calo_ampl->at(k));           //TH2 on amplitude and charge spectra with cut
            spectre_charge.Fill(om_num, -calo_charge->at(k));
            Result_tree.Fill();

            if (calo_nohits < 3) {
              spectre_ampl_single.Fill(om_num, -calo_ampl->at(k));
              spectre_charge_single.Fill(om_num, -calo_charge->at(k));
            }
          }
          // if (calo_side->at(k) == 1) return;
        }
        cut_track_side->clear();
        cut_track_layer->clear();
        cut_track_column->clear();
        last_column->clear();
        last_column_side->clear();
      }
    // }
    calo_compteur_nohits = 0;
    tracker_compteur_nohits = 0;

  }


  newfile->cd();
  Result_tree.Write();
  spectre_ampl.Write();
  spectre_charge.Write();
  spectre_ampl_single.Write();
  spectre_charge_single.Write();
  spectre_ampl_full.Write();
  spectre_charge_full.Write();
  newfile->Close();
  cout << "OK" << endl;
}

// void new_track_cutter() {
//   TFile *file = new TFile("snemo_run-902_udd.root", "READ");
//   std::vector<vector<long>> *R0 = new std::vector<vector<long>>;
//   std::vector<vector<long>> *R5 = new std::vector<vector<long>>;
//   std::vector<vector<long>> *R6 = new std::vector<vector<long>>;
//   std::vector<int> *calo_wall = new std::vector<int>;
//   std::vector<int> *calo_side = new std::vector<int>;
//   std::vector<int> *calo_column = new std::vector<int>;
//   std::vector<int> *calo_row = new std::vector<int>;
//   std::vector<int> *calo_charge = new std::vector<int>;
//   std::vector<int> *calo_ampl = new std::vector<int>;
//   std::vector<int> *tracker_side = new std::vector<int>;
//   std::vector<int> *tracker_column = new std::vector<int>;
//   std::vector<int> *tracker_layer = new std::vector<int>;
//
//   std::vector<int> *tracker_side_buffer = new std::vector<int>;
//   std::vector<int> *tracker_column_buffer = new std::vector<int>;
//   std::vector<int> *tracker_layer_buffer = new std::vector<int>;
//
//   int eventnumber, calo_nohits, tracker_nohits;
//   double backtoback[2] = {0,0};
//   TTree* tree = (TTree*)file->Get("SimData");
//   tree->SetBranchStatus("*",0);
//   tree->SetBranchStatus("header.eventnumber",1);
//   tree->SetBranchAddress("header.eventnumber", &eventnumber);
//   tree->SetBranchStatus("digicalo.nohits",1);
//   tree->SetBranchAddress("digicalo.nohits", &calo_nohits);
//   tree->SetBranchStatus("digicalo.wall",1);
//   tree->SetBranchAddress("digicalo.wall", &calo_wall);
//   tree->SetBranchStatus("digicalo.side",1);
//   tree->SetBranchAddress("digicalo.side", &calo_side);
//   tree->SetBranchStatus("digicalo.column",1);
//   tree->SetBranchAddress("digicalo.column", &calo_column);
//   tree->SetBranchStatus("digicalo.row",1);
//   tree->SetBranchAddress("digicalo.row", &calo_row);
//   tree->SetBranchStatus("digicalo.charge",1);
//   tree->SetBranchAddress("digicalo.charge", &calo_charge);
//   tree->SetBranchStatus("digicalo.peakamplitude",1);
//   tree->SetBranchAddress("digicalo.peakamplitude", &calo_ampl);
//   tree->SetBranchStatus("digitracker.nohits",1);
//   tree->SetBranchAddress("digitracker.nohits", &tracker_nohits);
//   tree->SetBranchStatus("digitracker.side",1);
//   tree->SetBranchAddress("digitracker.side", &tracker_side);
//   tree->SetBranchStatus("digitracker.layer",1);
//   tree->SetBranchAddress("digitracker.layer", &tracker_layer);
//   tree->SetBranchStatus("digitracker.column",1);
//   tree->SetBranchAddress("digitracker.column", &tracker_column);
//   tree->SetBranchStatus("digitracker.anodetimestampR0",1);
//   tree->SetBranchAddress("digitracker.anodetimestampR0", &R0);
//   tree->SetBranchStatus("digitracker.bottomcathodetimestamp",1);
//   tree->SetBranchAddress("digitracker.bottomcathodetimestamp", &R5);
//   tree->SetBranchStatus("digitracker.topcathodetimestamp",1);
//   tree->SetBranchAddress("digitracker.topcathodetimestamp", &R6);
//
//   TH2D spectre_ampl_full("spectre_ampl_full","spectre_ampl_full", 520, 0, 520, 1000, 0, 20000);
//   TH2D spectre_charge_full("spectre_charge_full","spectre_charge_full", 520, 0, 520, 1000, 0, 60000);
//   TH2D spectre_ampl("spectre_ampl","spectre_ampl", 520, 0, 520, 1000, 0, 20000);
//   TH2D spectre_charge("spectre_charge","spectre_charge", 520, 0, 520, 1000, 0, 60000);
//   TH2D spectre_ampl_backtoback("spectre_amp_backtobackl","spectre_ampl_backtoback", 520, 0, 520, 1000, 0, 20000);
//   TH2D spectre_charge_backtoback("spectre_charge_backtoback","spectre_charge_backtoback", 520, 0, 520, 1000, 0, 60000);
//   int compteur = 0;
//   double calo_hit, event;
//   double tracker_hit_side;
//   double tracker_hit_layer;
//   double tracker_hit_column;
//
//   TFile *newfile = new TFile("new_track.root", "RECREATE");
//
//   TTree Result_tree("Result_tree","");
//   Result_tree.Branch("eventnumber", &eventnumber);
//   Result_tree.Branch("calo_hit", &calo_hit);
//   Result_tree.Branch("tracker_hit_side", &tracker_hit_side);
//   Result_tree.Branch("tracker_hit_layer", &tracker_hit_layer);
//   Result_tree.Branch("tracker_hit_column", &tracker_hit_column);
//
//
//   std::vector<int> *charge_buffer = new std::vector<int>;
//   std::vector<int> *ampl_buffer = new std::vector<int>;
//   std::vector<int> *om = new std::vector<int>;
//
//   for (int i = 0; i < tree->GetEntries(); i++) {      //loop on event number
//     tree->GetEntry(i);
//
//     for (int k = 0; k < calo_nohits; k++) {      //k : loop on calo hit number
//       int om_num = calo_side->at(k)*260 + calo_column->at(k)*13 + calo_row->at(k);
//       calo_hit = om_num; // for sndisplay
//
//       if (-calo_ampl->at(k) > 200 && om_num < 520 && om_num % 13 != 0 && om_num % 13 != 12) {      // condition to cut small charge and keep only MW OM
//
//         spectre_ampl_full.Fill(om_num, -calo_ampl->at(k));            //TH2 on amplitude and charge spectra without cut
//         spectre_charge_full.Fill(om_num, -calo_charge->at(k));
//
//         for (int j = 0; j < tracker_nohits; j++) {        //j  loop on tracker hit number
//           if (tracker_side->at(j) == calo_side->at(k)) {
//             tracker_layer_buffer->push_back(tracker_layer->at(j));
//             tracker_column_buffer->push_back(tracker_column->at(j));
//           }
//         }
//
//         if(source_track(*tracker_column_buffer, *tracker_layer_buffer, calo_column->at(k)) == 1){
//           if (calo_track(*tracker_column_buffer, *tracker_layer_buffer, calo_column->at(k)) == 1) {
//             if (calo_side->at(k) == 0) backtoback[0] = 1;
//             if (calo_side->at(k) == 1) backtoback[1] = 1;
//             charge_buffer->push_back(-calo_charge->at(k));
//             ampl_buffer->push_back(-calo_ampl->at(k));
//             spectre_ampl.Fill(om_num, -calo_ampl->at(k));           //TH2 on amplitude and charge spectra with cut
//             spectre_charge.Fill(om_num, -calo_charge->at(k));
//
//             om->push_back(om_num);
//                 //
//             for (int l = 0; l < tracker_nohits; l++) {
//               tracker_hit_side = tracker_side->at(l);
//               tracker_hit_layer = tracker_layer->at(l);
//               tracker_hit_column = tracker_column->at(l);
//               Result_tree.Fill();
//             }
//
//             compteur++;
//           }
//         }
//       }
//     }
//
//     // if (backtoback[0] == 1 && backtoback[1] == 1) {
//     //   for (int f = 0; f < charge_buffer->size(); f++) {
//     //     spectre_ampl_backtoback.Fill(om->at(f), ampl_buffer->at(f));           //TH2 on amplitude and charge spectra with cut
//     //     spectre_charge_backtoback.Fill(om->at(f), charge_buffer->at(f));
//     //   }
//     //         // for sndisplay
//     //   for (int k = 0; k < calo_nohits; k++) {      //k : loop on calo hit number
//     //     for (int l = 0; l < tracker_nohits; l++) {
//     //       int om_num = calo_side->at(k)*260 + calo_column->at(k)*13 + calo_row->at(k);
//     //       calo_hit = om_num; // for sndisplay
//     //       tracker_hit_side = tracker_side->at(l);
//     //       tracker_hit_layer = tracker_layer->at(l);
//     //       tracker_hit_column = tracker_column->at(l);
//     //       Result_tree.Fill();
//     //     }
//     //   }
//     // }
//     // backtoback[0] = 0;
//     // backtoback[1] = 0;
//
//     // if(compteur == 1) sndisplay_demonstrator(calo_hit, tracker_hit, i);
//     charge_buffer->clear();
//     ampl_buffer->clear();
//     om->clear();
//
//   }
//
//
//   newfile->cd();
//   Result_tree.Write();
//   spectre_ampl.Write();
//   spectre_charge.Write();
//   spectre_ampl_backtoback.Write();
//   spectre_charge_backtoback.Write();
//   spectre_ampl_full.Write();
//   spectre_charge_full.Write();
//   newfile->Close();
//   cout << "OK" << endl;
// }
//

void backtoback() {
  TFile *file = new TFile("test.root", "READ");
  double calo_hit, tracker_hit_side, tracker_hit_layer, tracker_hit_column;
  int eventnumber;

  TTree* tree = (TTree*)file->Get("Result_tree");
  tree->SetBranchStatus("*",0);
  tree->SetBranchStatus("calo_hit",1);
  tree->SetBranchAddress("calo_hit", &calo_hit);
  tree->SetBranchStatus("tracker_hit_side",1);
  tree->SetBranchAddress("tracker_hit_side", &tracker_hit_side);
  tree->SetBranchStatus("tracker_hit_side",1);
  tree->SetBranchAddress("tracker_hit_side", &tracker_hit_side);
  tree->SetBranchStatus("tracker_hit_column",1);
  tree->SetBranchAddress("tracker_hit_column", &tracker_hit_column);
  tree->SetBranchStatus("eventnumber",1);
  tree->SetBranchAddress("eventnumber", &eventnumber);

  TH2D spectre_ampl_backtoback("spectre_amp_backtobackl","spectre_ampl_backtoback", 520, 0, 520, 1000, 0, 20000);
  TH2D spectre_charge_backtoback("spectre_charge_backtoback","spectre_charge_backtoback", 520, 0, 520, 1000, 0, 60000);

  TFile *newfile = new TFile("test_back.root", "RECREATE");

  TTree Result_tree("Result_tree","");
  Result_tree.Branch("eventnumber", &eventnumber);
  Result_tree.Branch("calo_hit", &calo_hit);
  Result_tree.Branch("tracker_hit_side", &tracker_hit_side);
  Result_tree.Branch("tracker_hit_layer", &tracker_hit_layer);
  Result_tree.Branch("tracker_hit_column", &tracker_hit_column);

  std::vector<int> *event_vec = new std::vector<int>;

  double prev_calo_hit = -1;
  int prev_eventnumber = -1;

  for (int i = 0; i < 200; i++) {
    tree->GetEntry(i);
    if (prev_calo_hit != calo_hit && prev_eventnumber == eventnumber) {
      event_vec->push_back(eventnumber);
      // cout << "ok" << endl;
    }
    // cout << "om = " << calo_hit << " and prev = " << prev_calo_hit << "  <<   ev n = " << eventnumber  << " and prev = " << prev_eventnumber << endl;
    prev_calo_hit = calo_hit;
    prev_eventnumber = eventnumber;
  }
  // return;
  int j = 0;
  int compteur = 0;
  for (int i = 0; i < 200; i++) {
    tree->GetEntry(i);

    cout <<"ev n = " << eventnumber  << " and prev = " << event_vec->at(j) << endl;
    if (eventnumber == event_vec->at(j)) {
      compteur = 1;
      Result_tree.Fill();
    }
    else if(compteur == 1) {
      j++;
      i-=1;
    }
    compteur = 0;
  }
  newfile->cd();
  Result_tree.Write();
  newfile->Close();
}

int main() {
  track_cutter();
  return 0;
}

// if ((calo_column == 0 || calo_column == 1 || calo_column == 2) && compteur == 0)
//   for (int i = 0; i < track_column.size(); i++)
//     if(track_column.at(i) >= 6 && track_column.at(i) <= 11 && track_layer.at(i) >= 0 && track_layer.at(i) < 3) compteur++;
// if (calo_column == 3 && compteur == 0)
//   for (int i = 0; i < track_column.size(); i++)
//     if(((track_column.at(i) >= 6 && track_column.at(i) <= 11) || (track_column.at(i) >= 25 && track_column.at(i) <= 30)) && track_layer.at(i) >= 0 && track_layer.at(i) < 3) compteur++;
// if ((calo_column == 4 ||calo_column == 5) && compteur == 0)
//   for (int i = 0; i < track_column.size(); i++)
//     if(track_column.at(i) >= 25 && track_column.at(i) <= 30 && track_layer.at(i) >= 0 && track_layer.at(i) < 3) compteur++;
// if (calo_column == 6 && compteur == 0)
//   for (int i = 0; i < track_column.size(); i++)
//     if(((track_column.at(i) >= 44 && track_column.at(i) <= 49) || (track_column.at(i) >= 25 && track_column.at(i) <= 30)) && track_layer.at(i) >= 0 && track_layer.at(i) < 3) compteur++;
// if ((calo_column == 7 ||calo_column == 8) && compteur == 0)
//   for (int i = 0; i < track_column.size(); i++)
//     if(track_column.at(i) >= 44 && track_column.at(i) <= 49 && track_layer.at(i) >= 0 && track_layer.at(i) < 3) compteur++;
// if ((calo_column == 9 || calo_column == 10) && compteur == 0)
//   for (int i = 0; i < track_column.size(); i++)
//     if(((track_column.at(i) >= 44 && track_column.at(i) <= 49) || (track_column.at(i) >= 63 && track_column.at(i) <= 68)) && track_layer.at(i) >= 0 && track_layer.at(i) < 3) compteur++;
// if ((calo_column == 11 ||calo_column == 12) && compteur == 0)
//   for (int i = 0; i < track_column.size(); i++)
//     if(track_column.at(i) >= 63 && track_column.at(i) <= 68 && track_layer.at(i) >= 0 && track_layer.at(i) < 3) compteur++;
// if (calo_column == 13 && compteur == 0)
//   for (int i = 0; i < track_column.size(); i++)
//     if(((track_column.at(i) >= 82 && track_column.at(i) <= 87) || (track_column.at(i) >= 63 && track_column.at(i) <= 68)) && track_layer.at(i) >= 0 && track_layer.at(i) < 3) compteur++;
// if ((calo_column == 14 || calo_column == 15) && compteur == 0)
//   for (int i = 0; i < track_column.size(); i++)
//     if(track_column.at(i) >= 82 && track_column.at(i) <= 87 && track_layer.at(i) >= 0 && track_layer.at(i) < 3) compteur++;
// if (calo_column == 16 && compteur == 0)
//   for (int i = 0; i < track_column.size(); i++)
//     if(((track_column.at(i) >= 82 && track_column.at(i) <= 87) || (track_column.at(i) >= 101 && track_column.at(i) <= 106)) && track_layer.at(i) >= 0 && track_layer.at(i) < 3) compteur++;
// if ((calo_column == 17 || calo_column == 18 || calo_column == 19) && compteur == 0)
//   for (int i = 0; i < track_column.size(); i++)
//     if(track_column.at(i) >= 101 && track_column.at(i) <= 106 && track_layer.at(i) >= 0 && track_layer.at(i) < 3) compteur++;
//
// if (compteur > 3) return 1;
// else return 0;
// }

// int calo_source_track(int compteur_source[][6], int compteur_calo[][20], int calo_column, int calo_side){        /// source localiser
//   if (compteur_calo[calo_side][calo_column] < 3) return 0;
//   if (calo_column == 0 && compteur_source[calo_side][0] > 2) return 1;
//   if (calo_column == 1 && compteur_source[calo_side][0] > 2) return 1;
//   if (calo_column == 2 && compteur_source[calo_side][0] > 2) return 1;
//   if (calo_column == 3 && (compteur_source[calo_side][0] > 2 || compteur_source[calo_side][1] > 2)) return 1;
//   if (calo_column == 4 && compteur_source[calo_side][1] > 2) return 1;
//   if (calo_column == 5 && compteur_source[calo_side][1] > 2) return 1;
//   if (calo_column == 6 && (compteur_source[calo_side][1] > 2 || compteur_source[calo_side][2] > 2)) return 1;
//   if (calo_column == 7 && compteur_source[calo_side][2] > 2) return 1;
//   if (calo_column == 8 && compteur_source[calo_side][2] > 2) return 1;
//   if (calo_column == 9 && (compteur_source[calo_side][2] > 2 || compteur_source[calo_side][3] > 2)) return 1;
//   if (calo_column == 10 && (compteur_source[calo_side][2] > 2 || compteur_source[calo_side][3] > 2)) return 1;
//   if (calo_column == 11 && compteur_source[calo_side][3] > 2) return 1;
//   if (calo_column == 12 && compteur_source[calo_side][3] > 2) return 1;
//   if (calo_column == 13 && (compteur_source[calo_side][3] > 2 || compteur_source[calo_side][4] > 2)) return 1;
//   if (calo_column == 14 && compteur_source[calo_side][4] > 2) return 1;
//   if (calo_column == 15 && compteur_source[calo_side][4] > 2) return 1;
//   if (calo_column == 16 && (compteur_source[calo_side][4] > 2 || compteur_source[calo_side][5] > 2)) return 1;
//   if (calo_column == 17 && compteur_source[calo_side][5] > 2) return 1;
//   if (calo_column == 18 && compteur_source[calo_side][5] > 2) return 1;
//   if (calo_column == 19 && compteur_source[calo_side][5] > 2) return 1;
//   else return 0;
// }
// void sndisplay_demonstrator(std::vector<short> calo_hit, std::vector<vector<short>> tracker_hit, int hit_number) {
//   sndisplay::demonstrator *sndemonstrator = new sndisplay::demonstrator ("demonstrator_test");
//
//   for (int i = 0; i < calo_hit.size(); i++) {
//     sndemonstrator->setomcontent(calo_hit.at(i),  1.0);
//   }
//
//   for (int i = 0; i < tracker_hit.at(0).size(); i++) {
//     sndemonstrator->setggcontent(tracker_hit.at(0).at(i), tracker_hit.at(1).at(i), tracker_hit.at(2).at(i), 1);
//   }
//
//   sndemonstrator->settitle("RUN 609 // TRIGGER 9");
//   sndemonstrator->draw_top();
//
//   // disable cells in area != 0
//   for (int side=0; side<2; ++side)
//     for (int row=15; row<113; ++row)
//       for (int layer=0; layer<9; ++layer)
// 	sndemonstrator->setggcolor(side, row, layer, kGray+1);
//
//   sndemonstrator->canvas->SaveAs(Form("events/event_n%d.png", hit_number));
//
//
// }
