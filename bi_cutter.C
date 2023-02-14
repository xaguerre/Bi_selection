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

int gg_counter(double timestamp, vector<vector<long>> R0, int ngg){
  int gg_counter = 0;
  for (int i = 0; i < R0.size(); i++) {
    if ((2*R0.at(i).at(0) - timestamp)*6.25e-9 > (-0.2e-6) && (2*R0.at(i).at(0) - timestamp)*6.25e-9 < (5e-6)) {
      gg_counter++;
    }
  }
  return gg_counter;
}


// void backtoback() {
//   TFile *file = new TFile("test.root", "READ");
//   double charge;
//   int eventnumber, om_number, calo_nohits;
//   vector<int> tracker_side;
//   vector<int> tracker_layer;
//   vector<int> tracker_column;
//
//   TTree* tree = (TTree*)file->Get("Result_tree");
//   tree->SetBranchStatus("*",0);
//   tree->SetBranchStatus("om_number",1);
//   tree->SetBranchAddress("om_number", &om_number);
//   tree->SetBranchStatus("calo_nohits",1);
//   tree->SetBranchAddress("calo_nohits", &calo_nohits);
//   tree->SetBranchStatus("charge",1);
//   tree->SetBranchAddress("charge", &charge);
//   tree->SetBranchStatus("tracker_side",1);
//   tree->SetBranchAddress("tracker_side", &tracker_side);
//   tree->SetBranchStatus("tracker_side",1);
//   tree->SetBranchAddress("tracker_side", &tracker_side);
//   tree->SetBranchStatus("tracker_column",1);
//   tree->SetBranchAddress("tracker_column", &tracker_column);
//   tree->SetBranchStatus("eventnumber",1);
//   tree->SetBranchAddress("eventnumber", &eventnumber);
//
//   TH2D spectre_ampl_backtoback("spectre_amp_backtobackl","spectre_ampl_backtoback", 520, 0, 520, 1000, 0, 20000);
//   TH2D spectre_charge_backtoback("spectre_charge_backtoback","spectre_charge_backtoback", 520, 0, 520, 1000, 0, 60000);
//
//   TFile *newfile = new TFile("test_back.root", "RECREATE");
//
//   TTree Result_tree("Result_tree","");
//   Result_tree.Branch("eventnumber", &eventnumber);
//   Result_tree.Branch("om_number", &om_number);
//   Result_tree.Branch("charge", &charge);
//   Result_tree.Branch("tracker_hit_side", &tracker_hit_side);
//   Result_tree.Branch("tracker_hit_layer", &tracker_hit_layer);
//   Result_tree.Branch("tracker_hit_column", &tracker_hit_column);
//
//   std::vector<int> *event_vec = new std::vector<int>;
//
//   double prev_calo_hit = -1;
//   int prev_eventnumber = -1;
//
//   for (int i = 0; i < 200; i++) {
//     tree->GetEntry(i);
//     if (prev_calo_hit != calo_hit && prev_eventnumber == eventnumber) {
//       event_vec->push_back(eventnumber);
//       // cout << "ok" << endl;
//     }
//     // cout << "om = " << calo_hit << " and prev = " << prev_calo_hit << "  <<   ev n = " << eventnumber  << " and prev = " << prev_eventnumber << endl;
//     prev_calo_hit = calo_hit;
//     prev_eventnumber = eventnumber;
//   }
//   // return;
//   int j = 0;
//   int compteur = 0;
//   for (int i = 0; i < 200; i++) {
//     tree->GetEntry(i);
//
//     cout <<"ev n = " << eventnumber  << " and prev = " << event_vec->at(j) << endl;
//     if (eventnumber == event_vec->at(j)) {
//       compteur = 1;
//       Result_tree.Fill();
//     }
//     else if(compteur == 1) {
//       j++;
//       i-=1;
//     }
//     compteur = 0;
//   }
//   newfile->cd();
//   Result_tree.Write();
//   newfile->Close();
// }

void track_cutter() {
  TFile *file = new TFile("data/snemo_run-902_udd.root", "READ");
  std::vector<vector<short>> *wave = new std::vector<vector<short>>;
  std::vector<vector<long>> *R0 = new std::vector<vector<long>>;
  std::vector<vector<long>> *R5 = new std::vector<vector<long>>;
  std::vector<vector<long>> *R6 = new std::vector<vector<long>>;
  std::vector<long> *timestamp = new std::vector<long>;
  std::vector<int> *calo_type = new std::vector<int>;
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
  tree->SetBranchStatus("digicalo.type",1);
  tree->SetBranchAddress("digicalo.type", &calo_type);
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

  int compteur = 0;
  double tracker_hit_side;
  double tracker_hit_layer;
  double tracker_hit_column;

  int compteur_source[2][6] = {0};
  int compteur_calo[2][20] = {0};
  int om_num, nelec, ngamma;
  TFile *newfile = new TFile("test2.root", "RECREATE");

  int calo_compteur_nohits = 0;
  int event = 0;

  std::vector<int> om_number;
  std::vector<int> charge;
  std::vector<int> amplitude;
  std::vector<int> calo_tdc;
  std::vector<int> associated_track;
  std::vector<int> source;
  std::vector<int> vec_nelec;
  std::vector<int> vec_ngamma;
  std::vector<int> e_event;

  TTree Result_tree("Result_tree","");
  Result_tree.Branch("eventnumber", &eventnumber);
  Result_tree.Branch("ncalo_tot", &calo_compteur_nohits);
  Result_tree.Branch("ngg_tot", &tracker_nohits);
  Result_tree.Branch("nelec", &nelec);
  Result_tree.Branch("ngamma", &ngamma);
  Result_tree.Branch("e_event", &e_event);
  Result_tree.Branch("om_number", &om_number);
  Result_tree.Branch("charge", &charge);
  Result_tree.Branch("amplitude", &amplitude);
  Result_tree.Branch("calo_tdc", &calo_tdc);
  Result_tree.Branch("gg_associated", &associated_track);
  Result_tree.Branch("emitting_source", &source);



  // for (int i = 0; i < tree->GetEntries(); i++) {      //loop on event number
  for (int i = 30; i < 31; i++) {      //loop on event number
    tree->GetEntry(i);
    int validator = 0;
    for (int k = 0; k < calo_nohits; k++) {      //k : loop on calo hit number
      if (calo_type.at(k) == 0) {
        om_num = calo_side->at(k)*260 + calo_column->at(k)*13 + calo_row->at(k);
        cout << om_num << endl;
        if (-calo_ampl->at(k) > 200 && om_num < 520 && om_num % 13 != 0 && om_num % 13 != 12) {      // condition to cut small charge and keep only MW OM
          int timed_gg;

          memset(compteur_source, 0, sizeof(int) * 2 * 6);
          memset(compteur_calo, 0, sizeof(int) * 2 * 20);
          source_track_selectionner(*tracker_side, *tracker_layer, *tracker_column, compteur_source, *R0, timestamp->at(k));
          calo_track_selectionner(*tracker_side, *tracker_layer, *tracker_column, compteur_calo, *last_column, *last_column_side, *R0, timestamp->at(k), *cut_track_side, *cut_track_layer, *cut_track_column);

          timed_gg = gg_counter(timestamp->at(k), *R0, tracker_nohits);
          calo_compteur_nohits++;
          if (timed_gg > 5 && timed_gg < 16){
            if (calo_source_track(compteur_source, compteur_calo, calo_column->at(k), calo_side->at(k)) == 1 && last_column->size() > 0) {
              vec_nelec.push_back(1);
              e_event.push_back(1);
              associated_track.push_back(timed_gg);
              validator = 1;
            }
            else {
              vec_ngamma.push_back(1);
              e_event.push_back(0);
              associated_track.push_back(0);
            }
          }
          else {
            vec_ngamma.push_back(1);
            e_event.push_back(0);
            associated_track.push_back(0);
          }
          om_number.push_back(om_num);
          charge.push_back(-calo_charge->at(k));
          amplitude.push_back(-calo_ampl->at(k));
          calo_tdc.push_back(timestamp->at(k));

          // source.push_back
        }
      }
    }
    nelec = 0;
    ngamma = 0;
    if (validator == 1) {
      for (int j = 0; j < vec_nelec.size(); j++) {nelec+=vec_nelec.at(j);}
      for (int l = 0; l < vec_ngamma.size(); l++) {ngamma+=vec_ngamma.at(l);}
      Result_tree.Fill();
      // event++;
    }
    validator = 0;
    vec_nelec.clear();
    vec_ngamma.clear();
    om_number.clear();
    charge.clear();
    amplitude.clear();
    calo_tdc.clear();
    associated_track.clear();
    source.clear();
    e_event.clear();
    calo_compteur_nohits = 0;
  }


  newfile->cd();
  Result_tree.Write();
  newfile->Close();
  cout << "OK" << endl;
}

int main() {
  track_cutter();
  return 0;
}
