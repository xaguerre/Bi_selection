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

double z_calculator_gg(double R0, double R5, double R6){
  double z_gg;
  double t_top = (R6 - R0)*12.5E-3;                                 // Put it in µs
  double t_bottom = (R5 - R0)*12.5E-3;


  if (t_bottom > 0 && t_bottom < 1.1E+17 && t_top > 0 && t_top < 1.1E+17) z_gg =  ((t_bottom - t_top)/(t_top + t_bottom));

  else z_gg = -9999;
  return z_gg;
}

void calo_track_selectionner(vector<int> track_side, vector<int> track_layer, vector<int> track_column, int compteur_calo[][20], vector<vector<long>> R0, vector<vector<long>> R5, vector<vector<long>> R6, double timestamp, double *z_column){        /// comparison between the calo column and the tracker layer
  vector<long> last_R5;
  vector<long> last_R6;
  vector<long> last_R0;
  vector<double> last_column;
  vector<int> last_side;
  vector<double> z;
  vector<double> col;

  int tab_bas[20] = {0,4,8,14,21,26,31,38,44,49,56,62,67,72,79,85,91,97,102,108};
  int tab_haut[20] = {5,10,16,22,28,33,39,45,51,57,63,67,73,80,86,92,98,104,109,112};


  for (int i = 0; i < track_side.size(); i++) {
    if ((2*R0.at(i).at(0) - timestamp)*6.25e-9 > (-0.2e-6) && (2*R0.at(i).at(0) - timestamp)*6.25e-9 < (5e-6)) {
      if (track_layer.at(i) > 4) {
        for (int j = 0; j < 20; j++) {
          if (track_column.at(i) >= tab_bas[j] && track_column.at(i) <= tab_haut[j]) compteur_calo[track_side.at(i)][j]++;
        }
        if (track_layer.at(i) == 8) {
          last_column.push_back(track_column.at(i));
          last_side.push_back(track_side.at(i));
          last_R5.push_back(R5.at(i).at(0));
          last_R6.push_back(R6.at(i).at(0));
          last_R0.push_back(R0.at(i).at(0));
        }
      }
    }
  }

  for (int i = 0; i < last_column.size(); i++) {
    for (int j = 0; j < 20; j++) {
      if (last_column.at(i) >= tab_bas[j] && last_column.at(i) <= tab_haut[j] && compteur_calo[last_side.at(i)][j] > 3){
        // cout << "R0 = " << last_R0.at(i) << "  and R5 = " << last_R5.at(i) << "  and R6 = " << last_R6.at(i) << endl;
        z.push_back(z_calculator_gg(last_R0.at(i), last_R5.at(i), last_R6.at(i)));
        // cout << "z = " << z.back() << endl;

        col.push_back(last_column.at(i));
      }
    }
  }

  double column;
  double z_gg;
  if (z.size() > 0) {
    if (z.size() > 1) {
      for (int i = 0; i < z.size(); i++) {
          // cout << "z = " << z << endl;
        z_gg += z.at(i);
        column += col.at(i);
      }
      z_column[0] = z_gg / z.size();
      z_column[1] = column / z.size();

    }
    else{
      z_column[0] = z.at(0);
      z_column[1] = col.at(0);
    }
  }
  else {
    z_column[0] = -1000;
    z_column[1] = -1000;
  }
}

void source_track_selectionner(vector<int> track_side, vector<int> track_layer, vector<int> track_column, int compteur_source[][6], vector<vector<long>> R0, vector<vector<long>> R5, vector<vector<long>> R6, double timestamp, double *z_column) {
  vector<long> first_R5;
  vector<long> first_R6;
  vector<long> first_R0;
  vector<double> first_column;
  vector<int> first_side;
  vector<double> z;
  vector<double> col;

  for (int i = 0; i < track_layer.size(); i++) {
    if ((2*R0.at(i).at(0) - timestamp)*6.25e-9 > (-0.2e-6) && (2*R0.at(i).at(0) - timestamp)*6.25e-9 < (5e-6)) {
      if (track_layer.at(i) < 4) {
        for (int j = 0; j < 6; j++) {
          if (abs(track_column.at(i) - (8.5 + j * 19)) < 3) {
            compteur_source[track_side.at(i)][j]++;
            if (track_layer.at(i) == 0) {
              first_side.push_back(track_side.at(i));
              first_column.push_back(track_column.at(i));
              first_R5.push_back(R5.at(i).at(0));
              first_R6.push_back(R6.at(i).at(0));
              first_R0.push_back(R0.at(i).at(0));
            }
          }
        }
      }
    }
  }

  for (int i = 0; i < first_column.size(); i++) {
    for (int j = 0; j < 6; j++) {
      if (abs(first_column.at(i)- (8.5 + j * 19)) < 3  && compteur_source[first_side.at(i)][j] > 3){
        z.push_back(z_calculator_gg(first_R0.at(i), first_R5.at(i), first_R6.at(i)));
        col.push_back(first_column.at(i));
      }
    }
  }

  double column;
  double z_gg;
  if (z.size() > 0) {
    if (z.size() > 1) {
      for (int i = 0; i < z.size(); i++) {
          // cout << "z = " << z << endl;
        z_gg += z.at(i);
        column += col.at(i);
      }
      z_column[0] = z_gg / z.size();
      z_column[1] = column / z.size();

    }
    else{
      z_column[0] = z.at(0);
      z_column[1] = col.at(0);
    }
  }
  else {
    z_column[0] = -1000;
    z_column[1] = -1000;
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

void track_cutter() {
  TFile *file = new TFile("data/snemo_run-902_udd.root", "READ");
  std::vector<vector<short>> *wave = new std::vector<vector<short>>;
  std::vector<vector<long>> *R0 = new std::vector<vector<long>>;
  std::vector<vector<long>> *R5 = new std::vector<vector<long>>;
  std::vector<vector<long>> *R6 = new std::vector<vector<long>>;
  std::vector<long> *timestamp = new std::vector<long>;
  std::vector<int> *calo_type = new std::vector<int>;
  std::vector<int> *calo_side = new std::vector<int>;
  std::vector<int> *calo_wall = new std::vector<int>;
  std::vector<int> *calo_column = new std::vector<int>;
  std::vector<int> *calo_row = new std::vector<int>;
  std::vector<int> *calo_charge = new std::vector<int>;
  std::vector<int> *calo_ampl = new std::vector<int>;
  std::vector<int> *tracker_side = new std::vector<int>;
  std::vector<int> *tracker_column = new std::vector<int>;
  std::vector<int> *tracker_layer = new std::vector<int>;
  std::vector<int> *cut_track_side = new std::vector<int>;
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
  tree->SetBranchStatus("digicalo.wall",1);
  tree->SetBranchAddress("digicalo.wall", &calo_wall);
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

  int compteur_source[2][6] = {0};
  int compteur_calo[2][20] = {0};
  double* z_column_last = new double[2];
  double* z_column_first = new double[2];
  int om_num, nelec, ngamma;
  TFile *newfile = new TFile("cutted_bi.root", "RECREATE");

  int calo_compteur_nohits = 0;
  double z_last_gg, z_first_gg;
  std::vector<int> om_number;
  std::vector<int> charge;
  std::vector<int> amplitude;
  std::vector<int> calo_tdc;
  std::vector<int> associated_track;
  std::vector<int> source;
  std::vector<int> vec_nelec;
  std::vector<int> vec_ngamma;
  std::vector<int> e_event;
  std::vector<double> vec_z_last_gg;
  std::vector<double> last_column;
  std::vector<double> vec_z_first_gg;
  std::vector<double> first_column;

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
  Result_tree.Branch("z_last_gg", &vec_z_last_gg);
  Result_tree.Branch("last_column", &last_column);
  Result_tree.Branch("z_first_gg", &vec_z_first_gg);
  Result_tree.Branch("first_column", &first_column);
  Result_tree.Branch("emitting_source", &source);

  for (int i = 0; i < tree->GetEntries(); i++) {      //loop on event number
  // for (int i = 0; i < 1000; i++) {      //loop on event number
  if (i % 100000 == 0) {
    cout << i << endl;
  }
    tree->GetEntry(i);
    z_column_last[0] = 0;
    z_column_last[1] = 0;
    z_column_first[0] = 0;
    z_column_first[1] = 0;
    // cout << "eventnumber = " << eventnumber << endl;
    int validator = 0;
    for (int k = 0; k < calo_nohits; k++) {      //k : loop on calo hit number
      // cout << "calo hit = " << k << endl;
      if (calo_type->at(k) == 0) om_num = calo_side->at(k)*260 + calo_column->at(k)*13 + calo_row->at(k);
      if (calo_type->at(k) == 1) om_num = 520 + calo_side->at(k)*64 +  calo_wall->at(k)*32 + calo_column->at(k)*16 + calo_row->at(k);
      if (calo_type->at(k) == 2) om_num = 520 + 128 + calo_side->at(k)*32 + calo_wall->at(k)*16 + calo_column->at(k);

      if (-calo_ampl->at(k) > 200) {      // condition to cut small charge and keep only MW OM
        int timed_gg;
        memset(compteur_source, 0, sizeof(int) * 2 * 6);
        memset(compteur_calo, 0, sizeof(int) * 2 * 20);

        source_track_selectionner(*tracker_side, *tracker_layer, *tracker_column, compteur_source, *R0, *R5, *R6, timestamp->at(k), z_column_first);
        calo_track_selectionner(*tracker_side, *tracker_layer, *tracker_column, compteur_calo, *R0, *R5, *R6, timestamp->at(k), z_column_last);

        z_last_gg = z_column_last[0];
        z_first_gg = z_column_first[0];
        timed_gg = gg_counter(timestamp->at(k), *R0, tracker_nohits);
        calo_compteur_nohits++;
        if (timed_gg > 5 && timed_gg < 16 && om_num < 520 && om_num % 13 != 0 && om_num % 13 != 12){
          if (calo_source_track(compteur_source, compteur_calo, calo_column->at(k), calo_side->at(k)) == 1 && z_last_gg <= 0.0926 + 0.1852*((om_num%13)-6)  && z_last_gg >= -0.093 + 0.1852*((om_num%13)-6)) {       //// Delta z moyen = 0.268589cm, soit 0.1852 de -1 à 1
            // cout << "side =" << calo_side->at(k) << endl;
            vec_z_last_gg.push_back(z_last_gg);
            vec_z_first_gg.push_back(z_first_gg);
            // cout << "z = " << z_last_gg << endl;
            vec_nelec.push_back(1);
            last_column.push_back(z_column_last[1]);
            first_column.push_back(z_column_first[1]);
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
    last_column.clear();
    vec_z_last_gg.clear();
    first_column.clear();
    vec_z_first_gg.clear();
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
