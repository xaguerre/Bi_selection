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

double source_column_row_mean[6][7];
double source_column_row_sigma[6][7];

double z_calculator_gg(double R0, double R5, double R6){
  double z_gg;
  double t_top = (R6 - R0)*12.5E-3;                                 // Put it in µs
  double t_bottom = (R5 - R0)*12.5E-3;

  // cout << "tb =" << t_bottom << " and tt = " << t_top << endl;
  if (t_bottom > 0 && t_bottom < 1.1E+17 && t_top > 0 && t_top < 1.1E+17) z_gg = ((t_bottom - t_top)/(t_top + t_bottom));
  else z_gg = -9999;
  return z_gg;
}

int row_to_source(int column){
  int colonne_gauche[5] = {4, 22, 42, 60, 78};
  int colonne_droite[5] = {13, 33, 52, 72, 91};

  for (int i = 0; i < 5; i++) {
    if (column > colonne_gauche[i] && column < colonne_droite[i]) {
      return i;
    }
  }
  return 10000;
}

void first_z_selectionner() {

  double mean, sigma;

  TFile *file = new TFile("first_z/z_distrib.root", "READ");
  TTree* tree = (TTree*)file->Get("Result_tree");
  tree->SetBranchStatus("*",0);
  tree->SetBranchStatus("mean",1);
  tree->SetBranchAddress("mean", &mean);
  tree->SetBranchStatus("sigma",1);
  tree->SetBranchAddress("sigma", &sigma);

  int compteur = 0;
  for(int i = 0; i < 4; i++){                   //loop on source column
    for(int j = 1; j < 6; j++){                 //loop on source row
      tree->GetEntry(compteur);
      source_column_row_mean[i][j] = mean;
      source_column_row_sigma[i][j] = sigma;
      // cout << "source " << i*100 + j << " mean = " << mean << " +- " << 3*sigma << endl;
      compteur++;
    }
  }
}

int source_numberer(double z, int column){
  // if (z < -0.5 && z > -1) cout << "z = " << z << " and column = " << column << endl;
  int source_number = -10;
  if (column > -1 && column < 6) {
  for(int j = 1; j < 6; j++){ // loop on the source row
      // if (z < -0.5 && z > -1)cout << " j = " << j << " and min z = " << source_column_row_mean[column][j] - 3*source_column_row_sigma[column][j] << " and max z = " << source_column_row_mean[column][j] + 3*source_column_row_sigma[column][j] << endl;
      if (z > source_column_row_mean[column][j] - 3*source_column_row_sigma[column][j] && z < source_column_row_mean[column][j] + 3*source_column_row_sigma[column][j]) {
        source_number = column*100 + j;
        return source_number;
      }
      else source_number = -9999;
    }
  }
  return source_number;
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


  z_column[0] = -1000;
  z_column[1] = -1000;

  double column;
  double z_gg;
  if (z.size() > 0) {
    if (z.size() > 1) {      // if two cell on the las layer in front of the OM
      for (int i = 0; i < z.size(); i++) {
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
}

void calo_track_selectionner(vector<int> track_side, vector<int> track_layer, vector<int> track_column, int compteur_calo[][20], vector<vector<long>> R0, vector<vector<long>> R5, vector<vector<long>> R6, double timestamp, double *z_column, int *flag6){        /// comparison between the calo column and the tracker layer
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
          *flag6 = 1;
          last_column.push_back(track_column.at(i));
          last_side.push_back(track_side.at(i));
          last_R5.push_back(R5.at(i).at(0));
          last_R6.push_back(R6.at(i).at(0));
          last_R0.push_back(R0.at(i).at(0));

        }
      }
    }
  }

  if (*flag6 == 9) *flag6 = 0;
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

  z_column[0] = -1000;
  z_column[1] = -1000;


  double column;
  double z_gg;
  if (z.size() > 0) {
    if (z.size() > 1) {
      for (int i = 0; i < z.size(); i++) {
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
}

int calo_source_track(int compteur_source[][6], int compteur_calo[][20], int calo_column, int calo_side, int* flag4, int *flag5){        /// source localiser
  if (compteur_calo[calo_side][calo_column] >= 3) *flag4 = 1;
  else *flag4 = 0;
  if (compteur_calo[calo_side][calo_column] >= 3) {
    for (int i = 0; i < 20; i++) {
      if (compteur_source[calo_side][i] >= 3){
        *flag5 = 1;
        return 1;
      }
      else *flag5 = 0;
    }
  }
  else *flag5 = 0;
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
  first_z_selectionner();
  // return;
  TFile *file = new TFile("data/snemo_run-974_udd.root", "READ");
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

  std::vector<int> flag;
  std::vector<int> flag0;
  std::vector<int> flag1;
  std::vector<int> flag2;
  std::vector<int> flag3;
  std::vector<int> flag4;
  std::vector<int> flag5;
  std::vector<int> flag6;
  std::vector<int> flag7;

  int compteur_source[2][6] = {0};
  int compteur_calo[2][20] = {0};
  double* z_column_last = new double[2];
  double* z_column_first = new double[2];
  int om_num, nelec, ngamma;

  TFile *newfile = new TFile("cut_974.root", "RECREATE");
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
  Result_tree.Branch("om_number", &om_number);
  Result_tree.Branch("charge", &charge);
  Result_tree.Branch("amplitude", &amplitude);
  Result_tree.Branch("calo_tdc", &calo_tdc);
  Result_tree.Branch("gg_associated", &associated_track);
  Result_tree.Branch("z_last_gg", &vec_z_last_gg);
  Result_tree.Branch("last_column", &last_column);
  Result_tree.Branch("z_first_gg", &vec_z_first_gg);
  Result_tree.Branch("first_column", &first_column);
  Result_tree.Branch("source_number", &source);
  Result_tree.Branch("flag", &flag);
  Result_tree.Branch("flag_e_event", &flag0);
  Result_tree.Branch("flag_charge", &flag1);
  Result_tree.Branch("flag_associated_nohit", &flag2);
  Result_tree.Branch("flag_MW", &flag3);
  Result_tree.Branch("flag_calo_square", &flag4);
  Result_tree.Branch("flag_source_square", &flag5);
  Result_tree.Branch("flag_last_column", &flag6);
  Result_tree.Branch("flag_last_z", &flag7);

  for (int i = 0; i < tree->GetEntries(); i++) {      //loop on event number
  // for (int i = 0; i < 3000; i++) {      //loop on event number

  if (i % 100000 == 0) {
    cout << 100.*i/tree->GetEntries() << "%" << endl;
  }
    tree->GetEntry(i);
    z_column_last[0] = 0;
    z_column_last[1] = 0;
    z_column_first[0] = 0;
    z_column_first[1] = 0;
    int validator = 0;
    for (int k = 0; k < calo_nohits; k++) {      //k : loop on calo hit number
      int flag4_int = 0;
      int flag5_int = 0;
      int flag6_int = 0;
      flag0.push_back(9);
      flag1.push_back(9);
      flag2.push_back(9);
      flag3.push_back(9);
      flag4.push_back(9);
      flag5.push_back(9);
      flag6.push_back(9);
      flag7.push_back(9);
      if (calo_type->at(k) == 0) om_num = calo_side->at(k)*260 + calo_column->at(k)*13 + calo_row->at(k);
      if (calo_type->at(k) == 1) om_num = 520 + calo_side->at(k)*64 +  calo_wall->at(k)*32 + calo_column->at(k)*16 + calo_row->at(k);
      if (calo_type->at(k) == 2) om_num = 520 + 128 + calo_side->at(k)*32 + calo_wall->at(k)*16 + calo_column->at(k);
      om_number.push_back(om_num);
      charge.push_back(-calo_charge->at(k));
      amplitude.push_back(-calo_ampl->at(k));
      calo_tdc.push_back(timestamp->at(k));


      if (-calo_ampl->at(k) > 200) {      // condition to cut small charge and keep only MW OM
        flag1.pop_back();
        flag1.push_back(1);
        int timed_gg;
        memset(compteur_source, 0, sizeof(int) * 2 * 6);
        memset(compteur_calo, 0, sizeof(int) * 2 * 20);

        source_track_selectionner(*tracker_side, *tracker_layer, *tracker_column, compteur_source, *R0, *R5, *R6, timestamp->at(k), z_column_first);
        calo_track_selectionner(*tracker_side, *tracker_layer, *tracker_column, compteur_calo, *R0, *R5, *R6, timestamp->at(k), z_column_last, &flag6_int);

        flag6.pop_back();
        flag6.push_back(flag6_int);
        z_last_gg = z_column_last[0];
        z_first_gg = z_column_first[0];
        timed_gg = gg_counter(timestamp->at(k), *R0, tracker_nohits);
        calo_compteur_nohits++;
        int first_gg_row = row_to_source(z_column_first[1]);
        int source_number = source_numberer(z_column_first[0], first_gg_row);
        // if (z_column_first[0] < 1 && z_column_first[0] > -1) cout << "z = " << z_column_first[0] << " and row  = "<< z_column_first[1] << " and row = " << first_gg_row << endl;
        // if (z_column_first[0] < 1 && z_column_first[0] > -1) cout << "-> source : " << source_number << endl;
        // cout << "calo_hit = " << om_num << endl;
        // cout << "calo source track = " << calo_source_track(compteur_source, compteur_calo, calo_column->at(k), calo_side->at(k)) << endl;
        // cout << "zlastgg = " << z_last_gg << endl;
        // cout << "" << endl;
        if (timed_gg > 5 && timed_gg < 16) {
          flag2.pop_back();
          flag2.push_back(1);
        }
        else if (flag2.at(k) == 9) {
          flag2.pop_back();
          flag2.push_back(0);
        }
        if (om_num < 520 && om_num % 13 != 0 && om_num % 13 != 12) {
          flag3.pop_back();
          flag3.push_back(1);
        }
        else if (flag3.at(k) == 9) {
          flag3.pop_back();
          flag3.push_back(0);
        }

        if (timed_gg > 5 && timed_gg < 16 && om_num < 520 && om_num % 13 != 0 && om_num % 13 != 12){
          if (z_last_gg <= 0.1 + 0.2*((om_num%13)-6)  && z_last_gg >= -0.1 + 0.2*((om_num%13)-6)){
            flag7.pop_back();
            flag7.push_back(1);
          }
          else if (flag7.at(k) == 9) {
            flag7.pop_back();
            flag7.push_back(0);
          }

          double calo_row2z = -15;
          double z_last_gg_min = -15;
          double z_last_gg_max = -15;

          if (calo_type->at(k) == 0){
            calo_row2z = -1.1165 +0.18714*calo_row->at(k);
            z_last_gg_min = calo_row2z - 0.1;
            z_last_gg_max = calo_row2z + 0.1;
          }
          if (calo_source_track(compteur_source, compteur_calo, calo_column->at(k), calo_side->at(k), &flag4_int, &flag5_int) == 1 && z_last_gg <= z_last_gg_max && z_last_gg >= z_last_gg_min){//  && source_number >= 0) {       //// Delta z moyen = 0.268589cm 0.29 à 3 sigma, soit 0.1852 (1)de -1 à 1 de moyenneur gategauss
            vec_z_last_gg.push_back(z_last_gg);
            vec_z_first_gg.push_back(z_first_gg);
            // cout << "first col = " << z_column_first[1] << endl;
            vec_nelec.push_back(1);
            last_column.push_back(z_column_last[1]);
            first_column.push_back(z_column_first[1]);
            associated_track.push_back(timed_gg);
            source.push_back(source_number);

            flag0.pop_back();
            flag0.push_back(1);
          }
          else {
            vec_ngamma.push_back(1);
            associated_track.push_back(0);
            source.push_back(-9999);
            flag0.pop_back();
            flag0.push_back(0);
          }
          flag4.pop_back();
          flag4.push_back(flag4_int);
          flag5.pop_back();
          flag5.push_back(flag5_int);
        }
        else {
          vec_ngamma.push_back(1);
          associated_track.push_back(0);
          source.push_back(-9999);
          flag0.pop_back();
          flag0.push_back(0);
        }


        // source.push_back
      }
      else {
        flag1.pop_back();
        flag1.push_back(0);
      }
      flag.push_back(flag0.at(k) +flag1.at(k)*10 +flag2.at(k)*100 + flag3.at(k)*1000 + flag4.at(k)*10000 + flag5.at(k)*100000 + flag6.at(k)*1000000 + flag7.at(k)*10000000);
    }
    nelec = 0;
    ngamma = 0;

    for (int j = 0; j < vec_nelec.size(); j++) {nelec+=vec_nelec.at(j);}
    for (int l = 0; l < vec_ngamma.size(); l++) {ngamma+=vec_ngamma.at(l);}

    Result_tree.Fill();
    flag0.clear();
    flag1.clear();
    flag2.clear();
    flag3.clear();
    flag4.clear();
    flag5.clear();
    flag6.clear();
    flag7.clear();
    flag.clear();
    validator = 0;
    vec_nelec.clear();
    vec_ngamma.clear();
    om_number.clear();
    charge.clear();
    amplitude.clear();
    calo_tdc.clear();
    associated_track.clear();
    source.clear();
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

void tri() {
  gStyle->SetOptFit(1);
  gStyle->SetOptStat(1);


  TFile *file = new TFile("data/snemo_run-974_udd.root", "READ");
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
  std::vector<long> *R0_new = new std::vector<long>;
  std::vector<long> *R5_new = new std::vector<long>;
  std::vector<long> *R6_new = new std::vector<long>;
  std::vector<int> *cell_number = new std::vector<int>;


  int eventnumber, calo_nohits, tracker_nohits;

  TTree* tree = (TTree*)file->Get("SimData");
  tree->SetBranchStatus("*",0);
  tree->SetBranchStatus("header.eventnumber",1);
  tree->SetBranchAddress("header.eventnumber", &eventnumber);
  tree->SetBranchStatus("digicalo.nohits",1);
  tree->SetBranchAddress("digicalo.nohits", &calo_nohits);
  tree->SetBranchStatus("digicalo.timestamp",1);
  tree->SetBranchAddress("digicalo.timestamp", &timestamp);
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


  TFile *newfile = new TFile("data/test_timestamp.root", "RECREATE");
  TTree Result_tree("Result_tree","");
  Result_tree.Branch("eventnumber", &eventnumber);
  Result_tree.Branch("tracker_side", &tracker_side);
  Result_tree.Branch("tracker_layer", &tracker_layer);
  Result_tree.Branch("tracker_column", &tracker_column);
  Result_tree.Branch("R0", &R0_new);
  Result_tree.Branch("R5", &R5_new);
  Result_tree.Branch("R6", &R6_new);
  Result_tree.Branch("cell_number", &cell_number);

  double om_520 [520];
  TH2D* om = new TH2D("om", "om", 520, 0, 520, 1000, 0, 50000);
  TH2D* R0_map = new TH2D("R0_map", "R0_map", 113, 0, 112, 19, 0, 18);
  TH2D* R0R5_map = new TH2D("R0R5_map", "R0R5_map", 113, 0, 112, 19, 0, 18);
  TH2D* R0R6_map = new TH2D("R0R6_map", "R0R6_map", 113, 0, 112, 19, 0, 18);

  for (int i = 0; i < tree->GetEntries(); i++) {      //loop on event number
  // for (int i = 0; i < 3000; i++) {      //loop on event number

    if (i % 100000 == 0) {
      cout << 100.*i/tree->GetEntries() << "%" << endl;
    }
    tree->GetEntry(i);
    // for (int j = 0; j < calo_charge->size(); j++) {
    //   if (calo_type == 0) {
    //     om_520[calo_side->at(j)*260 + calo_column->at(j)*13 + calo_row->at(j)] ++;
    //   }
    // }

    for (int j = 0; j < calo_charge->size(); j++) {
      for (int k = 0; k < tracker_side->size(); k++) {
        if (calo_charge->at(j) > 200) {
          if (R0->at(k).at(0) > 0 ) {
            if ((2*R0->at(k).at(0) - timestamp->at(j))*6.25e-9 > (-0.2e-6) && (2*R0->at(k).at(0) - timestamp->at(j))*6.25e-9 < (5e-6)){
              R0_new->push_back(R0->at(k).at(0));
              R5_new->push_back(R5->at(k).at(0));
              R6_new->push_back(R6->at(k).at(0));
              cell_number->push_back(tracker_side->at(k)*1017 +9*tracker_column->at(k)+tracker_layer->at(k));
              R0_map->Fill(tracker_column->at(k), tracker_side->at(k)*10 + tracker_layer->at(k));
              if (R5->at(k).at(0) > 0) {
                R0R5_map->Fill(tracker_column->at(k), tracker_side->at(k)*10 + tracker_layer->at(k));
              }
              if (R6->at(k).at(0) > 0) {
                R0R6_map->Fill(tracker_column->at(k), tracker_side->at(k)*10 + tracker_layer->at(k));
              }
            }
          }
          Result_tree.Fill();
          R0_new->clear();
          R5_new->clear();
          R6_new->clear();
          cell_number->clear();

        }
      }

    }
  }


  // for (int i = 0; i < 520; i++) {
  //   cout << "spectre om " << i << " -> entry = " << om_520[i] << endl;
  // }


  newfile->cd();
  Result_tree.Write();
  om->Write();
  R0_map->Write();
  R0R5_map->Write();
  R0R6_map->Write();
  newfile->Close();
}
//
void R_calculator() {

TFile *file = new TFile("data/test_timestamp.root", "READ");

  std::vector<int> *tracker_side = new std::vector<int>;
  std::vector<int> *tracker_column = new std::vector<int>;
  std::vector<int> *tracker_layer = new std::vector<int>;
  std::vector<long> *R0 = new std::vector<long>;
  std::vector<long> *R5 = new std::vector<long>;
  std::vector<long> *R6 = new std::vector<long>;
  std::vector<int> *cell_number = new std::vector<int>;

  double R0_new, R5_new, R6_new;
  int cell;

  int eventnumber, calo_nohits, tracker_nohits;

  TTree* tree = (TTree*)file->Get("Result_tree");
  tree->SetBranchStatus("*",0);
  tree->SetBranchStatus("eventnumber",1);
  tree->SetBranchAddress("eventnumber", &eventnumber);
  tree->SetBranchStatus("tracker_side",1);
  tree->SetBranchAddress("tracker_side", &tracker_side);
  tree->SetBranchStatus("tracker_layer",1);
  tree->SetBranchAddress("tracker_layer", &tracker_layer);
  tree->SetBranchStatus("tracker_column",1);
  tree->SetBranchAddress("tracker_column", &tracker_column);
  tree->SetBranchStatus("R0",1);
  tree->SetBranchAddress("R0", &R0);
  tree->SetBranchStatus("R5",1);
  tree->SetBranchAddress("R5", &R5);
  tree->SetBranchStatus("R6",1);
  tree->SetBranchAddress("R6", &R6);
  tree->SetBranchStatus("cell_number",1);
  tree->SetBranchAddress("cell_number", &cell_number);

  TFile *newfile = new TFile("tried_R.root", "RECREATE");
  TTree Result_tree("Result_tree","");
  Result_tree.Branch("R0", &R0_new);
  Result_tree.Branch("R5", &R5_new);
  Result_tree.Branch("R6", &R6_new);
  Result_tree.Branch("cell_number", &cell);

  double R0_tab[2034];
  double R5_tab[2034];
  double R6_tab[2034];


  for (int i = 0; i < tree->GetEntries(); i++) {
    tree->GetEntry(i);
    if (i % 100000 == 0) {
      cout << 100.*i/tree->GetEntries() << "%" << endl;
    }
    for (int j = 0; j < R0->size(); j++) {
      R0_tab[cell_number->at(j)] ++;
      if (R5->at(j) > 0) {
        R5_tab[cell_number->at(j)] ++;
      }
      if (R6->at(j) > 0) {
        R6_tab[cell_number->at(j)] ++;
      }
    }
  }

  for (int i = 0; i < 2034; i++) {
    R0_new = R0_tab[i];
    R5_new = R5_tab[i];
    R6_new = R6_tab[i];
    cell = i;
    Result_tree.Fill();
  }

  newfile->cd();
  Result_tree.Write();
  newfile->Close();

}

int main() {
  track_cutter();
  // tri();
  // R_calculator();
  return 0;
}
