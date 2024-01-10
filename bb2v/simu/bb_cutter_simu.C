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

bool dead_cell_num[2034];
bool cathodless_cell_num[2034];
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
  int colonne_gauche[6] = {4, 22, 42, 60, 78, 97};
  int colonne_droite[6] = {13, 33, 52, 72, 91, 110};

  for (int i = 0; i < 6; i++) {
    if (column > colonne_gauche[i] && column < colonne_droite[i]) {
      return i;
    }
  }
  return 10000;
}

void first_z_selectionner() {

  double mean, sigma;

  TFile *file = new TFile("../../first_z/z_distrib.root", "READ");
  TTree* tree = (TTree*)file->Get("Result_tree");
  tree->SetBranchStatus("*",0);
  tree->SetBranchStatus("mean",1);
  tree->SetBranchAddress("mean", &mean);
  tree->SetBranchStatus("sigma",1);
  tree->SetBranchAddress("sigma", &sigma);

  int compteur = 0;
  for(int i = 0; i < 6; i++){                   //loop on source column
    for(int j = 1; j < 8; j++){                 //loop on source row
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
  if (column > -1 && column < 8) {
    for(int j = 1; j < 8; j++){ // loop on the source row
      // if (z < -0.5 && z > -1)cout << " j = " << j << " and min z = " << source_column_row_mean[column][j] - 3*source_column_row_sigma[column][j] << " and max z = " << source_column_row_mean[column][j] + 3*source_column_row_sigma[column][j] << endl;
      if (z > source_column_row_mean[column][j] - 3*source_column_row_sigma[column][j] && z < source_column_row_mean[column][j] + 3*source_column_row_sigma[column][j]) {
        source_number = column*100 + j;
        return source_number;
      }
      else source_number = -10;
    }
  }
  return source_number;
}

void source_track_selectionner(vector<double> track_side, vector<double> track_layer, vector<double> track_column, int compteur_source[], vector<double> z_ini, double timestamp, double *z_column) {

  vector<long> first_z, second_z;
  vector<double> first_column, second_column;
  vector<int> first_side, second_side;

  for (int i = 0; i < track_layer.size(); i++) {
    int cellnum = track_side.at(i)*1017 + track_layer.at(i) + track_column.at(i)*9;
    if (dead_cell_num[cellnum] == false) continue;
    if (track_layer.at(i) < 4) {
      compteur_source[int (track_side.at(i))]++;
      if (track_layer.at(i) == 0 && cathodless_cell_num[cellnum] == true) {
        first_side.push_back(track_side.at(i));
        first_column.push_back(track_column.at(i));
        first_z.push_back(z_ini.at(i));
      }
      if (track_layer.at(i) == 1 && cathodless_cell_num[cellnum] == true) {
        second_side.push_back(track_side.at(i));
        second_column.push_back(track_column.at(i));
        second_z.push_back(z_ini.at(i));

      }
    }
  }


  z_column[0] = -10000;
  z_column[1] = -10000;
  z_column[2] = -10000;
  z_column[3] = -10000;

  double column;
  double z_gg;
  if (first_column.size() > 0) {
    if (first_column.size() > 1) {      // if two cell on the first layer in front of the OM
      for (int i = 0; i < first_column.size(); i++) {
        z_gg += first_z.at(i);
        column += first_column.at(i);
      }
      z_column[0] = z_gg / first_column.size();
      z_column[1] = column / first_column.size();
    }
    else{
      z_column[0] = first_z.at(0);
      z_column[1] = first_column.at(0);
    }
  }

  if (second_column.size() > 0) {
    if (second_column.size() > 1) {      // if two cell on the second layer in front of the OM
      for (int i = 0; i < second_column.size(); i++) {
        z_gg += second_z.at(i);
        column += second_column.at(i);
      }
      z_column[2] = z_gg / second_column.size();
      z_column[3] = column / second_column.size();
    }
    else{
      z_column[2] = second_z.at(0);
      z_column[3] = second_column.at(0);
    }
  }


}

void calo_track_selectionner(vector<double> track_side, vector<double> track_layer, vector<double> track_column, int compteur_calo[][20], vector<double> z_ini, double timestamp, double *z_column, int *flag6){        /// comparison between the calo column and the tracker layer
  vector<double> last_z;
  vector<double> last_column;
  vector<int> last_side;
  vector<double> z;
  vector<double> col;

  int tab_bas[20] = {0,4,8,14,21,26,31,38,44,49,56,62,67,72,79,85,91,97,102,108};
  int tab_haut[20] = {5,10,16,22,28,33,39,45,51,57,63,67,73,80,86,92,98,104,109,112};


  for (int i = 0; i < track_side.size(); i++) {
      int cellnum = track_side.at(i)*1017 + track_layer.at(i) + track_column.at(i)*9;
      if (dead_cell_num[cellnum] == false) continue;
      if (track_layer.at(i) > 4) {
        for (int j = 0; j < 20; j++) {
          if (track_column.at(i) >= tab_bas[j] && track_column.at(i) <= tab_haut[j]) compteur_calo[int (track_side.at(i))][j]++;
        }
        // cout << "layer : " << track_layer.at(i) << endl;
        if (track_layer.at(i) == 8 && cathodless_cell_num[cellnum] == true) {
          *flag6 = 1;
          last_column.push_back(track_column.at(i));
          last_side.push_back(track_side.at(i));
          last_z.push_back(z_ini.at(i));
        }
      }
    }

  // out << "source : "<< *flag6 << endl;


  if (*flag6 == 9) *flag6 = 0;
  for (int i = 0; i < last_column.size(); i++) {
    for (int j = 0; j < 20; j++) {
      if (last_column.at(i) >= tab_bas[j] && last_column.at(i) <= tab_haut[j] && compteur_calo[last_side.at(i)][j] > 3){
        z.push_back(last_z.at(i));
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

int calo_source_track(int compteur_source[], int compteur_calo[][20], int om_num, int* flag4, int *flag5){        /// source localiser

  int calo_side, calo_column;
  if (om_num < 260) calo_side = 0;
  else calo_side = 1;
  calo_column = (om_num - calo_side*260) / 13;
  if (compteur_calo[calo_side][calo_column] >= 3) *flag4 = 1;
  else if(compteur_calo[calo_side][calo_column] >= 1) *flag4 = 2;
  else *flag4 = 0;
  if (compteur_calo[calo_side][calo_column] >= 3) {
    for (int i = 0; i < 20; i++) {
      if (compteur_source[calo_side] >= 2){
        *flag5 = 1;
        return 1;
      }
      else *flag5 = 0;
    }
  }
  else *flag5 = 0;
  return 0;
}

int gg_counter(double timestamp, vector<vector<long>> R0, vector<int> tracker_side, vector<int> tracker_layer, vector<int> tracker_column, int om_num){
  int gg_counter = 0;
  for (int i = 0; i < R0.size(); i++) {
    int cellnum = tracker_side.at(i)*1017 + tracker_layer.at(i) + tracker_column.at(i)*9;
    if (dead_cell_num[cellnum] == false)  continue;
    if ((2*R0.at(i).at(0) - timestamp)*6.25e-9 > (-0.2e-6) && (2*R0.at(i).at(0) - timestamp)*6.25e-9 < (50e-6) && tracker_side.at(i) == om_num/260) {
      gg_counter++;
    }
  }
  return gg_counter;
}

void dead_catohdless_cell(){

  double R0, R5, R6;
  int cell_number;
  TFile *newfile = new TFile("../../tried_R_1058.root", "READ");
  TTree* tree = (TTree*)newfile->Get("Result_tree");
  tree->SetBranchStatus("*",0);
  tree->SetBranchStatus("R0",1);
  tree->SetBranchAddress("R0", &R0);
  tree->SetBranchStatus("R5",1);
  tree->SetBranchAddress("R5", &R5);
  tree->SetBranchStatus("R6",1);
  tree->SetBranchAddress("R6", &R6);
  tree->SetBranchStatus("cell_number",1);
  tree->SetBranchAddress("cell_number", &cell_number);

  for (int i = 0; i < tree->GetEntries(); i++) {
    tree->GetEntry(i);
    if (R0 < 2000) {
      dead_cell_num[i] = false;
    }
    dead_cell_num[2006] = false;
    dead_cell_num[2009] = false;
    if (R5/R0 < 0.85 || R6/R0 < 0.85 ){
      cathodless_cell_num[i] = false;
    }
    cathodless_cell_num[2006] = false;
  }
}

void track_cutter() {

  first_z_selectionner();

  memset (cathodless_cell_num, true, 2034);
  memset (dead_cell_num, true, 2034);
  dead_catohdless_cell();

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
  std::vector<double> *z_gg = new std::vector<double>;
  std::vector<int> *calo_nohit_om_time = new std::vector<int>;

  TFile *file = new TFile("simu_Tl_208_10G.root", "READ");
  int eventnumber;
  TTree* tree = (TTree*)file->Get("Result_tree");
  tree->SetBranchStatus("*",0);
  tree->SetBranchStatus("calo_nohits",1);
  tree->SetBranchAddress("calo_nohits", &calo_nohits);
  tree->SetBranchStatus("om_id",1);
  tree->SetBranchAddress("om_id", &om_num);
  tree->SetBranchStatus("energyvis_ubc",1);
  tree->SetBranchAddress("energyvis_ubc", &energy_ubc);
  tree->SetBranchStatus("energyvis",1);
  tree->SetBranchAddress("energyvis", &energy);
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
  tree->SetBranchStatus("z_gg",1);
  tree->SetBranchAddress("z_gg", &z_gg);
  std::vector<int> flag;
  std::vector<int> flag0;
  std::vector<int> flag1;
  std::vector<int> flag3;
  std::vector<int> flag4;
  std::vector<int> flag5;
  std::vector<int> flag6;
  std::vector<int> flag7;

  int compteur_source[2] = {0};
  int compteur_calo[2][20] = {0};
  double* z_column_last = new double[2];
  double* z_column_first = new double[4];
  TFile *newfile = new TFile("cut_Tl_simu.root", "RECREATE");
  int calo_compteur_nohits = 0;
  double z_last_gg, z_first_gg, z_second_gg;
  std::vector<int> om_number;
  std::vector<int> charge;
  std::vector<int> amplitude;
  std::vector<int> calo_tdc;
  std::vector<int> associated_track;
  std::vector<int> source;
  std::vector<double> vec_z_last_gg;
  std::vector<double> last_column;
  std::vector<double> vec_z_first_gg;
  std::vector<double> first_column;
  std::vector<double> vec_z_second_gg;
  std::vector<double> second_column;

  TTree Result_tree("Result_tree","");
  Result_tree.Branch("eventnumber", &eventnumber);
  Result_tree.Branch("ncalo_tot", &calo_compteur_nohits);
  Result_tree.Branch("ngg_tot", &tracker_nohits);
  Result_tree.Branch("om_number", &om_num);
  Result_tree.Branch("energyvis", &energy);
  Result_tree.Branch("energyvis_u", &energy_u);
  Result_tree.Branch("energyvis_bc", &energy_bc);
  Result_tree.Branch("energyvis_ubc", &energy_ubc);
  Result_tree.Branch("calo_tdc", &calo_tdc);
  Result_tree.Branch("gg_associated", &associated_track);
  Result_tree.Branch("z_last_gg", &vec_z_last_gg);
  Result_tree.Branch("last_column", &last_column);
  Result_tree.Branch("z_second_gg", &vec_z_second_gg);
  Result_tree.Branch("second_column", &second_column);
  Result_tree.Branch("z_first_gg", &vec_z_first_gg);
  Result_tree.Branch("first_column", &first_column);
  Result_tree.Branch("flag", &flag);
  Result_tree.Branch("flag_e_event", &flag0);
  Result_tree.Branch("flag_charge", &flag1);
  Result_tree.Branch("flag_MW", &flag3);
  Result_tree.Branch("flag_calo_square", &flag4);
  Result_tree.Branch("flag_source_square", &flag5);
  Result_tree.Branch("flag_last_column", &flag6);
  Result_tree.Branch("flag_last_z", &flag7);
  Result_tree.Branch("calo_nohit_om_time", &calo_nohit_om_time);

  for (int i = 0; i < tree->GetEntries(); i++) {      //loop on event number
  // for (int i = 14034; i < 14035; i++) {      //loop on event number

  if (i % 100000 == 0) {
    cout << 100.*i/tree->GetEntries() << "%" << endl;
  }
    tree->GetEntry(i);

    for (int k = 0; k < calo_nohits->at(0); k++) {      //k : loop on calo hit number
      z_column_last[0] = 0;
      z_column_last[1] = 0;
      z_column_first[0] = 0;
      z_column_first[1] = 0;
      z_column_first[2] = 0;
      z_column_first[3] = 0;

      last_column.push_back(-10000);
      first_column.push_back(-10000);
      second_column.push_back(-10000);
      associated_track.push_back(-10000);
      vec_z_last_gg.push_back(-10000);
      vec_z_first_gg.push_back(-10000);
      vec_z_second_gg.push_back(-10000);

      int flag4_int = 0;
      int flag5_int = 0;
      int flag6_int = 0;

      flag0.push_back(9);
      flag1.push_back(9);
      flag3.push_back(9);
      flag4.push_back(9);
      flag5.push_back(9);
      flag6.push_back(9);
      flag7.push_back(9);
      calo_nohit_om_time->push_back(0);
      calo_tdc.push_back(calo_timestamp->at(k));


      for (int k2 = 0; k2 < calo_nohits->at(0); k2++) {
        if (k2 == k) continue;
              if (abs(calo_timestamp->at(k2) - calo_timestamp->at(k)) < 16) calo_nohit_om_time->back()++;
      }

      // if (energy->at(k) > 0.1){
      if (energy->at(k) > 0.1 && om_num->at(k) != 74 && om_num->at(k) != 98 && om_num->at(k) != 80 && om_num->at(k) != 119 && om_num->at(k) != 136 && om_num->at(k) != 146 && om_num->at(k) != 318 && om_num->at(k) != 369 && om_num->at(k) != 461) {      // condition to cut small charge and keep only MW OM
        flag1.back()= 1;
        int timed_gg;
        memset(compteur_source, 0, sizeof(int) * 2);
        memset(compteur_calo, 0, sizeof(int) * 2 * 20);

        source_track_selectionner(*tracker_side, *tracker_layer, *tracker_column, compteur_source, *z_gg, calo_timestamp->at(k), z_column_first);
        calo_track_selectionner(*tracker_side, *tracker_layer, *tracker_column, compteur_calo, *z_gg, calo_timestamp->at(k), z_column_last, &flag6_int);
        calo_source_track(compteur_source, compteur_calo, om_num->at(k), &flag4_int, &flag5_int);
        flag6.back() = flag6_int;
        z_last_gg = z_column_last[0]/1500.;
        z_first_gg = z_column_first[0]/1500.;
        z_second_gg  = z_column_first[3]/1500.;
        calo_compteur_nohits++;

        if (om_num->at(k) < 520 && om_num->at(k) % 13 != 0 && om_num->at(k) % 13 != 12) flag3.back() = 1;
        else if (flag3.at(k) == 9) flag3.back() = 1;

        if (om_num->at(k) < 520 && om_num->at(k) % 13 != 0 && om_num->at(k) % 13 != 12){

          double calo_row2z = -15;
          double z_last_gg_min = -15;
          double z_last_gg_max = -15;

          double tab_z[11] = {0.86, 0.68, 0.52, 0.33, 0.17, 0, 0.18, 0.34, 0.52, 0.67, 0.84};

          if (flag3.back() == 1) {
            calo_row2z = 0.17*(om_num->at(k)%13) - 1.02;
            z_last_gg_min = calo_row2z - 0.091;
            z_last_gg_max = calo_row2z + 0.091;
          }



          last_column.back() = z_column_last[1];
          first_column.back() = z_column_first[2];
          second_column.back() = z_column_first[3];
          associated_track.back() = timed_gg;
          vec_z_last_gg.back() = z_last_gg;
          vec_z_first_gg.back() = z_first_gg;
          vec_z_second_gg.back() = z_second_gg;

          if (flag7.at(k) == 9) flag7.back() = 0;
          if (0.1 >= abs(abs(z_last_gg) - tab_z[om_num->at(k)]) ) {
            flag7.back() = 1;
          }
          if (calo_source_track(compteur_source, compteur_calo, om_num->at(k), &flag4_int, &flag5_int) == 1 && z_last_gg <= z_last_gg_max && z_last_gg >= z_last_gg_min && flag6_int == 1){//  && source_number >= 0) {       //// Delta z moyen = 0.268589cm 0.29 à 3 sigma, soit 0.1852 (1)de -1 à 1 de moyenneur gategauss
            // cout << "ok" << endl;
            flag0.back() = 1;
            // cout << flag0.back() << endl;;
          }
          else {
            flag0.back() = 0;
          }
          flag4.back() = flag4_int;
          flag5.back() = flag5_int;
        }
        else {
          flag0.back() = 0;
        }

        // source.push_back
      }
      else {
        flag1.back() = 0;
      }
      flag.push_back(flag0.at(k) +flag1.at(k)*10 + flag3.at(k)*1000 + flag4.at(k)*10000 + flag5.at(k)*100000 + flag6.at(k)*1000000 + flag7.at(k)*10000000);
    }

    // for (size_t j = 0; j < flag6.size(); j++) {
    //   cout << "event : " << j << " flag 6 " << flag6.at(j) << endl;
    // }


    Result_tree.Fill();
    flag0.clear();
    flag1.clear();
    flag3.clear();
    flag4.clear();
    flag5.clear();
    flag6.clear();
    flag7.clear();
    flag.clear();
    om_number.clear();
    charge.clear();
    amplitude.clear();
    calo_tdc.clear();
    associated_track.clear();
    last_column.clear();
    vec_z_last_gg.clear();
    first_column.clear();
    vec_z_first_gg.clear();
    second_column.clear();
    vec_z_second_gg.clear();
    calo_compteur_nohits = 0;
    calo_nohit_om_time->clear();
    energy->clear();
  }


  newfile->cd();
  Result_tree.Write();
  newfile->Close();
  cout << "OK" << endl;
}

void tri(string run) {
  gStyle->SetOptFit(1);
  gStyle->SetOptStat(1);


  TFile *file = new TFile(Form("data/snemo_run-%s_udd.root",run.c_str()), "READ");
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


  TFile *newfile = new TFile(Form("data/new_timestamp_%s.root",run.c_str()), "RECREATE");
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
              R0_map->Fill(tracker_column->at(k), tracker_side->at(k)*9 + tracker_layer->at(k));
              if (R5->at(k).at(0) > 0) {
                R0R5_map->Fill(tracker_column->at(k), tracker_side->at(k)*9 + tracker_layer->at(k));
              }
              if (R6->at(k).at(0) > 0) {
                R0R6_map->Fill(tracker_column->at(k), tracker_side->at(k)*9 + tracker_layer->at(k));
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
void R_calculator(string run) {

TFile *file = new TFile(Form("data/new_timestamp_%s.root",run.c_str()), "READ");

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

  TFile *newfile = new TFile(Form("tried_R_%s.root",run.c_str()), "RECREATE");
  TTree Result_tree("Result_tree","");
  Result_tree.Branch("R0", &R0_new);
  Result_tree.Branch("R5", &R5_new);
  Result_tree.Branch("R6", &R6_new);
  Result_tree.Branch("cell_number", &cell);

  double R0_tab[2034];
  double R5_tab[2034];
  double R6_tab[2034];
  memset (R0_tab, 0, 2034*sizeof(double));
  memset (R5_tab, 0, 2034*sizeof(double));
  memset (R6_tab, 0, 2034*sizeof(double));


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


int main(int argc, char const *argv[]) {

  string run = "";
  int full = 0;
  for(int i = 0; i<argc; i++){
    if (std::string(argv[i]) == "-r" ){
      run = std::string(argv[i+1]);
    }
    if (std::string(argv[i]) == "-f" ){
      full = 1;
    }
  }

  track_cutter();

  return 0;
}
