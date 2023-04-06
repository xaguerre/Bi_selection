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
bool dead_cell_num[2034];
bool cathodless_cell_num[2034];

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

  TFile *file = new TFile("../first_z/z_distrib.root", "READ");
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

void source_track_selectionner(vector<double> track_side, vector<double> track_layer, vector<double> track_column, int compteur_source[][6], vector<double> x_ini, vector<double> y_ini, vector<double> z_ini, double timestamp, double *column_position) {
  vector<double> first_x;
  vector<double> first_y;
  vector<double> first_z;
  vector<double> first_column;
  vector<int> first_side;
  vector<double> col;
  vector<double> x;
  vector<double> y;
  vector<double> z;


  for (int i = 0; i < track_layer.size(); i++) {
    int cellnum = track_side.at(i)*1017 + track_layer.at(i) + track_column.at(i)*9;
    if (dead_cell_num[cellnum] == false) continue;
    if (track_layer.at(i) < 4) {
      for (int j = 0; j < 6; j++) {
        if (abs(track_column.at(i) - (8.5 + j * 19)) < 3) {
          compteur_source[int (track_side.at(i))][j]++;
          if (track_layer.at(i) == 0 && cathodless_cell_num[cellnum] == true) {
            first_side.push_back(track_side.at(i));
            first_column.push_back(track_column.at(i));
            first_x.push_back(x_ini.at(i));
            first_y.push_back(y_ini.at(i));
            first_z.push_back(z_ini.at(i));
          }
        }
      }
    }
  }

  for (int i = 0; i < first_column.size(); i++) {
    for (int j = 0; j < 6; j++) {
      if (abs(first_column.at(i)- (8.5 + j * 19)) < 3  && compteur_source[int(first_side.at(i))][j] > 3){
        z.push_back(first_z.at(i));
        y.push_back(first_y.at(i));
        x.push_back(first_x.at(i));
        col.push_back(first_column.at(i));
      }
    }
  }

  column_position[0] = -1000;
  column_position[1] = -1000;
  column_position[2] = -1000;
  column_position[3] = -1000;

  double column;
  double z_gg;
  double y_gg;
  double x_gg;

  if (z.size() > 0) {
    if (z.size() > 1) {      // if two cell on the last layer in front of the OM
      for (int i = 0; i < z.size(); i++) {
        z_gg += z.at(i);
        y_gg += y.at(i);
        x_gg += x.at(i);
        column += col.at(i);
      }
      column_position[0] = x_gg / z.size();
      column_position[1] = y_gg / z.size();
      column_position[2] = z_gg / z.size();
      column_position[3] = column / z.size();
    }
    else{
      column_position[0] = x.at(0);
      column_position[1] = y.at(0);
      column_position[2] = z.at(0);
      column_position[3] = col.at(0);
    }
  }
}

void calo_track_selectionner(vector<double> track_side, vector<double> track_layer, vector<double> track_column, int compteur_calo[][20], vector<double> x_ini, vector<double> y_ini, vector<double> z_ini, double timestamp, double *column_position, int *flag6){        /// comparison between the calo column and the tracker layer

  vector<double> last_x;
  vector<double> last_y;
  vector<double> last_z;
  vector<double> last_column;
  vector<int> last_side;
  vector<double> col;
  vector<double> x;
  vector<double> y;
  vector<double> z;

  int tab_bas[20] = {0,4,8,14,21,26,31,38,44,49,56,62,67,72,79,85,91,97,102,108};
  int tab_haut[20] = {5,10,16,22,28,33,39,45,51,57,63,67,73,80,86,92,98,104,109,112};


  for (int i = 0; i < track_layer.size(); i++) {
    int cellnum = track_side.at(i)*1017 + track_layer.at(i) + track_column.at(i)*9;
    if (dead_cell_num[cellnum] == false) continue;
    if (track_layer.at(i) > 4) {
      for (int j = 0; j < 20; j++) {
        if (track_column.at(i) >= tab_bas[j] && track_column.at(i) <= tab_haut[j]) compteur_calo[int (track_side.at(i))][j]++;
      }
      // if (track_layer.at(i) == 8) cout << "cell num = " << cellnum << " and cathode = " << cathodless_cell_num[cellnum]
      if (track_layer.at(i) == 8 && cathodless_cell_num[cellnum] == true) {
        *flag6 = 1;
        last_side.push_back(track_side.at(i));
        last_column.push_back(track_column.at(i));
        last_x.push_back(x_ini.at(i));
        last_y.push_back(y_ini.at(i));
        last_z.push_back(z_ini.at(i));
      }
    }
  }


  if (*flag6 == 9) *flag6 = 0;
  for (int i = 0; i < last_column.size(); i++) {
    for (int j = 0; j < 20; j++) {
      if (last_column.at(i) >= tab_bas[j] && last_column.at(i) <= tab_haut[j] && compteur_calo[last_side.at(i)][j] > 3){
        z.push_back(last_z.at(i));
        y.push_back(last_y.at(i));
        x.push_back(last_x.at(i));
        col.push_back(last_column.at(i));
      }
    }
  }

  column_position[0] = -1000;
  column_position[1] = -1000;
  column_position[2] = -1000;
  column_position[3] = -1000;

  double column;
  double z_gg;
  double y_gg;
  double x_gg;

  if (z.size() > 0) {
    if (z.size() > 1) {      // if two cell on the las layer in front of the OM
      for (int i = 0; i < z.size(); i++) {
        z_gg += z.at(i);
        y_gg += y.at(i);
        x_gg += x.at(i);
        column += col.at(i);
      }
      column_position[0] = x_gg / z.size();
      column_position[1] = y_gg / z.size();
      column_position[2] = z_gg / z.size();
      column_position[3] = column / z.size();
    }
    else{
      column_position[0] = x.at(0);
      column_position[1] = y.at(0);
      column_position[2] = z.at(0);
      column_position[3] = col.at(0);
    }
  }
}

int calo_source_track(int compteur_source[][6], int compteur_calo[][20], int om_num, int* flag4, int *flag5){        /// source localiser

  int calo_side, calo_column;
  if (om_num < 260) calo_side = 0;
  else calo_side = 1;
  calo_column = (om_num - calo_side*260) / 13;
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

int gg_counter(double timestamp, vector<double> tracker_side, vector<double> tracker_layer, vector<double> tracker_column, int om_num){
  int gg_counter = 0;
  for (int i = 0; i < tracker_side.size(); i++) {
    int cellnum = tracker_side.at(i)*1017 + tracker_layer.at(i) + tracker_column.at(i)*9;
    if (dead_cell_num[cellnum] == false)  continue;
    if (tracker_side.at(i) == om_num/260) {
      gg_counter++;
    }
  }
  return gg_counter;
}

void dead_catohdless_cell(){

  double R0, R5, R6;
  int cell_number;

  TFile *newfile = new TFile("tried_R_974_fake.root", "READ");
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
    if (R0 < 10000) {
      dead_cell_num[i] = false;
    }
    if (R5/R0 < 0.85 || R6/R0 < 0.85 ){
      cathodless_cell_num[i] = false;
    }
  }
}


void track_cutter() {
  first_z_selectionner();
  // return;
  memset (dead_cell_num, true, 2034);
  memset (cathodless_cell_num, true, 2034);
  dead_catohdless_cell();

  TFile *file = new TFile("data/new_simu.root", "READ");

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

  int eventnumber;
  TTree* tree = (TTree*)file->Get("Result_tree");
  tree->SetBranchStatus("*",0);
  tree->SetBranchStatus("calo_nohits",1);
  tree->SetBranchAddress("calo_nohits", &calo_nohits);
  tree->SetBranchStatus("om_id",1);
  tree->SetBranchAddress("om_id", &om_num);
  tree->SetBranchStatus("energyvis_ubc",1);
  tree->SetBranchAddress("energyvis_ubc", &energy_ubc);
  tree->SetBranchStatus("energyvis_u",1);
  tree->SetBranchAddress("energyvis_u", &energy_u);
  tree->SetBranchStatus("energyvis_bc",1);
  tree->SetBranchAddress("energyvis_bc", &energy_bc);
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
  tree->SetBranchStatus("x_gg",1);
  tree->SetBranchAddress("x_gg", &x_gg);
  tree->SetBranchStatus("y_gg",1);
  tree->SetBranchAddress("y_gg", &y_gg);
  tree->SetBranchStatus("z_gg",1);
  tree->SetBranchAddress("z_gg", &z_gg);

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
  double* column_position_last = new double[4];
  double* column_position_first = new double[4];
  int  nelec, ngamma;
  TFile *newfile = new TFile("cut_bi_cells_fake_323.root", "RECREATE");
  int calo_compteur_nohits = 0;
  double x_last_gg, x_first_gg, y_last_gg, y_first_gg, z_last_gg, z_first_gg;
  std::vector<int> charge;
  std::vector<int> amplitude;
  std::vector<int> calo_tdc;
  std::vector<int> associated_track;
  std::vector<int> source;
  std::vector<int> vec_nelec;
  std::vector<int> vec_ngamma;
  std::vector<int> e_event;
  std::vector<double> vec_x_last_gg;
  std::vector<double> vec_y_last_gg;
  std::vector<double> vec_z_last_gg;
  std::vector<double> last_column;
  std::vector<double> vec_x_first_gg;
  std::vector<double> vec_y_first_gg;
  std::vector<double> vec_z_first_gg;
  std::vector<double> first_column;
  std::vector<double> cell_323;
  std::vector<double> cell_1963;
  std::vector<double> cell_412;
  std::vector<double> cell_1070;

  TTree Result_tree("Result_tree","");
  Result_tree.Branch("eventnumber", &eventnumber);
  Result_tree.Branch("ncalo_tot", &calo_compteur_nohits);
  Result_tree.Branch("ngg_tot", &tracker_nohits);
  Result_tree.Branch("nelec", &nelec);
  Result_tree.Branch("ngamma", &ngamma);
  Result_tree.Branch("e_event", &e_event);
  Result_tree.Branch("om_num", &om_num);
  Result_tree.Branch("energyvis", &energy);
  Result_tree.Branch("energyvis_u", &energy_u);
  Result_tree.Branch("energyvis_bc", &energy_bc);
  Result_tree.Branch("energyvis_ubc", &energy_ubc);
  Result_tree.Branch("calo_tdc", &calo_tdc);
  Result_tree.Branch("gg_associated", &associated_track);
  Result_tree.Branch("x_last_gg", &vec_x_last_gg);
  Result_tree.Branch("y_last_gg", &vec_y_last_gg);
  Result_tree.Branch("z_last_gg", &vec_z_last_gg);
  Result_tree.Branch("last_column", &last_column);
  Result_tree.Branch("x_first_gg", &vec_x_first_gg);
  Result_tree.Branch("y_first_gg", &vec_y_first_gg);
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
  Result_tree.Branch("cell_323", &cell_323);
  Result_tree.Branch("cell_1963", &cell_1963);
  Result_tree.Branch("cell_412", &cell_412);
  Result_tree.Branch("cell_1070", &cell_1070);

  for (int i = 0; i < tree->GetEntries()/10; i++) {      //loop on event number
    // for (int i = 0; i < 10000000; i++) {      //loop on event number

    if (i % 100000 == 0) {
      cout << 100.*i/tree->GetEntries() << "%" << endl;
    }
    tree->GetEntry(i);

    for (int k = 0; k < calo_nohits->at(0); k++) {      //k : loop on calo hit number
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

      calo_tdc.push_back(calo_timestamp->at(k));
      int cellcompteur_323 = 0;
      int cellcompteur_1963 = 0;
      int cellcompteur_412 = 0;
      int cellcompteur_1070 = 0;

      for (int l = 0; l < tracker_side->size(); l++) {
        int cellnum = tracker_side->at(l)*1017 + tracker_layer->at(l) + tracker_column->at(l)*9;
        if ( cellnum== 323 && dead_cell_num[cellnum] == true) cell_323.push_back(1);
        if (tracker_side->at(l)*1017 + tracker_layer->at(l) + tracker_column->at(l)*9 == 1963) cell_1963.push_back(1);
        if (tracker_side->at(l)*1017 + tracker_layer->at(l) + tracker_column->at(l)*9 == 412) cell_412.push_back(1);
        if (tracker_side->at(l)*1017 + tracker_layer->at(l) + tracker_column->at(l)*9 == 1070) cell_1070.push_back(1);
        cellcompteur_323++;
        cellcompteur_1963++;
        cellcompteur_412++;
        cellcompteur_1070++;
      }
      if (cellcompteur_323 == 0) cell_323.push_back(0);
      if (cellcompteur_1963 == 0) cell_1963.push_back(0);
      if (cellcompteur_412 == 0) cell_412.push_back(0);
      if (cellcompteur_1070 == 0) cell_1070.push_back(0);

      // if (energy->at(k) > 0.1){
      if (energy->at(k) > 0.1 && om_num->at(k) != 74 && om_num->at(k) != 98 && om_num->at(k) != 80 && om_num->at(k) != 119 && om_num->at(k) != 136 && om_num->at(k) != 146 && om_num->at(k) != 318 && om_num->at(k) != 369 && om_num->at(k) != 461) {      // condition to cut small charge and keep only MW OM
        flag1.pop_back();
        flag1.push_back(1);
        int timed_gg;
        memset(compteur_source, 0, sizeof(int) * 2 * 6);
        memset(compteur_calo, 0, sizeof(int) * 2 * 20);
        memset(column_position_last, 0, sizeof(double) * 4);
        memset(column_position_first, 0, sizeof(double) * 4);

        source_track_selectionner(*tracker_side, *tracker_layer, *tracker_column, compteur_source, *x_gg, *y_gg, *z_gg, calo_timestamp->at(k), column_position_first);
        calo_track_selectionner(*tracker_side, *tracker_layer, *tracker_column, compteur_calo, *x_gg, *y_gg, *z_gg, calo_timestamp->at(k), column_position_last, &flag6_int);

        flag6.pop_back();
        flag6.push_back(flag6_int);
        x_last_gg = column_position_last[0];
        x_first_gg = column_position_first[0];
        y_last_gg = column_position_last[1];
        y_first_gg = column_position_first[1];
        z_last_gg = column_position_last[2]/1500.;
        z_first_gg = column_position_first[2]/1500.;
        timed_gg = gg_counter(calo_timestamp->at(k), *tracker_side, *tracker_layer, *tracker_column, om_num->at(k));
        calo_compteur_nohits++;
        int first_gg_row = row_to_source(column_position_first[3]);
        int source_number = source_numberer(column_position_first[2], first_gg_row);

        // cout << "timed_gg = " << timed_gg << endl;

        if (timed_gg > 5 && timed_gg < 16) {
          flag2.pop_back();
          flag2.push_back(1);
        }
        else if (flag2.at(k) == 9) {
          flag2.pop_back();
          flag2.push_back(0);
        }
        if (om_num->at(k) < 520 && om_num->at(k) % 13 != 0 && om_num->at(k) % 13 != 12) {
          flag3.pop_back();
          flag3.push_back(1);
        }
        else if (flag3.at(k) == 9) {
          flag3.pop_back();
          flag3.push_back(0);
        }

        double calo_row2z = -15;
        double z_last_gg_min = -15;
        double z_last_gg_max = -15;
        if (om_num->at(k) < 520){
          calo_row2z = -1.1165 +0.18714*((om_num->at(k)%13));
          z_last_gg_min = calo_row2z - 0.1;
          z_last_gg_max = calo_row2z + 0.1;
        }
        if (timed_gg > 5 && timed_gg < 16 && om_num->at(k) < 520 && om_num->at(k) % 13 != 0 && om_num->at(k) % 13 != 12){
          if (z_last_gg <= z_last_gg_max && z_last_gg >= z_last_gg_min) {
            flag7.pop_back();
            flag7.push_back(1);
          }
          else if (flag7.at(k) == 9) {
            flag7.pop_back();
            flag7.push_back(0);
          }
          if (calo_source_track(compteur_source, compteur_calo, om_num->at(k), &flag4_int, &flag5_int) == 1 && z_last_gg <= z_last_gg_max && z_last_gg >= z_last_gg_min){//  && source_number >= 0) {       //// Delta z moyen = 0.268589cm 0.29 à 3 sigma, soit 0.1852 (1)de -1 à 1 de moyenneur gategauss

            vec_x_last_gg.push_back(x_last_gg);
            vec_x_first_gg.push_back(x_first_gg);
            vec_y_last_gg.push_back(y_last_gg);
            vec_y_first_gg.push_back(y_first_gg);
            vec_z_last_gg.push_back(z_last_gg);
            vec_z_first_gg.push_back(z_first_gg);
            // cout << "first col = " << column_position_first[1] << endl;
            vec_nelec.push_back(1);
            last_column.push_back(column_position_last[3]);
            first_column.push_back(column_position_first[3]);
            e_event.push_back(1);
            associated_track.push_back(timed_gg);
            source.push_back(source_number);

            flag0.pop_back();
            flag0.push_back(1);

          }
          else {
            vec_ngamma.push_back(1);
            e_event.push_back(0);
            associated_track.push_back(0);
            source.push_back(-9999);
            flag0.pop_back();
            flag0.push_back(0);
          }
          flag4.pop_back();
          flag4.push_back(flag4_int);
          flag5.pop_back();
          flag5.push_back(flag5_int);
          // cout << "flag 6 = " << flag6.at(k) << endl;
        }
        else {
          vec_ngamma.push_back(1);
          e_event.push_back(0);
          associated_track.push_back(0);
          source.push_back(-9999);
          flag0.pop_back();
          flag0.push_back(0);
        }
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
    vec_nelec.clear();
    vec_ngamma.clear();
    om_num->clear();
    charge.clear();
    amplitude.clear();
    calo_tdc.clear();
    associated_track.clear();
    source.clear();
    e_event.clear();
    last_column.clear();
    vec_x_last_gg.clear();
    vec_y_last_gg.clear();
    vec_z_last_gg.clear();
    first_column.clear();
    vec_x_first_gg.clear();
    vec_y_first_gg.clear();
    vec_z_first_gg.clear();
    cell_323.clear();
    cell_1963.clear();
    cell_412.clear();
    cell_1070.clear();

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
