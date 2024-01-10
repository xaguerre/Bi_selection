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


double energy_convertor[520];


void energy_convertor_filler() {

  int om;
  double mean, Chi2;

  TFile *file4 = new TFile("../Simu/Bi_fit/fitted_bi_gas.root", "READ");
  TTree* tree4 = (TTree*)file4->Get("Result_tree");
  tree4->SetBranchStatus("*",0);
  tree4->SetBranchStatus("om_number",1);
  tree4->SetBranchAddress("om_number", &om);
  tree4->SetBranchStatus("mean",1);
  tree4->SetBranchAddress("mean", &mean);

  double mean_tab_gas[520];
  memset (mean_tab_gas, 0, 520*sizeof(double));

  for (int i = 0; i < tree4->GetEntries(); i++) {
    tree4->GetEntry(i);
    if (mean > 0) {
      mean_tab_gas[om] = mean;
    }
  }

  TFile tree_file("../Bi_fit/fitted_bi_1059.root", "READ");
  int om_number, position;
  TTree* tree = (TTree*)tree_file.Get("Result_tree");
  tree->SetBranchStatus("*",0);
  tree->SetBranchStatus("mean",1);
  tree->SetBranchAddress("mean", &mean);
  tree->SetBranchStatus("om_number",1);
  tree->SetBranchAddress("om_number", &om_number);
  tree->SetBranchStatus("position",1);
  tree->SetBranchAddress("position", &position);
  tree->SetBranchStatus("Chi2_ndf",1);
  tree->SetBranchAddress("Chi2_ndf", &Chi2);

  for (int i = 0; i < tree->GetEntries(); i++) {
    tree->GetEntry(i);
    if (Chi2 > 0 && Chi2 < 5) {
      energy_convertor[om_number] = mean_tab_gas[om_number]/mean;
    }
  }
  return;
}

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

  TFile *file = new TFile("../first_z/z_distrib.root", "READ");
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

void source_track_selectionner(vector<int> track_side, vector<int> track_layer, vector<int> track_column, int compteur_source[], vector<vector<long>> R0, vector<vector<long>> R5, vector<vector<long>> R6, double timestamp, double *z_column) {
  vector<long> first_R5, second_R5;
  vector<long> first_R6, second_R6;
  vector<long> first_R0, second_R0;
  vector<double> first_column, second_column;
  vector<int> first_side, second_side;

  for (int i = 0; i < track_layer.size(); i++) {
    if ((2*R0.at(i).at(0) - timestamp)*6.25e-9 > (-0.2e-6) && (2*R0.at(i).at(0) - timestamp)*6.25e-9 < (50e-6)) {
      int cellnum = track_side.at(i)*1017 + track_layer.at(i) + track_column.at(i)*9;
      if (dead_cell_num[cellnum] == false) continue;
      if (track_layer.at(i) < 4) {
        compteur_source[track_side.at(i)]++;
        if (track_layer.at(i) == 0 && cathodless_cell_num[cellnum] == true) {
          first_side.push_back(track_side.at(i));
          first_column.push_back(track_column.at(i));
          first_R5.push_back(R5.at(i).at(0));
          first_R6.push_back(R6.at(i).at(0));
          first_R0.push_back(R0.at(i).at(0));
        }
        if (track_layer.at(i) == 1 && cathodless_cell_num[cellnum] == true) {
          second_side.push_back(track_side.at(i));
          second_column.push_back(track_column.at(i));
          second_R5.push_back(R5.at(i).at(0));
          second_R6.push_back(R6.at(i).at(0));
          second_R0.push_back(R0.at(i).at(0));
        }
      }
    }
  }

  z_column[0] = -1000;
  z_column[1] = -1000;
  z_column[2] = -1000;
  z_column[3] = -1000;

  double column;
  double z_gg;
  if (first_column.size() > 0) {
    if (first_column.size() > 1) {      // if two cell on the first layer in front of the OM
      for (int i = 0; i < first_column.size(); i++) {
        z_gg += z_calculator_gg(first_R0.at(i), first_R5.at(i), first_R6.at(i));
        column += first_column.at(i);
      }
      z_column[0] = z_gg / first_column.size();
      z_column[1] = column / first_column.size();
    }
    else{
      z_column[0] = z_calculator_gg(first_R0.at(0), first_R5.at(0), first_R6.at(0));
      z_column[1] = first_column.at(0);
    }
  }

  if (second_column.size() > 0) {
    if (second_column.size() > 1) {      // if two cell on the second layer in front of the OM
      for (int i = 0; i < second_column.size(); i++) {
        z_gg += z_calculator_gg(second_R0.at(i), second_R5.at(i), second_R6.at(i));
        column += second_column.at(i);
      }
      z_column[2] = z_gg / second_column.size();
      z_column[3] = column / second_column.size();
    }
    else{
      z_column[2] = z_calculator_gg(second_R0.at(0), second_R5.at(0), second_R6.at(0));
      z_column[3] = second_column.at(0);
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
    if ((2*R0.at(i).at(0) - timestamp)*6.25e-9 > (-0.2e-6) && (2*R0.at(i).at(0) - timestamp)*6.25e-9 < (50e-6)) {
      int cellnum = track_side.at(i)*1017 + track_layer.at(i) + track_column.at(i)*9;
      if (dead_cell_num[cellnum] == false) continue;
      if (track_layer.at(i) > 4) {
        for (int j = 0; j < 20; j++) {
          if (track_column.at(i) >= tab_bas[j] && track_column.at(i) <= tab_haut[j]) compteur_calo[track_side.at(i)][j]++;
        }
        // cout << "layer : " << track_layer.at(i) << endl;
        if (track_layer.at(i) == 8 && cathodless_cell_num[cellnum] == true) {
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
  // cout << "source : "<< *flag6 << endl;


  if (*flag6 == 9) *flag6 = 0;
  for (int i = 0; i < last_column.size(); i++) {
    for (int j = 0; j < 20; j++) {
      if (last_column.at(i) >= tab_bas[j] && last_column.at(i) <= tab_haut[j] && compteur_calo[last_side.at(i)][j] > 3){
        z.push_back(z_calculator_gg(last_R0.at(i), last_R5.at(i), last_R6.at(i)));
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

int calo_source_track(int compteur_source[], int compteur_calo[][20], int calo_column, int calo_side, int* flag4, int *flag5){        /// source localiser
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

void dead_catohdless_cell(string run){

  double R0, R5, R6;
  int cell_number;
  cout << run << endl;
  TFile *newfile = new TFile(Form("../tried_R_%s.root", run.c_str()), "READ");
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

void track_cutter(string run) {
  energy_convertor_filler();

  first_z_selectionner();

  memset (cathodless_cell_num, true, 2034);
  memset (dead_cell_num, true, 2034);
  dead_catohdless_cell(run);
  TFile *file = new TFile(Form("../data/snemo_run-%s_udd.root", run.c_str()), "READ");
  std::vector<vector<short>> *waveform = new std::vector<vector<short>>;
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
  std::vector<int> *calo_nohit_om_time = new std::vector<int>;
  std::vector<int> *rising_cell = new std::vector<int>;
  std::vector<int> *falling_cell = new std::vector<int>;

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
  tree->SetBranchAddress("digicalo.waveform", &waveform);
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
  tree->SetBranchStatus("digicalo.rising_cell",1);
  tree->SetBranchAddress("digicalo.rising_cell", &rising_cell);
  tree->SetBranchStatus("digicalo.falling_cell",1);
  tree->SetBranchAddress("digicalo.falling_cell", &falling_cell);
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
  int om_num;

  TFile *newfile = new TFile(Form("cut_%s_calibrated.root", run.c_str()), "RECREATE");
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
  std::vector<double> energy;

  TTree Result_tree("Result_tree","");
  Result_tree.Branch("eventnumber", &eventnumber);
  Result_tree.Branch("ncalo_tot", &calo_compteur_nohits);
  Result_tree.Branch("ngg_tot", &tracker_nohits);
  Result_tree.Branch("om_number", &om_number);
  Result_tree.Branch("charge", &charge);
  Result_tree.Branch("amplitude", &amplitude);
  Result_tree.Branch("calo_tdc", &calo_tdc);
  Result_tree.Branch("gg_associated", &associated_track);
  Result_tree.Branch("z_last_gg", &vec_z_last_gg);
  Result_tree.Branch("last_column", &last_column);
  Result_tree.Branch("z_first_gg", &vec_z_first_gg);
  Result_tree.Branch("first_column", &first_column);
  Result_tree.Branch("z_second_gg", &vec_z_second_gg);
  Result_tree.Branch("second_column", &second_column);
  Result_tree.Branch("flag", &flag);
  Result_tree.Branch("flag_e_event", &flag0);
  Result_tree.Branch("flag_charge", &flag1);
  Result_tree.Branch("flag_MW", &flag3);
  Result_tree.Branch("flag_calo_square", &flag4);
  Result_tree.Branch("flag_source_square", &flag5);
  Result_tree.Branch("flag_last_column", &flag6);
  Result_tree.Branch("flag_last_z", &flag7);
  Result_tree.Branch("calo_nohit_om_time", &calo_nohit_om_time);
  Result_tree.Branch("energy", &energy);
  Result_tree.Branch("calo_timestamp", &timestamp);
  Result_tree.Branch("rising_cell", &rising_cell);
  Result_tree.Branch("falling_cell", &falling_cell);
  Result_tree.Branch("waveform", &waveform);

  for (int i = 0; i < tree->GetEntries(); i++) {      //loop on event number
  // for (int i = 14034; i < 14035; i++) {      //loop on event number

  if (i % 100000 == 0) {
    cout << 100.*i/tree->GetEntries() << "%" << endl;
  }
    tree->GetEntry(i);

    for (int k = 0; k < calo_nohits; k++) {      //k : loop on calo hit number
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


      for (int k2 = 0; k2 < calo_nohits; k2++) {
        if (k2 == k) continue;
        if (abs(timestamp->at(k2) - timestamp->at(k)) < 16 && -calo_ampl->at(k2) > 200) calo_nohit_om_time->back()++;
      }

      if (calo_type->at(k) == 0) om_num = calo_side->at(k)*260 + calo_column->at(k)*13 + calo_row->at(k);
      if (calo_type->at(k) == 1) om_num = 520 + calo_side->at(k)*64 +  calo_wall->at(k)*32 + calo_column->at(k)*16 + calo_row->at(k);
      if (calo_type->at(k) == 2) om_num = 520 + 128 + calo_side->at(k)*32 + calo_wall->at(k)*16 + calo_column->at(k);
      om_number.push_back(om_num);
      energy.push_back(-calo_charge->at(k)*energy_convertor[om_num]);
      charge.push_back(-calo_charge->at(k));
      amplitude.push_back(-calo_ampl->at(k));
      calo_tdc.push_back(timestamp->at(k));


      if (-calo_ampl->at(k) > 200 && om_num != 74 && om_num != 98 && om_num != 80 && om_num != 119 && om_num != 136 && om_num != 146 && om_num != 318 && om_num != 369 && om_num != 461) {      // condition to cut small charge and keep only MW OM
      // if (-calo_ampl->at(k) > 200) {      // condition to cut small charge
        flag1.back ()= 1;
        int timed_gg;
        memset(compteur_source, 0, sizeof(int) * 2);
        memset(compteur_calo, 0, sizeof(int) * 2 * 20);

        source_track_selectionner(*tracker_side, *tracker_layer, *tracker_column, compteur_source, *R0, *R5, *R6, timestamp->at(k), z_column_first);
        calo_track_selectionner(*tracker_side, *tracker_layer, *tracker_column, compteur_calo, *R0, *R5, *R6, timestamp->at(k), z_column_last, &flag6_int);
        calo_source_track(compteur_source, compteur_calo, calo_column->at(k), calo_side->at(k), &flag4_int, &flag5_int);
        flag6.back() = flag6_int;
        z_last_gg = z_column_last[0];
        z_first_gg = z_column_first[0];
        z_second_gg  = z_column_first[3];
        timed_gg = gg_counter(timestamp->at(k), *R0, *tracker_side, *tracker_layer, *tracker_column, om_num);
        calo_compteur_nohits++;
        int first_gg_row = row_to_source(z_column_first[1]);

        if (om_num < 520 && om_num % 13 != 0 && om_num % 13 != 12) flag3.back() = 1;
        else if (flag3.at(k) == 9) flag3.back() = 1;

        if (om_num < 520 && om_num % 13 != 0 && om_num % 13 != 12){

          double calo_row2z = -15;
          double z_last_gg_min = -15;
          double z_last_gg_max = -15;

          if (calo_type->at(k) == 0){
            calo_row2z = -1.1165 +0.18714*calo_row->at(k);
            z_last_gg_min = calo_row2z - 0.1;
            z_last_gg_max = calo_row2z + 0.1;
          }

          last_column.back() = z_column_last[1];
          first_column.back() = z_column_first[1];
          second_column.back() = z_column_first[3];
          associated_track.back() = timed_gg;
          vec_z_last_gg.back() = z_last_gg;
          vec_z_first_gg.back() = z_first_gg;
          vec_z_second_gg.back() = z_second_gg;

          if (flag7.at(k) == 9) flag7.back() = 0;
          if (z_last_gg <= z_last_gg_max && z_last_gg >= z_last_gg_min) {
            flag7.back() = 1;
          }
          if (calo_source_track(compteur_source, compteur_calo, calo_column->at(k), calo_side->at(k), &flag4_int, &flag5_int) == 1 && z_last_gg <= z_last_gg_max && z_last_gg >= z_last_gg_min && flag6_int == 1){//  && source_number >= 0) {       //// Delta z moyen = 0.268589cm 0.29 à 3 sigma, soit 0.1852 (1)de -1 à 1 de moyenneur gategauss
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
          source.push_back(-10);
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
    energy.clear();
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

//
void waveformer() {

  int ncalo_tot;
  std::vector<int> *flag_e_event = new std::vector<int>;
  std::vector<int> *om_number = new std::vector<int>;
  std::vector<int> *flag_charge = new std::vector<int>;
  std::vector<int> *flag_associated_nohit = new std::vector<int>;
  std::vector<vector<short>> *waveform = new std::vector<vector<short>>;
  std::vector<double> *energy = new std::vector<double>;

  TFile *file = new TFile("cut_974_no.root", "READ");
  TTree* tree = (TTree*)file->Get("Result_tree");
  tree->SetBranchStatus("*",0);
  tree->SetBranchStatus("ncalo_tot",1);
  tree->SetBranchAddress("ncalo_tot", &ncalo_tot);
  tree->SetBranchStatus("om_number",1);
  tree->SetBranchAddress("om_number", &om_number);
  tree->SetBranchStatus("flag_e_event",1);
  tree->SetBranchAddress("flag_e_event", &flag_e_event);
  tree->SetBranchStatus("flag_associated_nohit",1);
  tree->SetBranchAddress("flag_associated_nohit", &flag_associated_nohit);
  tree->SetBranchStatus("waveform",1);
  tree->SetBranchAddress("waveform", &waveform);
  tree->SetBranchStatus("energy",1);
  tree->SetBranchAddress("energy", &energy);
  TH2D* wave = new TH2D("wave","wave", 520, 0, 520, 1024, 0, 1024);
  TH1D* wave_th1 = new TH1D("wave1d", "wave1d", 1024, 0, 1024);

  int compteurth1 = 0;
  float compteur[520] = {0};
  for (int i = 0; i < tree->GetEntries()/10; i++) {
    tree->GetEntry(i);

    for (int j = 0; j < om_number->size(); j++) {
      if (flag_e_event->at(j) == 1 && flag_associated_nohit->at(j) == 1 && energy->at(j) > 0.7) {
        compteur[om_number->at(j)]++;
        for (int k = 0; k < 1024; k++) {
          wave->Fill(om_number->at(j), k, waveform->at(j).at(k));
        }
      }
    }
  }

  float moyenne[520];
  memset(moyenne,0,sizeof(float)*520);
  for (int i = 0; i < 520; i++) {
    for (int k = 0; k < 1024; k++) {
      wave->SetBinContent(i+1,k+1, wave->GetBinContent(i+1, k+1)/compteur[i]);
      wave->SetBinError(i+1,k+1,0);
      if (k>0 && k < 193) {
        moyenne[i]+=wave->GetBinContent(i+1, k+1);
      }

    }

    moyenne[i] = moyenne[i]/192.;
    if (moyenne[i] == 0) moyenne[i] =1;
    // wave->SetBinContent(i,1024,wave->GetBinContent(i, 1023));
    // wave->SetBinError(i,1024,0);
  }

  for (int i = 0; i < 520; i++) {
    for (int k = 0; k < 1024; k++) {

      wave->SetBinContent(i+1,k+1,moyenne[i]-wave->GetBinContent(i+1,k+1));
    }
  }

  float peak_amplitude;
  int om;
  TFile *newfile = new TFile("test.root", "RECREATE");
  TTree Result_tree("Result_tree","");
  Result_tree.Branch("peak_amplitude", &peak_amplitude);
  Result_tree.Branch("om", &om);

  TH2F* mean_wave = new TH2F("mean_wave","mean_wave",1024,0,1024,500,0,500);

  for (int i = 0; i < 520; i++) {
    for (int k = 0; k < 1024; k++) {
      mean_wave->SetBinContent(k+1, round(wave->GetBinContent(i+1,k+1)), mean_wave->GetBinContent(i+1,k+1) + 1);
    }
    om = i;
    peak_amplitude = wave->ProjectionY("test",i+1,i+1)->GetMaximum();
    Result_tree.Fill();
  }

  newfile->cd();
  Result_tree.Write();
  wave->Write();
  mean_wave->Write();
  newfile->Close();
}

void biplot_eres_wave_mean() {

  float peak_amplitude, eres;
  int om;
  TFile *file = new TFile("mean_waveform.root", "READ");
  TTree* tree = (TTree*)file->Get("Result_tree");
  tree->SetBranchStatus("*",0);
  tree->SetBranchStatus("peak_amplitude",1);
  tree->SetBranchAddress("peak_amplitude", &peak_amplitude);
  tree->SetBranchStatus("om",1);
  tree->SetBranchAddress("om", &om);

  float amplitude[520];
  for (int i = 0; i < 520; i++) {
    tree->GetEntry(i);
    amplitude[om] = peak_amplitude;
  }

  TFile *file2 = new TFile("ERES.root", "READ");
  TTree* tree2 = (TTree*)file2->Get("Result_tree");
  tree2->SetBranchStatus("*",0);
  tree2->SetBranchStatus("eres",1);
  tree2->SetBranchAddress("eres", &eres);
  tree2->SetBranchStatus("om",1);
  tree2->SetBranchAddress("om", &om);

  float eres_tab[520];
  for (int i = 0; i < 520; i++) {
    tree2->GetEntry(i);
    eres_tab[om] = eres;
  }

  TH2F* biplot = new TH2F("mean_wave_eres", "mean_wave_eres", 150,250,400,150,10,25);
  gROOT->cd();
  for (size_t i = 0; i < 520; i++) {
    if (eres_tab[i] > 0 && amplitude[i] > 0) {
      biplot->SetBinContent(amplitude[i]-250, (eres_tab[i]-10)*10, biplot->GetBinContent(amplitude[i+1], eres_tab[i+1])+1);
    }
  }


  TFile *newfile = new TFile("test.root", "RECREATE");

  newfile->cd();
  biplot->Write();
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
  if (full == 1) {
    tri(run);
    R_calculator(run);
  }
  track_cutter(run);


  // waveformer();
  // biplot_eres_wave_mean();
  return 0;
}
