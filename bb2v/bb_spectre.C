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
// #include "sndisplay.cc

double energy_corrector[520];
double time_corrector[520];


int column_verificator(int *column) {
  int col = 0;
  if (abs(column[1]-column[0]) < 2) col = 1;
  if (abs(column[1]-column[2]) < 3) col = 1;
  if (abs(column[3]-column[0]) < 3) col = 1;
  if (abs(column[2]-column[3]) < 4) col = 1;
  return col;
}

int om_verificator(int om1, int om2){
  if (abs(om1/13 - om2/13) > 1 ) return 1;
  return 0;
}

int side_verificator(int om1, int om2){
  if (om1/260 == om2/260) return 1;
  return 0;
}

void energy_corrector_filler() {
  TFile tree_file("../Simu/Bi_fit/fitted_bi_good_gas_edep.root", "READ");
  double mean;
  int om_number;
  TTree* tree = (TTree*)tree_file.Get("Result_tree");
  tree->SetBranchStatus("*",0);
  tree->SetBranchStatus("mean",1);
  tree->SetBranchAddress("mean", &mean);
  tree->SetBranchStatus("om_number",1);
  tree->SetBranchAddress("om_number", &om_number);

  for (int i = 0; i < tree->GetEntries(); i++) {
    tree->GetEntry(i);
    energy_corrector[om_number] = mean;
  }
  return;
}

void time_filler() {
  memset (time_corrector, 0, 520*sizeof(double));
  double err_time;
  std::ifstream tempsit("tempsit.txt");
  int om_num;
  while (tempsit >> om_num)
  {
    tempsit >> time_corrector[om_num];
    tempsit >> err_time;
  }
  std::ifstream tempsfr("tempsfr.txt");
  while (tempsfr >> om_num)
  {
    tempsfr >> time_corrector[om_num];
    tempsfr >> err_time;
  }
  for (int om_number = 0; om_number < 520; om_number++) {
    if (om_number == 210 || om_number == 192 || om_number == 258 || om_number == 293 || om_number == 9 || om_number == 34 || om_number == 512) {
      time_corrector[om_number] = 0;
    }
  }
  return;
}

void spectre_bb(int run) {
  time_filler();
  gStyle->SetOptStat(1);
  gStyle->SetOptFit(1);
  energy_corrector_filler();
  int om;
  double mean, sigma;

  TFile *file3 = new TFile("../Simu/Bi_fit/fitted_bi_gas.root", "READ");
  TTree* tree3 = (TTree*)file3->Get("Result_tree");
  tree3->SetBranchStatus("*",0);
  tree3->SetBranchStatus("om_number",1);
  tree3->SetBranchAddress("om_number", &om);
  tree3->SetBranchStatus("mean",1);
  tree3->SetBranchAddress("mean", &mean);
  tree3->SetBranchStatus("sigma",1);
  tree3->SetBranchAddress("sigma", &sigma);


  double mean_tab[520];
  double sigma_tab[520];
  memset (mean_tab, 0, 520*sizeof(double));
  memset (sigma_tab, 0, 520*sizeof(double));

  for (int i = 0; i < tree3->GetEntries(); i++) {
    tree3->GetEntry(i);
    if (mean > 0) {
      mean_tab[om] = mean;
      sigma_tab[om] = sigma;
    }
  }

  TFile *file5 = new TFile("../Bi_fit/fitted_bi_1063.root", "READ");
  TTree* tree5 = (TTree*)file5->Get("Result_tree");
  tree5->SetBranchStatus("*",0);
  tree5->SetBranchStatus("om_number",1);
  tree5->SetBranchAddress("om_number", &om);
  tree5->SetBranchStatus("mean",1);
  tree5->SetBranchAddress("mean", &mean);
  tree5->SetBranchStatus("sigma",1);
  tree5->SetBranchAddress("sigma", &sigma);


  double mean_tab_1059[712];
  double sigma_tab_1059[712];
  memset (mean_tab_1059, 0, 712*sizeof(double));
  memset (sigma_tab_1059, 0, 712*sizeof(double));

  for (int i = 0; i < tree3->GetEntries(); i++) {
    tree5->GetEntry(i);
    if (mean > 0) {
      mean_tab_1059[om] = mean;
      sigma_tab_1059[om] = sigma;
    }
  }

  int ncalo_tot, ngg_tot, om_number_int, bdf, eventnumber;
  double p0, p1, p2, Chi2, mean_error, sigma_error;

  std::vector<int> *om_number = new std::vector<int>;
  std::vector<double> *energy = new std::vector<double>;
  std::vector<int> *flag_e_event = new std::vector<int>;
  std::vector<int> *flag_calo_square = new std::vector<int>;
  std::vector<int> *flag_source_square = new std::vector<int>;
  std::vector<int> *flag_last_column = new std::vector<int>;
  std::vector<int> *flag_MW = new std::vector<int>;
  std::vector<int> *flag_charge = new std::vector<int>;
  std::vector<int> *flag_last_z = new std::vector<int>;
  std::vector<int> *gg_associated = new std::vector<int>;
  std::vector<double> *z_last_gg = new std::vector<double>;
  std::vector<double> *last_column = new std::vector<double>;
  std::vector<double> *first_column = new std::vector<double>;
  std::vector<double> *second_column = new std::vector<double>;
  std::vector<int> *calo_nohit_om_time = new std::vector<int>;
  std::vector<int> *rising_cell = new std::vector<int>;
  std::vector<double> *calo_timestamp = new std::vector<double>;
  std::vector<int> *charge = new std::vector<int>;
  std::vector<int> *falling_cell = new std::vector<int>;

  TFile tree_file(Form("cut_%d.root",run), "READ");
  TTree* tree = (TTree*)tree_file.Get("Result_tree");
  tree->SetBranchStatus("*",0);
  tree->SetBranchStatus("eventnumber",1);
  tree->SetBranchAddress("eventnumber", &eventnumber);
  tree->SetBranchStatus("om_number",1);
  tree->SetBranchAddress("om_number", &om_number);
  tree->SetBranchStatus("flag_e_event",1);
  tree->SetBranchAddress("flag_e_event", &flag_e_event);
  tree->SetBranchStatus("flag_calo_square",1);
  tree->SetBranchAddress("flag_calo_square", &flag_calo_square);
  tree->SetBranchStatus("flag_last_column",1);
  tree->SetBranchAddress("flag_last_column", &flag_last_column);
  tree->SetBranchStatus("flag_source_square",1);
  tree->SetBranchAddress("flag_source_square", &flag_source_square);
  tree->SetBranchStatus("flag_last_z",1);
  tree->SetBranchAddress("flag_last_z", &flag_last_z);
  tree->SetBranchStatus("flag_MW",1);
  tree->SetBranchAddress("flag_MW", &flag_MW);
  tree->SetBranchStatus("flag_charge",1);
  tree->SetBranchAddress("flag_charge", &flag_charge);
  tree->SetBranchStatus("energy",1);
  tree->SetBranchAddress("energy", &energy);
  tree->SetBranchStatus("charge",1);
  tree->SetBranchAddress("charge", &charge);
  tree->SetBranchStatus("ncalo_tot",1);
  tree->SetBranchAddress("ncalo_tot", &ncalo_tot);
  tree->SetBranchStatus("ngg_tot",1);
  tree->SetBranchAddress("ngg_tot", &ngg_tot);
  tree->SetBranchStatus("z_last_gg",1);
  tree->SetBranchAddress("z_last_gg", &z_last_gg);
  tree->SetBranchStatus("last_column",1);
  tree->SetBranchAddress("last_column", &last_column);
  tree->SetBranchStatus("first_column",1);
  tree->SetBranchAddress("first_column", &first_column);
  tree->SetBranchStatus("second_column",1);
  tree->SetBranchAddress("second_column", &second_column);
  tree->SetBranchStatus("calo_timestamp",1);
  tree->SetBranchAddress("calo_timestamp", &calo_timestamp);
  tree->SetBranchStatus("calo_nohit_om_time",1);
  tree->SetBranchAddress("calo_nohit_om_time", &calo_nohit_om_time);
  tree->SetBranchStatus("digicalo.falling_cell",1);
  tree->SetBranchAddress("digicalo.falling_cell", &falling_cell);
  gROOT->cd();

  TH1D* spectre_back_to_back = new TH1D("spectre_back_to_back", "spectre_back_to_back", 400,0,10);
  TH1D* spectre_same_side = new TH1D("spectre_same_side", "spectre_same_side", 400,0,10);
  TH1D* spectre_full = new TH1D("spectrefull", "spectrefull", 300,0,10);
  TH1D* spectre_full_energy = new TH1D("spectrefull_energy", "spectrefull_energy", 400,0,10);
  TH1D* spectre_temps = new TH1D("spectretemps", "spectretemps", 200,-50,50);
  TH1D* spectre_temps2 = new TH1D("spectretemps2", "spectretemps2", 200,-50,50);

  int column[4];

  int compteur2 = 0;
  TH1D* distrib = new TH1D("distrib", "distrib", 10,0,10);
  TH2D* Correlation = new TH2D("Correlation", "Correlation", 200,0,2,200,0,2);
  int compteur = 0;
  int intj = 0;
  int compteur3 = 0;

  int compteur_side =0;
  int compteur_back = 0;

  // for (int i = 0; i < 1; i++) {
  for (int i = 0; i < tree->GetEntries(); i++) {
    tree->GetEntry(i);
    std::vector<int> *event = new std::vector<int>;
    std::vector<int> *eventt = new std::vector<int>;
    memset(column,0,4*sizeof(int));

    for (int j = 0; j < om_number->size(); j++) {
      if (flag_calo_square->at(j) == 1 && flag_charge ->at(j) == 1 && flag_MW->at(j) == 1 && flag_source_square->at(j) == 1 && flag_last_column->at(j) == 1 && energy->at(j)>0.2) {
        compteur++;
        event->push_back(j);
      }
      else if (flag_calo_square->at(j) == 1 && flag_charge ->at(j) == 1 && flag_MW->at(j) == 1 && flag_source_square->at(j) == 1 && flag_last_column->at(j) == 1){
        eventt->push_back(j);
        compteur3++;
      }
    }
    if (event->size() == 2) {
      for (int eventi = 0; eventi < event->size(); eventi++) {
        for (int eventj = 1+eventi; eventj < event->size(); eventj++) {
          // if (abs(calo_timestamp->at(event->at(eventi)) - calo_timestamp->at(event->at(eventj))) < 50) {
          column[0] = first_column->at(eventi);
          if (column[0] < 0) column[0] = -1000;
          column[1] = first_column->at(eventj);
          if (column[1] < 0) column[0] = -2000;
          column[2] = second_column->at(eventi);
          if (column[2] < 0) column[0] = -3000;
          column[3] = second_column->at(eventj);
          if (column[3] < 0) column[0] = -4000;

          if (column_verificator(column) ==1 && abs(calo_timestamp->at(event->at(eventi)) + falling_cell->at(event->at(eventi)) - time_corrector[om_number->at(event->at(eventi))] - calo_timestamp->at(event->at(eventj)) - falling_cell->at(event->at(eventj)) + time_corrector[om_number->at(event->at(eventj))]) < 50) {
            spectre_temps->Fill(calo_timestamp->at(event->at(eventi)) + falling_cell->at(event->at(eventi)) - time_corrector[om_number->at(event->at(eventi))] - calo_timestamp->at(event->at(eventj)) - falling_cell->at(event->at(eventj)) + time_corrector[om_number->at(event->at(eventj))]);
            // cout << abs(calo_timestamp->at(event->at(eventi)) - time_corrector[om_number->at(event->at(eventi))] - calo_timestamp->at(event->at(eventj)) - time_corrector[om_number->at(event->at(eventj))])<< endl;

            // spectre_temps2->Fill(abs(calo_timestamp->at(event->at(eventi)) + time_corrector[om_number->at(event->at(eventi))] - calo_timestamp->at(event->at(eventj))) + time_corrector[om_number->at(event->at(eventj))]);

            if (side_verificator(om_number->at(event->at(eventi)), om_number->at(event->at(eventj))) == 1) {
              if (om_verificator(om_number->at(event->at(eventi)), om_number->at(event->at(eventj))) == 1) {
                compteur_side++;
                spectre_same_side->Fill(charge->at(event->at(eventj))/mean_tab_1059[om_number->at(event->at(eventj))]*energy_corrector[om_number->at(event->at(eventj))] + charge->at(event->at(eventi))/mean_tab_1059[om_number->at(event->at(eventi))]*energy_corrector[om_number->at(event->at(eventi))]);
                spectre_full->Fill(charge->at(event->at(eventi))/mean_tab_1059[om_number->at(event->at(eventi))]);
                spectre_full->Fill(charge->at(event->at(eventj))/mean_tab_1059[om_number->at(event->at(eventj))]);
                spectre_full_energy->Fill(charge->at(event->at(eventj))/mean_tab_1059[om_number->at(event->at(eventj))]*energy_corrector[om_number->at(event->at(eventj))] + charge->at(event->at(eventi))/mean_tab_1059[om_number->at(event->at(eventi))]*energy_corrector[om_number->at(event->at(eventi))]);
                Correlation->Fill(charge->at(event->at(eventi))/mean_tab_1059[om_number->at(event->at(eventi))],charge->at(event->at(eventj))/mean_tab_1059[om_number->at(event->at(eventj))]);
                if (charge->at(event->at(eventj))/mean_tab_1059[om_number->at(event->at(eventj))]*energy_corrector[om_number->at(event->at(eventj))] + charge->at(event->at(eventi))/mean_tab_1059[om_number->at(event->at(eventi))]*energy_corrector[om_number->at(event->at(eventi))] > 2.8)
                  if (charge->at(event->at(eventj))/mean_tab_1059[om_number->at(event->at(eventj))]*energy_corrector[om_number->at(event->at(eventj))] + charge->at(event->at(eventi))/mean_tab_1059[om_number->at(event->at(eventi))]*energy_corrector[om_number->at(event->at(eventi))] < 3.2){
                    cout <<"energy" << charge->at(event->at(eventj))/mean_tab_1059[om_number->at(event->at(eventj))]*energy_corrector[om_number->at(event->at(eventj))] + charge->at(event->at(eventi))/mean_tab_1059[om_number->at(event->at(eventi))]*energy_corrector[om_number->at(event->at(eventi))] << " -> event number " << eventnumber << endl;
                    cout << "om : " << om_number->at(event->at(eventi)) << " : " << charge->at(event->at(eventj))/mean_tab_1059[om_number->at(event->at(eventj))]*energy_corrector[om_number->at(event->at(eventj))]<< " MeV et " << om_number->at(event->at(eventj)) << " : " << charge->at(event->at(eventi))/mean_tab_1059[om_number->at(event->at(eventi))]*energy_corrector[om_number->at(event->at(eventi))] << " MeV"<< endl;
                    for (int h = 0; h < eventt->size(); h++) {
                      if (abs(calo_timestamp->at(eventt->at(h)) - calo_timestamp->at(event->at(eventj))) < 50) {
                        cout << "om + : " << om_number->at(eventt->at(h)) << " : " << charge->at(eventt->at(h))/mean_tab_1059[om_number->at(eventt->at(h))]*energy_corrector[om_number->at(eventt->at(h))] << " MeV" << endl;
                      }
                    }
                  }
              }
            }
            else {
              compteur_back++;
              spectre_back_to_back->Fill(charge->at(event->at(eventj))/mean_tab_1059[om_number->at(event->at(eventj))]*energy_corrector[om_number->at(event->at(eventj))] + charge->at(event->at(eventi))/mean_tab_1059[om_number->at(event->at(eventi))]*energy_corrector[om_number->at(event->at(eventi))]);
              spectre_full->Fill(charge->at(event->at(eventi))/mean_tab_1059[om_number->at(event->at(eventi))]);
              spectre_full->Fill(charge->at(event->at(eventj))/mean_tab_1059[om_number->at(event->at(eventj))]);
              spectre_full_energy->Fill(charge->at(event->at(eventj))/mean_tab_1059[om_number->at(event->at(eventj))]*energy_corrector[om_number->at(event->at(eventj))] + charge->at(event->at(eventi))/mean_tab_1059[om_number->at(event->at(eventi))]*energy_corrector[om_number->at(event->at(eventi))]);
              Correlation->Fill(charge->at(event->at(eventi))/mean_tab_1059[om_number->at(event->at(eventi))],charge->at(event->at(eventj))/mean_tab_1059[om_number->at(event->at(eventj))]);
              if (charge->at(event->at(eventj))/mean_tab_1059[om_number->at(event->at(eventj))]*energy_corrector[om_number->at(event->at(eventj))] + charge->at(event->at(eventi))/mean_tab_1059[om_number->at(event->at(eventi))]*energy_corrector[om_number->at(event->at(eventi))] > 2.8)
                if (charge->at(event->at(eventj))/mean_tab_1059[om_number->at(event->at(eventj))]*energy_corrector[om_number->at(event->at(eventj))] + charge->at(event->at(eventi))/mean_tab_1059[om_number->at(event->at(eventi))]*energy_corrector[om_number->at(event->at(eventi))] < 3.2){
                cout <<"energy" << charge->at(event->at(eventj))/mean_tab_1059[om_number->at(event->at(eventj))]*energy_corrector[om_number->at(event->at(eventj))] + charge->at(event->at(eventi))/mean_tab_1059[om_number->at(event->at(eventi))]*energy_corrector[om_number->at(event->at(eventi))] << " -> event number " << eventnumber << endl;
                cout << "om : " << om_number->at(event->at(eventi)) << " : " << charge->at(event->at(eventj))/mean_tab_1059[om_number->at(event->at(eventj))]*energy_corrector[om_number->at(event->at(eventj))]<< "MeV et " << om_number->at(event->at(eventj)) << " : " << charge->at(event->at(eventi))/mean_tab_1059[om_number->at(event->at(eventi))]*energy_corrector[om_number->at(event->at(eventi))] << " MeV"<< endl;
                for (int h = 0; h < eventt->size(); h++) {
                  if (abs(calo_timestamp->at(eventt->at(h)) - calo_timestamp->at(event->at(eventj))) < 50) {
                    cout << "om + : " << om_number->at(eventt->at(h)) << " : " << charge->at(eventt->at(h))/mean_tab_1059[om_number->at(eventt->at(h))]*energy_corrector[om_number->at(eventt->at(h))] << " MeV" << endl;
                  }
                }
              }
            }

          // if (energy->at(eventi) < 0.75) {
            //   spectre_05->Fill(energy->at(eventi));
            // }
            // else if (energy->at(eventi) > 0.75) {
            //   spectre_1->Fill(energy->at(eventi));
            // }
          }
        }
      }
    }
    // Correlation->Fill(ene[0],ene[1]);
    compteur = 0;
    compteur2 = 0;
    compteur3 = 0;
    intj = 0;

    delete event;

  }
  cout << spectre_full_energy->GetEntries() << " : compteur " << compteur_back  + compteur_side << endl;
  cout << spectre_back_to_back->GetEntries() << " : compteur " << compteur_back<< endl;
  cout << spectre_same_side->GetEntries() << " : compteur " << compteur_side << endl;


  Correlation->Draw("colz");
  TCanvas *c = new TCanvas();
  // spectre_full->Draw();
  spectre_full_energy->Draw();
  spectre_back_to_back->Draw("same");
  spectre_same_side->Draw("same");

  spectre_back_to_back->SetLineColor(kRed);
  spectre_same_side->SetLineColor(kGreen+8);



  TLegend *legend = new TLegend(0.5,0.7,0.7,0.9);
  legend->AddEntry(spectre_full_energy,"spectre complet");
  legend->AddEntry(spectre_same_side,"spectre evenements meme cotes");
  legend->AddEntry(spectre_back_to_back,"spectre evenements cotes opposes");
  legend->Draw("same");
  TFile newfile("spectre_bb2_test.root","RECREATE");
  newfile.cd();
  spectre_full_energy->Write();
  spectre_full->Write();
  newfile.Close();

  TCanvas *d = new TCanvas();
  spectre_temps->Draw();
  // spectre_temps2->Draw("same");
}

//
// void tracer() {
//
//   int om;
//   double eres_1MeV;
//
//   double eres_tab[520];
//
//   TFile *file5 = new TFile("../Bi_fit/fit_histo/fitted_bi_974_edep_AI.root", "READ");
//   TTree* tree5 = (TTree*)file5->Get("Result_tree");
//   tree5->SetBranchStatus("*",0);
//   tree5->SetBranchStatus("om_number",1);
//   tree5->SetBranchAddress("om_number", &om);
//   tree5->SetBranchStatus("eres_1MeV",1);
//   tree5->SetBranchAddress("eres_1MeV", &eres_1MeV);
//
//   for (size_t i = 0; i < tree5->GetEntries(); i++) {
//     tree5->GetEntry(i);
//     eres_tab[om] = eres_1MeV;
//   }
//
//
//   TH1D* energyth = new TH1D("","",100,0.1,3.2);
//   TH2D* energyresth = new TH2D("","",100,5,25,100,5,25);
//
//   double energy_tab[60] = {0.584,2.271,2.031,0.768,2.665,0.343,0.955,1.993,2.367,0.586,2.625,0.467,0.995,2.010,1.432,1.622,0.718,2.133,2.527,0.300,1.294,1.587,1.293,1.625,1.347,1.619,1.528,1.336,0.742,2.174,0.742,2.180,0.762,2.265,1.161,1.939,2.058,0.747,0.855,2.201,1.007,1.934,0.767,2.172,1.474,1.592,1.347,1.469,1.156,1.696,1.710,1.100,1.920,1.007,0.907,1.940,2.461,0.351};
//   int om_tab[60] = {54,303,239,512,152,412,196,464,184,431,101,379,179,398,163,397,30,317,254,513,150,411,226,501,56,302,41,327,223,497,177,465,151,424,18,305,102,387,189,436,245,501,229,515,224,509,374,500,134,394,108,355,174,447,228,488,218,465};
//
//
//   for (int i = 0; i < 60; i++) {
//     energyth->Fill(energy_tab[i]);
//   }
//   for (int i = 0; i < 30; i++) {
//     energyresth->Fill(eres_tab[om_tab[i*2]],eres_tab[om_tab[i*2+1]]);
//     cout << i << " -> "<< eres_tab[om_tab[i*2]] << " : " << eres_tab[om_tab[i*2+1]] << endl;;
//   }
//
// energyresth->Draw();
//
// // energyth->Draw();
//
// }
