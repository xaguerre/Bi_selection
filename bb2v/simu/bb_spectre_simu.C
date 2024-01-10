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
double eres_corrector[520];

int column_verificator(int *column) {
  int col = 0;
  if (abs(column[1]-column[0]) < 2) col = 1;
  if (abs(column[1]-column[2]) < 3) col = 1;
  if (abs(column[3]-column[0]) < 3) col = 1;
  if (abs(column[2]-column[3]) < 4) col = 1;
  return col;
}

void energy_corrector_filler() {
  TFile tree_file("../../Simu/Bi_fit/fitted_bi_good_gas_edep.root", "READ");
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

void eres_corrector_filler() {
  TFile tree_file("../../Simu/Bi_fit/fitted_bi_good_gas_ecor.root", "READ");
  double eres_1MeV;
  int om_number;
  TTree* tree = (TTree*)tree_file.Get("Result_tree");
  tree->SetBranchStatus("*",0);
  tree->SetBranchStatus("eres_1MeV",1);
  tree->SetBranchAddress("eres_1MeV", &eres_1MeV);
  tree->SetBranchStatus("om_number",1);
  tree->SetBranchAddress("om_number", &om_number);

  for (int i = 0; i < tree->GetEntries(); i++) {
    tree->GetEntry(i);
    eres_corrector[om_number] = eres_1MeV;
  }
  return;
}

int om_verificator(int om1, int om2){
  if (abs(om1/13 - om2/13) > 1 ) return 1;
  return 0;
}

int side_verificator(int om1, int om2){
  if (om1/260 == om2/260) return 1;
  return 0;
}

void spectre_bb(string infile, string outfile) {
  energy_corrector_filler();
  eres_corrector_filler();
  gStyle->SetOptStat(1110);
  gStyle->SetOptFit(1);

  int om;
  double mean, sigma, eres_1MeV;

  TFile eres_file("../../Bi_fit/fit_histo/fitted_bi_974_edep_AI.root", "READ");
  TTree* eres_tree = (TTree*)eres_file.Get("Result_tree");
  eres_tree->SetBranchStatus("*",0);
  eres_tree->SetBranchStatus("eres_1MeV",1);
  eres_tree->SetBranchAddress("eres_1MeV", &eres_1MeV);
  eres_tree->SetBranchStatus("om_number",1);
  eres_tree->SetBranchAddress("om_number", &om);

  double FWHM_data[520];
  memset(FWHM_data, 0, 520*sizeof(double));
  for (int i = 0; i < eres_tree->GetEntries(); i++) {
    eres_tree->GetEntry(i);
    FWHM_data[om] = eres_1MeV;
  }
  TFile newfile(Form("%s_spectre_new.root", outfile.c_str()),"RECREATE");
  // TFile newfile("Tl_spectre.root","RECREATE");
  int ncalo_tot, ngg_tot, om_number_int, bdf, eventnumber;
  double p0, p1, p2, Chi2, mean_error, sigma_error;

  std::vector<int> *om_number = new std::vector<int>;
  std::vector<double> *energy = new std::vector<double>;
  std::vector<double> *energyvis = new std::vector<double>;
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

  TFile tree_file("cut_Tl_208_10G2.root", "READ");
  // TFile tree_file(Form("cut_%s.root", infile.c_str()), "READ");
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
  tree->SetBranchStatus("energy_ubc",1);
  tree->SetBranchAddress("energy_ubc", &energy);
  tree->SetBranchStatus("energyvis_ubc",1);
  tree->SetBranchAddress("energyvis_ubc", &energyvis);
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
  gROOT->cd();

  TH1D* spectre_back_to_back = new TH1D("spectre_back_to_back", "spectre_back_to_back", 400,0,10);
  TH1D* spectre_same_side = new TH1D("spectre_same_side", "spectre_same_side", 400,0,10);
  TH1D* spectre_full = new TH1D("spectrefull", "spectrefull", 400,0,10);
  TH1D* spectre_full_energy = new TH1D("spectrefull_energy", "spectrefull_energy", 200,0,10);
  TH1D* spectre_full_energy_intern = new TH1D("spectrefull_energy_intern", "spectrefull_energy_intern", 200,0,10);
  TH1D* spectre_full_energy_extern = new TH1D("spectrefull_energy_extern", "spectrefull_energy_extern", 200,0,10);

  TH1D* spectre_full_energyvis= new TH1D("spectrefull_energyvis", "spectrefull_energyvis", 400,0,10);

  int column[4];
  TRandom3 rando;

  int compteur2 = 0;
  TH1D* distrib = new TH1D("distrib", "distrib", 10,0,10);
  TH2D* Correlation = new TH2D("Correlation", "Correlation", 200,0,2,200,0,2);
  int compteur = 0;
  int intj = 0;
   int compteur_side =0;
  int compteur_back = 0;
  TH1D* distrib_temps = new TH1D("distrib_temps", "distrib_temps", 1000,-100,100);


  // for (int i = 0; i < 1; i++) {
  for (int i = 0; i < tree->GetEntries(); i++) {
    tree->GetEntry(i);
    std::vector<int> *event = new std::vector<int>;
    memset(column,0,4*sizeof(int));
    for (int j = 0; j < om_number->size(); j++) {
      if (flag_calo_square->at(j) == 1 && flag_charge ->at(j) == 1 && flag_MW->at(j) == 1 && flag_source_square->at(j) == 1 && flag_last_column->at(j) == 1 && energy->at(j)>0.2) {
        // if (flag_e_event->at(j) == 1 && energy->at(j)>0.2) {
        compteur++;
        event->push_back(j);
      }
    }
    if (event->size() == 2) {
      for (int eventi = 0; eventi < event->size(); eventi++) {
        for (int eventj = 1+eventi; eventj < event->size(); eventj++) {
          distrib_temps->Fill(calo_timestamp->at(event->at(eventi)) - calo_timestamp->at(event->at(eventj)));
            column[0] = first_column->at(eventi);
            if (column[0] < 0) column[0] = -1000;
            column[1] = first_column->at(eventj);
            if (column[1] < 0) column[0] = -2000;
            column[2] = second_column->at(eventi);
            if (column[2] < 0) column[0] = -3000;
            column[3] = second_column->at(eventj);
            if (column[3] < 0) column[0] = -4000;
            if (column_verificator(column) ==1 ) {
              double eresenergyi = rando.Gaus(energy->at(event->at(eventi)), ((FWHM_data[om_number->at(event->at(eventi))])/235.482)*sqrt(energy->at(event->at(eventi))))/sqrt(energy_corrector[om_number->at(eventj)]);
              double eresenergyj = rando.Gaus(energy->at(event->at(eventj)), ((FWHM_data[om_number->at(event->at(eventj))])/235.482)*sqrt(energy->at(event->at(eventi))))/sqrt(energy_corrector[om_number->at(eventj)]);
              if (abs(calo_timestamp->at(event->at(eventi)) - calo_timestamp->at(event->at(eventj))) < 2) {
                spectre_full_energy_intern->Fill(eresenergyj + eresenergyi);
                if (eresenergyj + eresenergyi > 3) {
                  cout << calo_timestamp->at(event->at(eventi)) - calo_timestamp->at(event->at(eventj)) << endl;
                }
              }
              else spectre_full_energy_extern->Fill(eresenergyj + eresenergyi);
              spectre_full_energy->Fill(eresenergyj + eresenergyi);
              // cout <<  i << endl;
              // energyvis_ubc->at(j) = rando.Gaus(energyvis_ubc->at(j), (FWHM_data[om_number->at(j)]-res_corrector[om_number->at(j)])/(235.482*FWHM_Bordeaux[om_number->at(j)])*sqrt(energyvis_ubc->at(j)));

              // if (side_verificator(om_number->at(event->at(eventi)), om_number->at(event->at(eventj))) == 1) {
              //   if (om_verificator(om_number->at(event->at(eventi)), om_number->at(event->at(eventj))) == 1) {
              //     spectre_same_side->Fill(eresenergyj + eresenergyi);
              //     spectre_full->Fill(eresenergyi);
              //     spectre_full->Fill(eresenergyj);
              //     spectre_full_energy->Fill(eresenergyj + eresenergyi);
              //     Correlation->Fill(eresenergyi,eresenergyj);
              //   }
              // }
              // else {
              //   spectre_back_to_back->Fill(eresenergyj + eresenergyi);
              //   spectre_full->Fill(eresenergyi);
              //   spectre_full->Fill(eresenergyj);
              //   spectre_full_energy->Fill(eresenergyj + eresenergyi);
              //   Correlation->Fill(eresenergyi,eresenergyj);
              // }


              // double eresenergyi = rando.Gaus(energy->at(event->at(eventi)), ((FWHM_data[om_number->at(event->at(eventi))])/235.482)*sqrt(energy->at(event->at(eventi))))/sqrt(energy_corrector[om_number->at(eventj)]);
              // double eresenergyj = rando.Gaus(energy->at(event->at(eventj)), ((FWHM_data[om_number->at(event->at(eventj))])/235.482)*sqrt(energy->at(event->at(eventi))))/sqrt(energy_corrector[om_number->at(eventj)]);
              // spectre_full->Fill(energy->at(event->at(eventi)));
              // spectre_full->Fill(energy->at(event->at(eventj)));
              // spectre_full_energy->Fill(eresenergyi + eresenergyj);
              // // spectre_full_energyvis->Fill(energyvis->at(event->at(eventj))/energy_corrector[om_number->at(eventj)] + energyvis->at(event->at(eventi))/energy_corrector[om_number->at(eventi)] );
              // Correlation->Fill(energy->at(event->at(eventi)),energy->at(event->at(eventj)));

          }
        }
      }
    }
    // Correlation->Fill(ene[0],ene[1]);
    compteur = 0;
    compteur2 = 0;
    intj = 0;

    delete event;

  }

  // cout << spectre_full_energy->GetEntries() << " : compteur " << compteur_back  + compteur_side << endl;
  // cout << spectre_back_to_back->GetEntries() << " : compteur " << compteur_back<< endl;
  // cout << spectre_same_side->GetEntries() << " : compteur " << compteur_side << endl;

  // spectre_full->Draw();
  // spectre_full_energy->Draw();
  // spectre_back_to_back->Draw("sames");
  // spectre_same_side->Draw("sames");
  //
  // spectre_back_to_back->SetLineColor(kRed);
  // spectre_same_side->SetLineColor(kGreen+8);

  TCanvas *c = new TCanvas();
  spectre_full_energy->Draw();
  spectre_full_energy->SetLineColor(kBlack);
  spectre_full_energy->SetLineWidth(2);
  spectre_full_energy_intern->Draw("sames");
  spectre_full_energy_intern->SetLineColor(kBlue);
  spectre_full_energy->SetLineWidth(2);
  spectre_full_energy_extern->Draw("sames");
  spectre_full_energy_extern->SetLineColor(kRed);
  spectre_full_energy->SetLineWidth(2);
  return;
  TLegend *legend = new TLegend(0.5,0.7,0.7,0.9);
  legend->AddEntry(spectre_full_energy,"spectre complet");
  legend->AddEntry(spectre_same_side,"spectre evenements meme cotes");
  legend->AddEntry(spectre_back_to_back,"spectre evenements cotes opposes");
  legend->Draw("same");
  TFile newfile2("spectre_bb2_test.root","RECREATE");
  newfile2.cd();
  spectre_full_energy->Write();
  spectre_full->Write();
  newfile2.Close();
}

void looper(/* arguments */) {
  spectre_bb("Tl_208_10G", "Tl");
  spectre_bb("Bi_214_20G", "Bi");
  spectre_bb("K_40_50G", "K");
  spectre_bb("bb_simu", "BB2");
}


void eff_gamma(string infile, string outfile) {
  energy_corrector_filler();
  eres_corrector_filler();
  gStyle->SetOptStat(1110);
  gStyle->SetOptFit(1);

  int om;
  double mean, sigma, eres_1MeV;

  TFile eres_file("../../Bi_fit/fit_histo/fitted_bi_974_edep_AI.root", "READ");
  TTree* eres_tree = (TTree*)eres_file.Get("Result_tree");
  eres_tree->SetBranchStatus("*",0);
  eres_tree->SetBranchStatus("eres_1MeV",1);
  eres_tree->SetBranchAddress("eres_1MeV", &eres_1MeV);
  eres_tree->SetBranchStatus("om_number",1);
  eres_tree->SetBranchAddress("om_number", &om);

  double FWHM_data[520];
  memset(FWHM_data, 0, 520*sizeof(double));
  for (int i = 0; i < eres_tree->GetEntries(); i++) {
    eres_tree->GetEntry(i);
    FWHM_data[om] = eres_1MeV;
  }
  // TFile newfile("Tl_spectre.root","RECREATE");
  int ncalo_tot, ngg_tot, om_number_int, bdf, eventnumber;
  double p0, p1, p2, Chi2, mean_error, sigma_error;

  std::vector<int> *om_number = new std::vector<int>;
  std::vector<double> *energy = new std::vector<double>;
  std::vector<double> *energyvis = new std::vector<double>;
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

  // TFile tree_file("cut_Tl_208_10G.root", "READ");
  TFile tree_file(Form("cut_%s.root", infile.c_str()), "READ");
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
  tree->SetBranchStatus("energy_ubc",1);
  tree->SetBranchAddress("energy_ubc", &energy);
  tree->SetBranchStatus("energyvis_ubc",1);
  tree->SetBranchAddress("energyvis_ubc", &energyvis);
  tree->SetBranchStatus("z_last_gg",1);
  tree->SetBranchAddress("z_last_gg", &z_last_gg);
  tree->SetBranchStatus("last_column",1);
  tree->SetBranchAddress("last_column", &last_column);
  tree->SetBranchStatus("first_column",1);
  tree->SetBranchAddress("first_column", &first_column);
  tree->SetBranchStatus("second_column",1);
  tree->SetBranchAddress("second_column", &second_column);
  gROOT->cd();


  TH2D* spectre_full = new TH2D("spectrefull", "spectrefull", 712,0,712, 1000,0,10);

  int column[4];
  TRandom3 rando;

  int compteur2 = 0;

  int compteur = 0;
  int intj = 0;

  // for (int i = 0; i < 1; i++) {
  for (int i = 0; i < tree->GetEntries(); i++) {
    tree->GetEntry(i);
    for (int j = 0; j < energyvis->size(); j++) {
      spectre_full->Fill(om_number->at(j), energy->at(j));
    }

  }
  TFile newfile("test","RECREATE");

  double eff_tot, eff_part;
  int om_num;

  TTree Result_tree4("Result_tree","");
  Result_tree4.Branch("eff_tot", &eff_tot);
  Result_tree4.Branch("eff_part", &eff_part);
  Result_tree4.Branch("om_num", &om_num);

  for (size_t i = 0; i < 712; i++) {
    om_num = i;
    eff_tot = spectre_full->ProjectionY("",i+1,i+1)->Integral()/1e10;
    eff_part = spectre_full->ProjectionY("",i+1,i+1)->Integral(200,1000)/1e10;
    Result_tree4.Fill();
  }




  spectre_full->Draw();
  newfile.cd();
  Result_tree4.Write();
  newfile.Close();

}
