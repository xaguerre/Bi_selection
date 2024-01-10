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


void stored_eres() {
  std::ifstream infile("stored_eres.txt");
  double FWHM[712];
  memset(FWHM,0, 712*sizeof(double));
  if (infile) cout << "file ok " << endl;
  else cout <<"file not ok " << endl;
  int om;
  double eres;
  int compteur = 0;
  string line;
  while(getline(infile,line)){
     FWHM[compteur] = stod(line);
    cout << FWHM[compteur] << endl;;
    compteur++;
  }

  TH1D *distrib = new TH1D("distrib","distrib", 100,5,15);

  TFile *file = new TFile("stored_eres.root", "RECREATE");
  TTree Result_tree("Result_tree","");
  Result_tree.Branch("FWHM", &eres);
  Result_tree.Branch("om", &om);

  for (size_t i = 0; i < 520; i++) {
    distrib->Fill(FWHM[i]);
    eres = FWHM[i];
    om = i;
    Result_tree.Fill();
  }

  distrib->Draw();

  file->cd();
  Result_tree.Write();
  file->Close();
}

void biplot_bdf() {

  int om_number;
  float ERES;
  float eres[520];
  float eres_bdf[520];
  memset(eres,0,520*sizeof(float));
  memset(eres_bdf,0,520*sizeof(float));

  TFile tree_file("Simu/eres_1055.root", "READ");
  TTree* tree = (TTree*)tree_file.Get("Result_tree");
  tree->SetBranchStatus("*",0);
  tree->SetBranchStatus("om",1);
  tree->SetBranchAddress("om", &om_number);
  tree->SetBranchStatus("eres",1);
  tree->SetBranchAddress("eres", &ERES);

  for (int i = 0; i < tree->GetEntries(); i++) {
      tree->GetEntry(i);
      eres[om_number] = ERES;
    }

  TFile tree_file2("Simu/eres_1055_bdf.root", "READ");
  TTree* tree2 = (TTree*)tree_file2.Get("Result_tree");
  tree2->SetBranchStatus("*",0);
  tree2->SetBranchStatus("om",1);
  tree2->SetBranchAddress("om", &om_number);
  tree2->SetBranchStatus("eres",1);
  tree2->SetBranchAddress("eres", &ERES);

  for (int i = 0; i < tree2->GetEntries(); i++) {
    tree2->GetEntry(i);
    eres_bdf[om_number] = ERES;
  }
  gROOT->cd();
  TH2D *biplot = new TH2D("biplot", "biplot",100,5,25,100,5,25);

  for (int i = 0; i < 520; i++) {
    if (eres[i] >0 && eres_bdf[i] > 0) {
      if (abs(eres[i] - eres_bdf[i]) > 1.5) {
        cout << i << " eres " <<eres[i]<< " bdf " <<eres_bdf[i]<<endl;
      }
      biplot->Fill(eres[i], eres_bdf[i]);
    }
  }
  biplot->Draw();

}

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

void z_comparator() {
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
  std::vector<double> *z_last_gg = new std::vector<double>;

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
  tree->SetBranchStatus("z_last_gg",1);
  tree->SetBranchAddress("z_last_gg", &z_last_gg);

  TH1D *spectre_simu = new TH1D ("z_last_gg_simu", "", 500, -1, 1);
  // tree->Project("spectre_simu", "om_num", "flag_last_column == 1 && ncalo_tot < 2 && energyvis_ubc > 0.6", "", tree->GetEntries()/10);
  // tree->Project("spectre_simu", "om_num", "flag_charge == 1 && flag_associated_nohit == 1 && flag_MW == 1 && flag_calo_square == 1 && flag_source_square == 1 && flag_last_column == 1 && ncalo_tot < 2 && energyvis_ubc > 0.6", "", tree->GetEntries()/10);
  tree->Project("z_last_gg_simu", "z_last_gg", "", "", tree->GetEntries()/10);
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
  tree2->SetBranchStatus("z_last_gg",1);
  tree2->SetBranchAddress("z_last_gg", &z_last_gg);

  TH1D *spectre_data = new TH1D ("z_last_gg_data", "", 500, -1, 1);
  // tree2->Project("spectre_data", "om_number", "flag_last_column == 1 && ncalo_tot < 2 && charge/22000. > 0.6", "", tree2->GetEntries()/10);
  // tree2->Project("spectre_data", "om_number", "flag_charge == 1 && flag_associated_nohit == 1 && flag_MW == 1 && flag_calo_square == 1 && flag_source_square == 1 && flag_last_column == 1 && ncalo_tot < 2 && charge/22000. > 0.6", "", tree2->GetEntries()/10);
  tree2->Project("z_last_gg_data", "z_last_gg", "", "", tree2->GetEntries()/10);
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

  spectre_simu->Draw("pl");
  spectre_simu->SetFillColorAlpha(kBlue, 0.35);
  spectre_data->Scale(1/1.78);
  spectre_data->Draw("plsame");
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

Double_t ScaleX(Double_t x, int om)
{
  Double_t v;
  // v =  x*energy_corrector[om] + 0.021; // "linear scaling" function example   ECORR
  // v =  x*energy_corrector[om] + 0.038; // "linear scaling" function example   ECORR
  // cout << energy_corrector[om] << endl;
  // v =  x*energy_corrector[om]; // Passer de Qnorm à Edep SI *(energy_corrector)
  // v =  x*(energy_corrector[om]-0.013) + 0.013; // Passer de Qnorm à Edep SI *(energy_corrector) puis à Edep AI* (1-b/energy_correctore)
  v =  x*(1-0.061) + 0.061; // Passer de Qnorm à Edep SI *(energy_corrector) puis à Edep AI* (1-b/energy_correctore)
  // v = 0.8882 * x + 0.0400; // "linear scaling" function example      EDEP
  return v;
}

void ScaleAxis(TAxis *a, Double_t (*Scale)(Double_t, int), int om)
{
  if (!a) return; // just a precaution
  if (a->GetXbins()->GetSize())
    {
      // an axis with variable bins
      // note: bins must remain in increasing order, hence the "Scale"
      // function must be strictly (monotonically) increasing
      TArrayD X(*(a->GetXbins()));
      for(Int_t i = 0; i < X.GetSize(); i++) X[i] = Scale(X[i], om);
      a->Set((X.GetSize() - 1), X.GetArray()); // new Xbins
    }
  else
    {
      // an axis with fix bins
      // note: we modify Xmin and Xmax only, hence the "Scale" function
      // must be linear (and Xmax must remain greater than Xmin)
      a->Set( a->GetNbins(),
              Scale(a->GetXmin(), om), // new Xmin
              Scale(a->GetXmax(), om) ); // new Xmax
    }
  return;
}

void ScaleXaxis(TH1 *h, Double_t (*Scale)(Double_t, int), int om)
{
  if (!h) return; // just a precaution
  ScaleAxis(h->GetXaxis(), Scale, om);
  return;
}


void spectre_comparator(){


  double simu_mean[520];
  memset (simu_mean, 0, 520*sizeof(double));
  double correction_tab[520];
  memset (correction_tab, 0, 520*sizeof(double));
  int om_num;
  double mean, correction;

  TFile *file = new TFile("Simu/Bi_fit/ubc_corrector.root", "READ");
  TTree* tree = (TTree*)file->Get("Result_tree");
  tree->SetBranchStatus("*",0);
  tree->SetBranchStatus("om_number",1);
  tree->SetBranchAddress("om_number", &om_num);
  tree->SetBranchStatus("mean",1);
  tree->SetBranchAddress("mean", &mean);
  tree->SetBranchStatus("correction",1);
  tree->SetBranchAddress("correction", &correction);

  for (int i = 0; i < tree->GetEntries(); i++) {
    tree->GetEntry(i);
    simu_mean[om_num] = mean;
    correction_tab[om_num] = correction;
  }

  int ncalo_tot;
  std::vector<int> *flag_e_event = new std::vector<int>;
  std::vector<int> *om_number = new std::vector<int>;
  std::vector<double> *energyvis_ubc = new std::vector<double>;
  std::vector<double> *energyvis = new std::vector<double>;
  std::vector<double> *charge = new std::vector<double>;
  std::vector<int> *flag_charge = new std::vector<int>;
  std::vector<int> *flag_associated_nohit = new std::vector<int>;
  std::vector<int> *flag_MW = new std::vector<int>;
  std::vector<int> *flag_calo_square = new std::vector<int>;
  std::vector<int> *flag_source_square = new std::vector<int>;
  std::vector<int> *flag_last_column = new std::vector<int>;
  std::vector<int> *flag_last_z = new std::vector<int>;
  std::vector<int> *calo_nohit_om_time = new std::vector<int>;

  TFile *file_simu = new TFile("Simu/cut_bi_new.root", "READ");
  TTree* tree_simu = (TTree*)file_simu->Get("Result_tree");
  tree_simu->SetBranchStatus("*",0);
  tree_simu->SetBranchStatus("ncalo_tot",1);
  tree_simu->SetBranchAddress("ncalo_tot", &ncalo_tot);
  tree_simu->SetBranchStatus("om_num",1);
  tree_simu->SetBranchAddress("om_num", &om_number);
  tree_simu->SetBranchStatus("flag_e_event",1);
  tree_simu->SetBranchAddress("flag_e_event", &flag_e_event);
  tree_simu->SetBranchStatus("flag_charge",1);
  tree_simu->SetBranchAddress("flag_charge", &flag_charge);
  tree_simu->SetBranchStatus("flag_associated_nohit",1);
  tree_simu->SetBranchAddress("flag_associated_nohit", &flag_associated_nohit);
  tree_simu->SetBranchStatus("flag_MW",1);
  tree_simu->SetBranchAddress("flag_MW", &flag_MW);
  tree_simu->SetBranchStatus("flag_calo_square",1);
  tree_simu->SetBranchAddress("flag_calo_square", &flag_calo_square);
  tree_simu->SetBranchStatus("flag_source_square",1);
  tree_simu->SetBranchAddress("flag_source_square", &flag_source_square);
  tree_simu->SetBranchStatus("flag_last_column",1);
  tree_simu->SetBranchAddress("flag_last_column", &flag_last_column);
  tree_simu->SetBranchStatus("flag_last_z",1);
  tree_simu->SetBranchAddress("flag_last_z", &flag_last_z);
  tree_simu->SetBranchStatus("energyvis_ubc",1);
  tree_simu->SetBranchAddress("energyvis_ubc", &energyvis_ubc);
  tree_simu->SetBranchStatus("energyvis",1);
  tree_simu->SetBranchAddress("energyvis", &energyvis);

  TH2D* simu = new TH2D("simu","simu",520,0,520,200,0,2);
  TH2D* simu_ubc = new TH2D("simu_ubc","simu_ubc",520,0,520,200,0,2);
  TH2D* data = new TH2D("data","data",520,0,520,200,0,2);

  for (int i = 0; i < tree_simu->GetEntries(); i++) {
    tree_simu->GetEntry(i);
    for (int j = 0; j < flag_e_event->size(); j++){
      if (flag_e_event->at(j) == 1 && ncalo_tot < 2 ) {
        simu->Fill(om_number->at(j),energyvis->at(j));
        simu_ubc->Fill(om_number->at(j),energyvis_ubc->at(j)*correction_tab[om_number->at(j)]);
      }
    }
  }

  cout << "simu one" << endl;

  TFile *file_data = new TFile("cut_974.root", "READ");
  TTree* tree_data = (TTree*)file_data->Get("Result_tree");
  tree_data->SetBranchStatus("*",0);
  tree_data->SetBranchStatus("ncalo_tot",1);
  tree_data->SetBranchAddress("ncalo_tot", &ncalo_tot);
  tree_data->SetBranchStatus("om_number",1);
  tree_data->SetBranchAddress("om_number", &om_number);
  tree_data->SetBranchStatus("flag_e_event",1);
  tree_data->SetBranchAddress("flag_e_event", &flag_e_event);
  tree_data->SetBranchStatus("flag_charge",1);
  tree_data->SetBranchAddress("flag_charge", &flag_charge);
  tree_data->SetBranchStatus("flag_associated_nohit",1);
  tree_data->SetBranchAddress("flag_associated_nohit", &flag_associated_nohit);
  tree_data->SetBranchStatus("flag_MW",1);
  tree_data->SetBranchAddress("flag_MW", &flag_MW);
  tree_data->SetBranchStatus("flag_calo_square",1);
  tree_data->SetBranchAddress("flag_calo_square", &flag_calo_square);
  tree_data->SetBranchStatus("flag_source_square",1);
  tree_data->SetBranchAddress("flag_source_square", &flag_source_square);
  tree_data->SetBranchStatus("flag_last_column",1);
  tree_data->SetBranchAddress("flag_last_column", &flag_last_column);
  tree_data->SetBranchStatus("flag_last_z",1);
  tree_data->SetBranchAddress("flag_last_z", &flag_last_z);
  tree_data->SetBranchStatus("energy",1);
  tree_data->SetBranchAddress("energy", &energyvis);
  tree_data->SetBranchStatus("calo_nohit_om_time",1);
  tree_data->SetBranchAddress("calo_nohit_om_time", &calo_nohit_om_time);


  for (int i = 0; i < tree_data->GetEntries(); i++) {
    tree_data->GetEntry(i);
    for (int j = 0; j < flag_e_event->size(); j++){
      if (flag_e_event->at(j) == 1 &&  calo_nohit_om_time->at(j) < 2 ) {
        data->Fill(om_number->at(j),energyvis->at(j)*simu_mean[om_number->at(j)]);
      }
    }
  }

  TH1D *spectre_data = new TH1D ("spectre_data", "spectre_data", 10, 0, 3) ;
  TH1D *spectre_simu = new TH1D ("spectre_simu", "spectre_simu", 10, 0, 3) ;
  TH1D *spectre_simu_ubc = new TH1D ("spectre_simu_ubc", "spectre_simu_ubc", 10, 0, 3) ;

  for (int om = 0; om < 520; om++) {

    TCanvas* c = new TCanvas;
    spectre_data = data->ProjectionY("Data", om+1, om+1);
    ScaleXaxis(spectre_data, ScaleX, om);
    spectre_simu = simu->ProjectionY("Simu", om+1, om+1);
    spectre_simu_ubc = simu_ubc->ProjectionY("Simu_ubc", om+1, om+1);

    for (int i = 0; i < 200; i++) {
      spectre_data->SetBinError(i,0);
    }

    TLegend *legend = new TLegend(0.3,0.8,0.45,0.9);
    legend->AddEntry(spectre_data, "data", "F");
    legend->AddEntry(spectre_simu, "simu", "F");
    legend->AddEntry(spectre_simu_ubc, "simu_ubc", "F");

    spectre_simu->Draw();
    spectre_simu->SetLineColor(kRed);
    spectre_simu->SetLineWidth(2);

    spectre_data->Draw("sames");
    spectre_data->SetLineWidth(2);
    spectre_data->SetLineColor(kBlack);
    spectre_data->Scale(1/1.78);

    spectre_simu_ubc->Draw("sames");
    spectre_simu_ubc->SetLineWidth(2);
    spectre_simu_ubc->SetLineColor(kGreen);

    legend->Draw("same");

    c->SaveAs(Form("Comparaison_correction/comparaison_OM_%d.png",om));
  }


}


void outfile() {
  int om_number;
  float FWHM;
  TFile tree_file2("Simu/eres_1059_bdf.root", "READ");
  TTree* tree2 = (TTree*)tree_file2.Get("Result_tree");
  tree2->SetBranchStatus("*",0);
  tree2->SetBranchStatus("om",1);
  tree2->SetBranchAddress("om", &om_number);
  tree2->SetBranchStatus("eres",1);
  tree2->SetBranchAddress("eres", &FWHM);
  std::ofstream outFile("FWHM_1059_bdf.txt");

  int compteur = 0;
  for (size_t i = 0; i < 520; i++) {
    tree2->GetEntry(compteur);
    if (om_number == i) {
      outFile << i << "\t"<< FWHM << endl;
      compteur++;
    }
    else outFile << i << "\t"<< 999 << endl;
  }


}

void time_part_biplot ()
{
  gStyle->SetOptFit(111);
  int om_number, position;
  double mean, sigma, mean_error, sigma_error;
  float eres;
  TFile *file = new TFile("Bi_fit/fitted_bi_time.root", "READ");
  TTree* tree = (TTree*)file->Get("Result_tree");
  tree->SetBranchStatus("*",0);
  tree->SetBranchStatus("om_number",1);
  tree->SetBranchAddress("om_number", &om_number);
  tree->SetBranchStatus("mean",1);
  tree->SetBranchAddress("mean", &mean);
  tree->SetBranchStatus("mean_error",1);
  tree->SetBranchAddress("mean_error", &mean_error);
  tree->SetBranchStatus("sigma",1);
  tree->SetBranchAddress("sigma", &sigma);
  tree->SetBranchStatus("sigma_error",1);
  tree->SetBranchAddress("sigma_error", &sigma_error);
  tree->SetBranchStatus("position",1);
  tree->SetBranchAddress("position", &position);

  double absolute_mean_time[9][520];
  memset (absolute_mean_time, 0, 9*520*sizeof(double));
  double absolute_mean_time_error[9][520];
  memset (absolute_mean_time_error, 0, 9*520*sizeof(double));
  double mean_time[9][520];
  memset (mean_time, 0, 9*520*sizeof(double));
  double mean_time_error[9][520];
  memset (mean_time_error, 0, 9*520*sizeof(double));

  for (int i = 0; i < tree->GetEntries(); i++) {
    tree->GetEntry(i);
    if (mean > 0) {
      mean_time[position][om_number] = 235*sigma/mean;

      mean_time_error[position][om_number] = mean_time[position][om_number]*(sigma_error/sigma + mean_error/mean);
      absolute_mean_time[position][om_number] = mean;
      absolute_mean_time_error[position][om_number] = mean_error;
    }
  }

  TFile *file2 = new TFile("Bi_fit/fitted_bi.root", "READ");
  TTree* tree2 = (TTree*)file2->Get("Result_tree");
  tree2->SetBranchStatus("*",0);
  tree2->SetBranchStatus("mean",1);
  tree2->SetBranchAddress("mean", &mean);
  tree2->SetBranchStatus("mean_error",1);
  tree2->SetBranchAddress("mean_error", &mean_error);
  tree2->SetBranchStatus("sigma",1);
  tree2->SetBranchAddress("sigma", &sigma);
  tree2->SetBranchStatus("sigma_error",1);
  tree2->SetBranchAddress("sigma_error", &sigma_error);
  tree2->SetBranchStatus("om_number",1);
  tree2->SetBranchAddress("om_number", &om_number);

  double eres_mean_tab[520];
  memset (eres_mean_tab, 0, 520*sizeof(double));
  double eres_mean_tab_error[520];
  memset (eres_mean_tab_error, 0, 520*sizeof(double));
  double eres_tab[520];
  memset (eres_tab, 0, 520*sizeof(double));
  double eres_tab_error[520];
  memset (eres_tab_error, 0, 520*sizeof(double));
  int compteur2 =0;
  for (int i = 0; i < tree2->GetEntries(); i++) {
    tree2->GetEntry(i);
    if (mean > 0) {
      eres_tab[om_number] = 235*sigma/mean;
      eres_tab_error[om_number] = eres_tab[om_number]*(sigma_error/sigma + mean_error/mean);
      eres_mean_tab[om_number] = mean;
      eres_mean_tab_error[om_number] = mean_error;
    }
  }

  TFile newfile("gain_evolution_974.root","RECREATE");
  int om;
  double p0, p1, p0_error, Chi2;
  TTree Result_tree("Result_tree","");
  Result_tree.Branch("om_number", &om);
  Result_tree.Branch("constant", &p0);
  Result_tree.Branch("constant_error", &p0_error);
  Result_tree.Branch("Chi2_ndf", &Chi2);


  TGraphErrors *biplot_delta_1 = new TGraphErrors;
  TGraphErrors *biplot_delta_2 = new TGraphErrors;
  TGraphErrors *biplot_delta_3 = new TGraphErrors;
  TGraphErrors *biplot_delta_4 = new TGraphErrors;
  gROOT->cd();

  for (int i = 0; i < 520; i++) {
    if ( i%13 !=0 && i%13 !=12 && eres_tab[i] > 0) {
      TGraphErrors *om_mean = new TGraphErrors;
      TGraphErrors *om_mean_normalized = new TGraphErrors;

      TCanvas *c = new TCanvas();
      for (int j = 0; j < 9; j++) {
        int point = om_mean->GetN();
        om_mean->SetPoint(point, j, absolute_mean_time[j][i]);
        om_mean->SetPointError(point, 0, absolute_mean_time_error[j][i]);
        om_mean_normalized->SetPoint(point, j, absolute_mean_time[j][i]/eres_mean_tab[i]);
        om_mean_normalized->SetPointError(point, 0, sqrt(pow(absolute_mean_time_error[j][i]/eres_mean_tab[i],2) + pow(eres_mean_tab_error[i]*absolute_mean_time[j][i]/pow(eres_mean_tab[i],2),2)));
      }

      om_mean->Draw();
      om_mean->GetXaxis()->SetTitle("Time (h)");
      om_mean->GetYaxis()->SetTitle("mean");
      om_mean->SetMarkerStyle(2);
      c->SaveAs(Form("Bi_fit/temps/individual_om_fit/fit_om_%d.png",i));
      delete om_mean;

      TCanvas *c2 = new TCanvas();
      om_mean_normalized->Draw();
      om_mean_normalized->GetXaxis()->SetTitle("Time (h)");
      om_mean_normalized->GetYaxis()->SetTitle("mean");
      om_mean_normalized->SetMarkerStyle(2);
      TF1 *f = new TF1("test","pol0");
      f->SetRange(0,9);
      f->Draw("same");
      om_mean_normalized->Fit(f, "RQ0");
      om_mean_normalized->GetYaxis()->SetRangeUser(0.9,1.1);

      p0 = f->GetParameter(0);
      p0_error = f->GetParError(0);
      Chi2 = f->GetChisquare()/f->GetNDF();
      om = i;
      Result_tree.Fill();
      c2->SaveAs(Form("Bi_fit/temps/individual_om_fit_normalized/fit_om_%d.root",i));
      c2->SaveAs(Form("Bi_fit/temps/individual_om_fit_normalized/fit_om_%d.png",i));
      delete om_mean_normalized;

    }
  }

  newfile.cd();
  Result_tree.Write();
  newfile.Close();



}

void time_full_biplot ()
{
  int om_number, position;
  double mean, sigma, mean_error, sigma_error;
  float eres;
  TFile *file = new TFile("Bi_fit/fitted_bi_time.root", "READ");
  TTree* tree = (TTree*)file->Get("Result_tree");
  tree->SetBranchStatus("*",0);
  tree->SetBranchStatus("om_number",1);
  tree->SetBranchAddress("om_number", &om_number);
  tree->SetBranchStatus("mean",1);
  tree->SetBranchAddress("mean", &mean);
  tree->SetBranchStatus("mean_error",1);
  tree->SetBranchAddress("mean_error", &mean_error);
  tree->SetBranchStatus("sigma",1);
  tree->SetBranchAddress("sigma", &sigma);
  tree->SetBranchStatus("sigma_error",1);
  tree->SetBranchAddress("sigma_error", &sigma_error);
  tree->SetBranchStatus("position",1);
  tree->SetBranchAddress("position", &position);

  double absolute_mean_time[4][520];
  memset (absolute_mean_time, 0, 4*520*sizeof(double));
  double absolute_mean_time_error[4][520];
  memset (absolute_mean_time_error, 0, 4*520*sizeof(double));
  double mean_time[4][520];
  memset (mean_time, 0, 4*520*sizeof(double));
  double mean_time_error[4][520];
  memset (mean_time_error, 0, 4*520*sizeof(double));

  for (int i = 0; i < tree->GetEntries(); i++) {
    tree->GetEntry(i);
    if (mean > 0) {
      mean_time[position][om_number] = 235*sigma/mean;

      mean_time_error[position][om_number] = mean_time[position][om_number]*(sigma_error/sigma + mean_error/mean);
      absolute_mean_time[position][om_number] = mean;
      absolute_mean_time_error[position][om_number] = mean_error;
    }
  }

  TFile *file2 = new TFile("Bi_fit/fitted_bi.root", "READ");
  TTree* tree2 = (TTree*)file2->Get("Result_tree");
  tree2->SetBranchStatus("*",0);
  tree2->SetBranchStatus("mean",1);
  tree2->SetBranchAddress("mean", &mean);
  tree2->SetBranchStatus("mean_error",1);
  tree2->SetBranchAddress("mean_error", &mean_error);
  tree2->SetBranchStatus("sigma",1);
  tree2->SetBranchAddress("sigma", &sigma);
  tree2->SetBranchStatus("sigma_error",1);
  tree2->SetBranchAddress("sigma_error", &sigma_error);
  tree2->SetBranchStatus("om_number",1);
  tree2->SetBranchAddress("om_number", &om_number);

  double eres_mean_tab[520];
  memset (eres_mean_tab, 0, 520*sizeof(double));
  double eres_mean_tab_error[520];
  memset (eres_mean_tab_error, 0, 520*sizeof(double));
  double eres_tab[520];
  memset (eres_tab, 0, 520*sizeof(double));
  double eres_tab_error[520];
  memset (eres_tab_error, 0, 520*sizeof(double));
  int compteur2 =0;
  for (int i = 0; i < tree2->GetEntries(); i++) {
    tree2->GetEntry(i);
    if (mean > 0) {
      eres_tab[om_number] = 235*sigma/mean;
      cout << "om : " << i << " mean time " <<  eres_tab[om_number]<< endl;
      eres_tab_error[om_number] = eres_tab[om_number]*(sigma_error/sigma + mean_error/mean);
      eres_mean_tab[om_number] = mean;
      eres_mean_tab_error[om_number] = mean_error;
    }
  }

  TGraphErrors *biplot_mean = new TGraphErrors;
  TGraphErrors *biplot_mean_1 = new TGraphErrors;
  TGraphErrors *biplot_mean_2 = new TGraphErrors;
  TGraphErrors *biplot_mean_3 = new TGraphErrors;
  TGraphErrors *biplot_mean_4 = new TGraphErrors;

  TGraphErrors *biplot_time_1 = new TGraphErrors;
  TGraphErrors *biplot_time_2 = new TGraphErrors;
  TGraphErrors *biplot_time_3 = new TGraphErrors;
  TGraphErrors *biplot_time_4 = new TGraphErrors;

  TGraphErrors *biplot_delta_1 = new TGraphErrors;
  TGraphErrors *biplot_delta_2 = new TGraphErrors;
  TGraphErrors *biplot_delta_3 = new TGraphErrors;
  TGraphErrors *biplot_delta_4 = new TGraphErrors;


  for (int i = 0; i < 520; i++) {
    if (i%13 != 11 && i%13 !=1 && i%13 !=0 && i%13 !=12 && eres_tab[i] > 0) {

      // int point1 = biplot_mean->GetN();
      // biplot_mean->SetPoint(point1, i, eres_mean_tab[i]);
      // biplot_mean->SetPointError(point1, 0, eres_mean_tab_error[i]);
      int point2 = biplot_mean_1->GetN();
      biplot_mean_1->SetPoint(point2, i, absolute_mean_time[0][i]/eres_mean_tab[i]);
      biplot_mean_1->SetPointError(point2, 0,sqrt(pow(absolute_mean_time_error[0][i]/eres_mean_tab[i],2) + pow(eres_mean_tab_error[i]*absolute_mean_time[0][i]/pow(eres_mean_tab[i],2),2)));
    }
  }


  biplot_time_1->GetXaxis()->SetTitle("eres full");
  biplot_time_1->GetYaxis()->SetTitle("eres debut/fin");
  biplot_time_1->SetMarkerStyle(8);
  biplot_time_1->SetMarkerColor(kRed);
  biplot_time_1->Draw("Ap");
  biplot_time_2->SetMarkerStyle(3);
  biplot_time_2->SetMarkerColor(kBlue);
  biplot_time_2->Draw("psame");
  biplot_time_3->SetMarkerStyle(4);
  biplot_time_3->SetMarkerColor(kBlack);
  biplot_time_3->Draw("psame");
  biplot_time_4->SetMarkerStyle(5);
  biplot_time_4->SetMarkerColor(8);
  biplot_time_4->Draw("psame");

  TLegend *legend = new TLegend(0.3,0.8,0.45,0.9);
  legend->AddEntry(biplot_time_1, "First quarter", "P");
  legend->AddEntry(biplot_time_2, "Second quarter", "P");
  legend->AddEntry(biplot_time_3, "Third quarter", "P");
  legend->AddEntry(biplot_time_4, "Fourth quarter", "P");

  legend->Draw("same");

  TCanvas *c = new TCanvas;
  c->cd();
  biplot_delta_1->GetXaxis()->SetTitle("OM");
  biplot_delta_1->GetYaxis()->SetTitle("Delta eres (full-debut/fin)");
  biplot_delta_1->SetMarkerStyle(2);
  biplot_delta_1->SetMarkerColor(kRed);
  biplot_delta_1->Draw("Ap");
  biplot_delta_2->SetMarkerStyle(3);
  biplot_delta_2->SetMarkerColor(kBlue);
  biplot_delta_2->Draw("psame");
  biplot_delta_3->SetMarkerStyle(4);
  biplot_delta_3->SetMarkerColor(kBlack);
  biplot_delta_3->Draw("psame");
  biplot_delta_4->SetMarkerStyle(5);
  biplot_delta_4->SetMarkerColor(8);
  biplot_delta_4->Draw("psame");

  TLegend *legend2 = new TLegend(0.3,0.8,0.45,0.9);
  legend2->AddEntry(biplot_delta_1, "First quarter", "P");
  legend2->AddEntry(biplot_delta_2, "Second quarter", "P");
  legend2->AddEntry(biplot_delta_3, "Third quarter", "P");
  legend2->AddEntry(biplot_delta_4, "Fourth quarter", "P");
  legend2->Draw("same");

  TCanvas *c2 = new TCanvas;
  c2->cd();
  biplot_mean_1->GetXaxis()->SetTitle("OM");
  biplot_mean_1->GetYaxis()->SetTitle("mean");
  biplot_mean_1->SetMarkerStyle(2);
  biplot_mean_1->SetMarkerColor(kRed);
  biplot_mean_1->Draw("Ap");
  biplot_mean_2->SetMarkerStyle(3);
  biplot_mean_2->SetMarkerColor(kBlue);
  biplot_mean_2->Draw("psame");
  biplot_mean_3->SetMarkerStyle(4);
  biplot_mean_3->SetMarkerColor(kBlack);
  biplot_mean_3->Draw("psame");
  biplot_mean_4->SetMarkerStyle(5);
  biplot_mean_4->SetMarkerColor(8);
  biplot_mean_4->Draw("psame");
  // biplot_mean->SetMarkerStyle(30);
  // biplot_mean->SetMarkerColor(6);
  // biplot_mean->Draw("psame");

  TLegend *legend3 = new TLegend(0.3,0.8,0.45,0.9);
  legend3->AddEntry(biplot_mean_1, "First quarter", "P");
  legend3->AddEntry(biplot_mean_2, "Second quarter", "P");
  legend3->AddEntry(biplot_mean_3, "Third quarter", "P");
  legend3->AddEntry(biplot_mean_4, "Fourth quarter", "P");
  legend3->AddEntry(biplot_mean, "Total", "P");
  legend3->Draw("same");



}

void time_biplot ()
{
  int om_number, position;
  double mean, sigma, mean_error, sigma_error;
  float eres;
  TFile *file = new TFile("Bi_fit/fitted_bi_time.root", "READ");
  TTree* tree = (TTree*)file->Get("Result_tree");
  tree->SetBranchStatus("*",0);
  tree->SetBranchStatus("om_number",1);
  tree->SetBranchAddress("om_number", &om_number);
  tree->SetBranchStatus("mean",1);
  tree->SetBranchAddress("mean", &mean);
  tree->SetBranchStatus("mean_error",1);
  tree->SetBranchAddress("mean_error", &mean_error);
  tree->SetBranchStatus("sigma",1);
  tree->SetBranchAddress("sigma", &sigma);
  tree->SetBranchStatus("sigma_error",1);
  tree->SetBranchAddress("sigma_error", &sigma_error);
  tree->SetBranchStatus("position",1);
  tree->SetBranchAddress("position", &position);

  double absolute_mean_time[4][520];
  memset (absolute_mean_time, 0, 4*520*sizeof(double));
  double absolute_mean_time_error[4][520];
  memset (absolute_mean_time_error, 0, 4*520*sizeof(double));
  double mean_time[4][520];
  memset (mean_time, 0, 4*520*sizeof(double));
  double mean_time_error[4][520];
  memset (mean_time_error, 0, 4*520*sizeof(double));

  for (int i = 0; i < tree->GetEntries(); i++) {
    tree->GetEntry(i);
    if (mean > 0) {
      mean_time[position][om_number] = 235*sigma/mean;

      mean_time_error[position][om_number] = mean_time[position][om_number]*(sigma_error/sigma + mean_error/mean);
      absolute_mean_time[position][om_number] = mean;
      absolute_mean_time_error[position][om_number] = mean_error;
    }
  }

  TFile *file2 = new TFile("Bi_fit/fitted_bi.root", "READ");
  TTree* tree2 = (TTree*)file2->Get("Result_tree");
  tree2->SetBranchStatus("*",0);
  tree2->SetBranchStatus("mean",1);
  tree2->SetBranchAddress("mean", &mean);
  tree2->SetBranchStatus("mean_error",1);
  tree2->SetBranchAddress("mean_error", &mean_error);
  tree2->SetBranchStatus("sigma",1);
  tree2->SetBranchAddress("sigma", &sigma);
  tree2->SetBranchStatus("sigma_error",1);
  tree2->SetBranchAddress("sigma_error", &sigma_error);
  tree2->SetBranchStatus("om_number",1);
  tree2->SetBranchAddress("om_number", &om_number);

  double eres_mean_tab[520];
  memset (eres_mean_tab, 0, 520*sizeof(double));
  double eres_mean_tab_error[520];
  memset (eres_mean_tab_error, 0, 520*sizeof(double));
  double eres_tab[520];
  memset (eres_tab, 0, 520*sizeof(double));
  double eres_tab_error[520];
  memset (eres_tab_error, 0, 520*sizeof(double));
  int compteur2 =0;
  for (int i = 0; i < tree2->GetEntries(); i++) {
    tree2->GetEntry(i);
    if (mean > 0) {
      eres_tab[om_number] = 235*sigma/mean;
      cout << "om : " << i << " mean time " <<  eres_tab[om_number]<< endl;
      eres_tab_error[om_number] = eres_tab[om_number]*(sigma_error/sigma + mean_error/mean);
      eres_mean_tab[om_number] = mean;
      eres_mean_tab_error[om_number] = mean_error;
    }
  }

  TGraphErrors *biplot_mean = new TGraphErrors;
  TGraphErrors *biplot_mean_1 = new TGraphErrors;
  TGraphErrors *biplot_mean_2 = new TGraphErrors;
  TGraphErrors *biplot_mean_3 = new TGraphErrors;
  TGraphErrors *biplot_mean_4 = new TGraphErrors;

  TGraphErrors *biplot_time_1 = new TGraphErrors;
  TGraphErrors *biplot_time_2 = new TGraphErrors;
  TGraphErrors *biplot_time_3 = new TGraphErrors;
  TGraphErrors *biplot_time_4 = new TGraphErrors;

  TGraphErrors *biplot_delta_1 = new TGraphErrors;
  TGraphErrors *biplot_delta_2 = new TGraphErrors;
  TGraphErrors *biplot_delta_3 = new TGraphErrors;
  TGraphErrors *biplot_delta_4 = new TGraphErrors;


  for (int i = 0; i < 520; i++) {
    if (i%13 != 11 && i%13 !=1 && i%13 !=0 && i%13 !=12 && eres_tab[i] > 0) {

      // int point1 = biplot_mean->GetN();
      // biplot_mean->SetPoint(point1, i, eres_mean_tab[i]);
      // biplot_mean->SetPointError(point1, 0, eres_mean_tab_error[i]);
      int point2 = biplot_mean_1->GetN();
      biplot_mean_1->SetPoint(point2, i, absolute_mean_time[0][i]/eres_mean_tab[i]);
      biplot_mean_1->SetPointError(point2, 0,sqrt(pow(absolute_mean_time_error[0][i]/eres_mean_tab[i],2) + pow(eres_mean_tab_error[i]*absolute_mean_time[0][i]/pow(eres_mean_tab[i],2),2)));
      int point3 = biplot_mean_2->GetN();
      biplot_mean_2->SetPoint(point3, i, absolute_mean_time[1][i]/eres_mean_tab[i]);
      biplot_mean_2->SetPointError(point3, 0,sqrt(pow(absolute_mean_time_error[1][i]/eres_mean_tab[i],2) + pow(eres_mean_tab_error[i]*absolute_mean_time[1][i]/pow(eres_mean_tab[i],2),2)));
      int point4 = biplot_mean_3->GetN();
      biplot_mean_3->SetPoint(point4, i, absolute_mean_time[2][i]/eres_mean_tab[i]);
      biplot_mean_3->SetPointError(point4, 0,sqrt(pow(absolute_mean_time_error[2][i]/eres_mean_tab[i],2) + pow(eres_mean_tab_error[i]*absolute_mean_time[2][i]/pow(eres_mean_tab[i],2),2)));
      int point = biplot_mean_4->GetN();
      biplot_mean_4->SetPoint(point, i, absolute_mean_time[3][i]/eres_mean_tab[i]);
      biplot_mean_4->SetPointError(point, 0,sqrt(pow(absolute_mean_time_error[3][i]/eres_mean_tab[i],2) + pow(eres_mean_tab_error[i]*absolute_mean_time[3][i]/pow(eres_mean_tab[i],2),2)));


      int pt1 = biplot_time_1->GetN();
      biplot_time_1->SetPoint(pt1, eres_tab[i], mean_time[0][i]);
      biplot_time_1->SetPointError(pt1, eres_tab_error[i], mean_time_error[0][i]);
      int pt2 = biplot_time_2->GetN();
      biplot_time_2->SetPoint(pt2, eres_tab[i], mean_time[1][i]);
      biplot_time_2->SetPointError(pt2, eres_tab_error[i], mean_time_error[1][i]);
      int pt3 = biplot_time_3->GetN();
      biplot_time_3->SetPoint(pt3, eres_tab[i], mean_time[2][i]);
      biplot_time_3->SetPointError(pt3, eres_tab_error[i], mean_time_error[2][i]);
      int pt4 = biplot_time_4->GetN();
      biplot_time_4->SetPoint(pt4, eres_tab[i], mean_time[3][i]);
      biplot_time_4->SetPointError(pt4, eres_tab_error[i], mean_time_error[3][i]);

      int PT1 = biplot_delta_1->GetN();
      biplot_delta_1->SetPoint(PT1, i, mean_time[0][i] - eres_tab[i]);
      biplot_delta_1->SetPointError(PT1, 0, sqrt(pow(mean_time_error[0][i],2) - pow(eres_tab_error[i],2)));
      int PT2 = biplot_delta_2->GetN();
      biplot_delta_2->SetPoint(PT2, i, mean_time[1][i] - eres_tab[i]);
      biplot_delta_2->SetPointError(PT2, 0, sqrt(pow(mean_time_error[1][i],2) - pow(eres_tab_error[i],2)));
      int PT3 = biplot_delta_3->GetN();
      biplot_delta_3->SetPoint(PT3, i, mean_time[2][i] - eres_tab[i]);
      biplot_delta_3->SetPointError(PT3, 0, sqrt(pow(mean_time_error[2][i],2) - pow(eres_tab_error[i],2)));
      int PT4 = biplot_delta_4->GetN();
      biplot_delta_4->SetPoint(PT4, i, mean_time[3][i] - eres_tab[i]);
      biplot_delta_4->SetPointError(PT4, 0, sqrt(pow(mean_time_error[3][i],2) - pow(eres_tab_error[i],2)));



    }
  }


  biplot_time_1->GetXaxis()->SetTitle("eres full");
  biplot_time_1->GetYaxis()->SetTitle("eres debut/fin");
  biplot_time_1->SetMarkerStyle(8);
  biplot_time_1->SetMarkerColor(kRed);
  biplot_time_1->Draw("Ap");
  biplot_time_2->SetMarkerStyle(3);
  biplot_time_2->SetMarkerColor(kBlue);
  biplot_time_2->Draw("psame");
  biplot_time_3->SetMarkerStyle(4);
  biplot_time_3->SetMarkerColor(kBlack);
  biplot_time_3->Draw("psame");
  biplot_time_4->SetMarkerStyle(5);
  biplot_time_4->SetMarkerColor(8);
  biplot_time_4->Draw("psame");

  TLegend *legend = new TLegend(0.3,0.8,0.45,0.9);
  legend->AddEntry(biplot_time_1, "First quarter", "P");
  legend->AddEntry(biplot_time_2, "Second quarter", "P");
  legend->AddEntry(biplot_time_3, "Third quarter", "P");
  legend->AddEntry(biplot_time_4, "Fourth quarter", "P");

  legend->Draw("same");

  TCanvas *c = new TCanvas;
  c->cd();
  biplot_delta_1->GetXaxis()->SetTitle("OM");
  biplot_delta_1->GetYaxis()->SetTitle("Delta eres (full-debut/fin)");
  biplot_delta_1->SetMarkerStyle(2);
  biplot_delta_1->SetMarkerColor(kRed);
  biplot_delta_1->Draw("Ap");
  biplot_delta_2->SetMarkerStyle(3);
  biplot_delta_2->SetMarkerColor(kBlue);
  biplot_delta_2->Draw("psame");
  biplot_delta_3->SetMarkerStyle(4);
  biplot_delta_3->SetMarkerColor(kBlack);
  biplot_delta_3->Draw("psame");
  biplot_delta_4->SetMarkerStyle(5);
  biplot_delta_4->SetMarkerColor(8);
  biplot_delta_4->Draw("psame");

  TLegend *legend2 = new TLegend(0.3,0.8,0.45,0.9);
  legend2->AddEntry(biplot_delta_1, "First quarter", "P");
  legend2->AddEntry(biplot_delta_2, "Second quarter", "P");
  legend2->AddEntry(biplot_delta_3, "Third quarter", "P");
  legend2->AddEntry(biplot_delta_4, "Fourth quarter", "P");
  legend2->Draw("same");

  TCanvas *c2 = new TCanvas;
  c2->cd();
  biplot_mean_1->GetXaxis()->SetTitle("OM");
  biplot_mean_1->GetYaxis()->SetTitle("mean");
  biplot_mean_1->SetMarkerStyle(2);
  biplot_mean_1->SetMarkerColor(kRed);
  biplot_mean_1->Draw("Ap");
  biplot_mean_2->SetMarkerStyle(3);
  biplot_mean_2->SetMarkerColor(kBlue);
  biplot_mean_2->Draw("psame");
  biplot_mean_3->SetMarkerStyle(4);
  biplot_mean_3->SetMarkerColor(kBlack);
  biplot_mean_3->Draw("psame");
  biplot_mean_4->SetMarkerStyle(5);
  biplot_mean_4->SetMarkerColor(8);
  biplot_mean_4->Draw("psame");
  // biplot_mean->SetMarkerStyle(30);
  // biplot_mean->SetMarkerColor(6);
  // biplot_mean->Draw("psame");

  TLegend *legend3 = new TLegend(0.3,0.8,0.45,0.9);
  legend3->AddEntry(biplot_mean_1, "First quarter", "P");
  legend3->AddEntry(biplot_mean_2, "Second quarter", "P");
  legend3->AddEntry(biplot_mean_3, "Third quarter", "P");
  legend3->AddEntry(biplot_mean_4, "Fourth quarter", "P");
  legend3->AddEntry(biplot_mean, "Total", "P");
  legend3->Draw("same");



}

void run_eres_biplot ()
{

  double mean, sigma, mean_error, sigma_error;
  int om_number;
  TFile *file4 = new TFile("Simu/Bi_fit/fitted_bi.root", "READ");
  TTree* tree4 = (TTree*)file4->Get("Result_tree");
  tree4->SetBranchStatus("*",0);
  tree4->SetBranchStatus("om_number",1);
  tree4->SetBranchAddress("om_number", &om_number);
  tree4->SetBranchStatus("mean",1);
  tree4->SetBranchAddress("mean", &mean);

  double means_tab[520];
  double sigmas_tab[520];
  memset (means_tab, 0, 520*sizeof(double));
  memset (sigmas_tab, 0, 520*sizeof(double));

  for (int i = 0; i < tree4->GetEntries(); i++) {
    tree4->GetEntry(i);
    if (mean > 0) {
      means_tab[om_number] = mean;
    }
  }

  int bad[34] = {2,5,84,118,120,121,122,123,131,132,133,134,192,292,316,330,340,348,366,368,370,458,470,471,491,495,508,510,511,512,513,514,515,516};
  int bad_ep[16] = {2,5,84,192,292,316,330,340,348,366,368,370,458,470,471,491};
  int bad_it[9] = {118,120,121,122,123,131,132,133,134};
  int bad_fr[9] = {495,508,510,511,512,513,514,515,516};
  float eres;
  TFile *file = new TFile("Bi_fit/fitted_bi_1048.root", "READ");
  TTree* tree = (TTree*)file->Get("Result_tree");
  tree->SetBranchStatus("*",0);
  tree->SetBranchStatus("om_number",1);
  tree->SetBranchAddress("om_number", &om_number);
  tree->SetBranchStatus("mean",1);
  tree->SetBranchAddress("mean", &mean);
  tree->SetBranchStatus("mean_error",1);
  tree->SetBranchAddress("mean_error", &mean_error);
  tree->SetBranchStatus("sigma",1);
  tree->SetBranchAddress("sigma", &sigma);
  tree->SetBranchStatus("sigma_error",1);
  tree->SetBranchAddress("sigma_error", &sigma_error);

  double absolute_mean_time[2][520];
  memset (absolute_mean_time, 0, 2*520*sizeof(double));
  double absolute_mean_time_error[2][520];
  memset (absolute_mean_time_error, 0, 2*520*sizeof(double));
  double mean_time[2][520];
  memset (mean_time, 0, 2*520*sizeof(double));
  double mean_time_error[2][520];
  memset (mean_time_error, 0, 2*520*sizeof(double));

  for (int i = 0; i < tree->GetEntries(); i++) {
    tree->GetEntry(i);
    if (mean > 0) {
      mean_time[0][om_number] = sqrt(means_tab[om_number])*235*sigma/mean;
      mean_time_error[0][om_number] = sqrt(means_tab[om_number])*mean_time[0][om_number]*(sigma_error/sigma + mean_error/mean);
      absolute_mean_time[0][om_number] = mean;
      absolute_mean_time_error[0][om_number] = mean_error;
    }
  }

  TFile *file3 = new TFile("Bi_fit/fitted_bi_1055.root", "READ");
  TTree* tree3 = (TTree*)file3->Get("Result_tree");
  tree3->SetBranchStatus("*",0);
  tree3->SetBranchStatus("om_number",1);
  tree3->SetBranchAddress("om_number", &om_number);
  tree3->SetBranchStatus("mean",1);
  tree3->SetBranchAddress("mean", &mean);
  tree3->SetBranchStatus("mean_error",1);
  tree3->SetBranchAddress("mean_error", &mean_error);
  tree3->SetBranchStatus("sigma",1);
  tree3->SetBranchAddress("sigma", &sigma);
  tree3->SetBranchStatus("sigma_error",1);
  tree3->SetBranchAddress("sigma_error", &sigma_error);

  for (int i = 0; i < tree->GetEntries(); i++) {
    tree3->GetEntry(i);
    if (mean > 0) {
      mean_time[1][om_number] = sqrt(means_tab[om_number])*235*sigma/mean;
      mean_time_error[1][om_number] = sqrt(means_tab[om_number])*mean_time[1][om_number]*(sigma_error/sigma + mean_error/mean);
      absolute_mean_time[1][om_number] = mean;
      absolute_mean_time_error[1][om_number] = mean_error;
    }
  }

  TFile *file2 = new TFile("Bi_fit/fitted_bi_974.root", "READ");
  TTree* tree2 = (TTree*)file2->Get("Result_tree");
  tree2->SetBranchStatus("*",0);
  tree2->SetBranchStatus("mean",1);
  tree2->SetBranchAddress("mean", &mean);
  tree2->SetBranchStatus("mean_error",1);
  tree2->SetBranchAddress("mean_error", &mean_error);
  tree2->SetBranchStatus("sigma",1);
  tree2->SetBranchAddress("sigma", &sigma);
  tree2->SetBranchStatus("sigma_error",1);
  tree2->SetBranchAddress("sigma_error", &sigma_error);
  tree2->SetBranchStatus("om_number",1);
  tree2->SetBranchAddress("om_number", &om_number);

  double eres_mean_tab[520];
  memset (eres_mean_tab, 0, 520*sizeof(double));
  double eres_mean_tab_error[520];
  memset (eres_mean_tab_error, 0, 520*sizeof(double));
  double eres_tab[520];
  memset (eres_tab, 0, 520*sizeof(double));
  double eres_tab_error[520];
  memset (eres_tab_error, 0, 520*sizeof(double));
  int compteur2 =0;
  for (int i = 0; i < tree2->GetEntries(); i++) {
    tree2->GetEntry(i);
    if (mean > 0) {
      eres_tab[om_number] = sqrt(means_tab[om_number])*235*sigma/mean;
      eres_tab_error[om_number] = sqrt(means_tab[om_number])*eres_tab[om_number]*(sigma_error/sigma + mean_error/mean);
      eres_mean_tab[om_number] = mean;
      eres_mean_tab_error[om_number] = mean_error;
    }
  }

  TGraphErrors *biplot_mean_1 = new TGraphErrors;
  TGraphErrors *biplot_mean_2 = new TGraphErrors;
  TGraphErrors *biplot_mean_1bad = new TGraphErrors;
  TGraphErrors *biplot_mean_2bad = new TGraphErrors;

  TGraphErrors *biplot_time_1 = new TGraphErrors;
  TGraphErrors *biplot_time_2 = new TGraphErrors;
  TGraphErrors *biplot_time_1bad = new TGraphErrors;
  TGraphErrors *biplot_time_2bad = new TGraphErrors;

  TGraphErrors *biplot_delta_1 = new TGraphErrors;
  TGraphErrors *biplot_delta_2 = new TGraphErrors;
  TGraphErrors *biplot_delta_1bad = new TGraphErrors;
  TGraphErrors *biplot_delta_2bad = new TGraphErrors;

  TGraphErrors *biplot_delta_new = new TGraphErrors;
  TGraphErrors *biplot_delta_newbad = new TGraphErrors;
  TGraphErrors *biplot_mean_new = new TGraphErrors;
  TGraphErrors *biplot_mean_newbad = new TGraphErrors;

  TH1D* dist_delta_eres_974_1048 = new TH1D("dist_delta_eres_974_1048", "dist_delta_eres_974_1048", 200,-10,10);
  TH1D* dist_delta_eres_974_1055 = new TH1D("dist_delta_eres_974_1055", "dist_delta_eres_974_1055", 200,-10,10);
  TH1D* dist_delta_eres_1048_1055 = new TH1D("dist_delta_eres_1048_1055", "dist_delta_eres_1048_1055", 200,-10,10);
  TH1D* dist_delta_eresbad_974_1048 = new TH1D("dist_delta_eres_974_1048", "dist_delta_eres_974_1048", 200,-10,10);
  TH1D* dist_delta_eresbad_974_1055 = new TH1D("dist_delta_eresbad_974_1055", "dist_delta_eresbad_974_1055", 200,-10,10);
  TH1D* dist_delta_eresbad_1048_1055 = new TH1D("dist_delta_eresbad_1048_1055", "dist_delta_eresbad_1048_1055", 200,-10,10);
  TH1D* dist_delta_eresbadit_974_1048 = new TH1D("dist_delta_eresit_974_1048", "dist_delta_eresit_974_1048", 200,-10,10);
  TH1D* dist_delta_eresbadit_974_1055 = new TH1D("dist_delta_eresbadit_974_1055", "dist_delta_eresbadit_974_1055", 200,-10,10);
  TH1D* dist_delta_eresbadit_1048_1055 = new TH1D("dist_delta_eresbadit_1048_1055", "dist_delta_eresbadit_1048_1055", 200,-10,10);
  TH1D* dist_delta_eresbadfr_974_1048 = new TH1D("dist_delta_eresfr_974_1048", "dist_delta_eresfr_974_1048", 200,-10,10);
  TH1D* dist_delta_eresbadfr_974_1055 = new TH1D("dist_delta_eresbadfr_974_1055", "dist_delta_eresbadfr_974_1055", 200,-10,10);
  TH1D* dist_delta_eresbadfr_1048_1055 = new TH1D("dist_delta_eresbadfr_1048_1055", "dist_delta_eresbadfr_1048_1055", 200,-10,10);

  int compteur = 0;
  int compteur_bad_it = 0;
  int compteur_bad_fr = 0;
  int compteur_bad_ep = 0;
  for (int i = 0; i < 520; i++) {
    if (i%13 !=0 && i%13 !=12 && eres_tab[i] > 0) {

      int point2 = biplot_mean_1->GetN();
      biplot_mean_1->SetPoint(point2, i, eres_mean_tab[i]/absolute_mean_time[0][i]);
      biplot_mean_1->SetPointError(point2, 0,sqrt(pow(eres_mean_tab_error[i]/absolute_mean_time[0][i],2) + pow(absolute_mean_time_error[0][i]*eres_mean_tab[i]/pow(absolute_mean_time[0][i],2),2)));
      int point3 = biplot_mean_2->GetN();
      biplot_mean_2->SetPoint(point3, i, eres_mean_tab[i]/absolute_mean_time[1][i]);
      biplot_mean_2->SetPointError(point3, 0,sqrt(pow(eres_mean_tab_error[i]/absolute_mean_time[1][i],2) + pow(absolute_mean_time_error[1][i]*eres_mean_tab[i]/pow(absolute_mean_time[1][i],2),2)));

      int pt1 = biplot_time_1->GetN();
      biplot_time_1->SetPoint(pt1, eres_tab[i], mean_time[0][i]);
      biplot_time_1->SetPointError(pt1, eres_tab_error[i], mean_time_error[0][i]);
      int pt2 = biplot_time_2->GetN();
      biplot_time_2->SetPoint(pt2, eres_tab[i], mean_time[1][i]);
      biplot_time_2->SetPointError(pt2, eres_tab_error[i], mean_time_error[1][i]);

      dist_delta_eres_974_1048->Fill(-mean_time[0][i] + eres_tab[i]);
      dist_delta_eres_974_1055->Fill(-mean_time[1][i] + eres_tab[i]);
      dist_delta_eres_1048_1055->Fill(-mean_time[1][i] + mean_time[0][i]);

      int PT1 = biplot_delta_1->GetN();
      biplot_delta_1->SetPoint(PT1, i, -mean_time[0][i] + eres_tab[i]);
      biplot_delta_1->SetPointError(PT1, 0, sqrt(pow(mean_time_error[0][i],2) + pow(eres_tab_error[i],2)));
      int PT2 = biplot_delta_2->GetN();
      biplot_delta_2->SetPoint(PT2, i, -mean_time[1][i] + eres_tab[i]);
      biplot_delta_2->SetPointError(PT2, 0, sqrt(pow(mean_time_error[1][i],2) + pow(eres_tab_error[i],2)));

      int PT1_new = biplot_delta_new->GetN();
      biplot_delta_new->SetPoint(PT1_new, i, -mean_time[0][i] + mean_time[1][i]);
      biplot_delta_new->SetPointError(PT1_new, 0, sqrt(pow(mean_time_error[0][i],2) + pow(mean_time_error[1][i],2)));

      int point2_new = biplot_mean_new->GetN();
      biplot_mean_new->SetPoint(point2_new, i, 100*(absolute_mean_time[0][i] - absolute_mean_time[1][i])/absolute_mean_time[0][i]);
      biplot_mean_new->SetPointError(point2_new, 0,100*sqrt(pow(absolute_mean_time_error[1][i]/absolute_mean_time[0][i],2) + pow(absolute_mean_time[1][i]*absolute_mean_time_error[0][i]/pow(absolute_mean_time[0][i],2),2)));

      if (i == bad_ep[compteur_bad_ep]) {
        compteur_bad_ep++;
        dist_delta_eresbad_974_1048->Fill(-mean_time[0][i] + eres_tab[i]);
        dist_delta_eresbad_974_1055->Fill(-mean_time[1][i] + eres_tab[i]);
        dist_delta_eresbad_1048_1055->Fill(-mean_time[1][i] + mean_time[0][i]);
      }
      if (i == bad_it[compteur_bad_it]) {
        compteur_bad_it++;
        dist_delta_eresbadit_974_1048->Fill(-mean_time[0][i] + eres_tab[i]);
        dist_delta_eresbadit_974_1055->Fill(-mean_time[1][i] + eres_tab[i]);
        dist_delta_eresbadit_1048_1055->Fill(-mean_time[1][i] + mean_time[0][i]);
      }
      if (i == bad_fr[compteur_bad_fr]) {
        compteur_bad_fr++;
        dist_delta_eresbadfr_974_1048->Fill(-mean_time[0][i] + eres_tab[i]);
        dist_delta_eresbadfr_974_1055->Fill(-mean_time[1][i] + eres_tab[i]);
        dist_delta_eresbadfr_1048_1055->Fill(-mean_time[1][i] + mean_time[0][i]);
      }

      if (i > bad[compteur]) compteur++;
      if (i == bad[compteur]) {
        compteur++;
        int point2bad = biplot_mean_1bad->GetN();
        biplot_mean_1bad->SetPoint(point2bad, i, eres_mean_tab[i]/absolute_mean_time[0][i]);
        biplot_mean_1bad->SetPointError(point2bad, 0,sqrt(pow(eres_mean_tab_error[i]/absolute_mean_time[0][i],2) + pow(absolute_mean_time_error[0][i]*eres_mean_tab[i]/pow(absolute_mean_time[0][i],2),2)));
        int point3bad = biplot_mean_2bad->GetN();
        biplot_mean_2bad->SetPoint(point3bad, i, eres_mean_tab[i]/absolute_mean_time[1][i]);
        biplot_mean_2bad->SetPointError(point3bad, 0,sqrt(pow(eres_mean_tab_error[i]/absolute_mean_time[1][i],2) + pow(absolute_mean_time_error[1][i]*eres_mean_tab[i]/pow(absolute_mean_time[1][i],2),2)));

        int pt1bad = biplot_time_1bad->GetN();
        biplot_time_1bad->SetPoint(pt1bad, eres_tab[i], mean_time[0][i]);
        biplot_time_1bad->SetPointError(pt1bad, eres_tab_error[i], mean_time_error[0][i]);
        int pt2bad = biplot_time_2bad->GetN();
        biplot_time_2bad->SetPoint(pt2bad, eres_tab[i], mean_time[1][i]);
        biplot_time_2bad->SetPointError(pt2bad, eres_tab_error[i], mean_time_error[1][i]);

        int PT1bad = biplot_delta_1bad->GetN();
        biplot_delta_1bad->SetPoint(PT1bad, i, -mean_time[0][i] + eres_tab[i]);
        biplot_delta_1bad->SetPointError(PT1bad, 0, sqrt(pow(mean_time_error[0][i],2) + pow(eres_tab_error[i],2)));
        int PT2bad = biplot_delta_2bad->GetN();
        biplot_delta_2bad->SetPoint(PT2bad, i, -mean_time[1][i] + eres_tab[i]);
        biplot_delta_2bad->SetPointError(PT2bad, 0, sqrt(pow(mean_time_error[1][i],2) + pow(eres_tab_error[i],2)));

        int PT1_newbad = biplot_delta_newbad->GetN();
        biplot_delta_newbad->SetPoint(PT1_newbad, i, -mean_time[0][i] + mean_time[1][i]);
        biplot_delta_newbad->SetPointError(PT1_newbad, 0, sqrt(pow(mean_time_error[0][i],2) + pow(mean_time_error[1][i],2)));

        int point2_newbad = biplot_mean_newbad->GetN();
        biplot_mean_newbad->SetPoint(point2_newbad, i, 100*(absolute_mean_time[0][i] - absolute_mean_time[1][i])/absolute_mean_time[0][i]);
        biplot_mean_newbad->SetPointError(point2_newbad, 0,100*sqrt(pow(absolute_mean_time_error[1][i]/absolute_mean_time[0][i],2) + pow(absolute_mean_time[1][i]*absolute_mean_time_error[0][i]/pow(absolute_mean_time[0][i],2),2)));
      }
    }
  }


  biplot_time_1->GetYaxis()->SetTitle("eres 1048 - 1055");
  biplot_time_1->GetXaxis()->SetTitle("eres 974");
  biplot_time_1->SetMarkerStyle(4);
  biplot_time_1->SetMarkerColor(17);
  biplot_time_1->Draw("Ap");
  biplot_time_2->SetMarkerStyle(25);
  biplot_time_2->SetMarkerColor(17);
  biplot_time_2->Draw("psame");
  biplot_time_1bad->SetMarkerStyle(8);
  biplot_time_1bad->SetMarkerColor(kRed);
  biplot_time_1bad->Draw("psame");
  biplot_time_2bad->SetMarkerStyle(21);
  biplot_time_2bad->SetMarkerColor(kBlue);
  biplot_time_2bad->Draw("psame");


  TLegend *legend = new TLegend(0.3,0.8,0.45,0.9);
  legend->AddEntry(biplot_time_1, "Eres 978/1048", "P");
  legend->AddEntry(biplot_time_1bad, "Eres 978/1048 bad OM", "P");
  legend->AddEntry(biplot_time_2, "Eres 978/1055", "P");
  legend->AddEntry(biplot_time_2bad, "Eres 978/1055 bad OM", "P");

  legend->Draw("same");

  TCanvas *c = new TCanvas;
  c->cd();
  biplot_delta_1->GetXaxis()->SetTitle("OM");
  biplot_delta_1->GetYaxis()->SetTitle("Delta eres ");
  biplot_delta_1->SetMarkerStyle(4);
  biplot_delta_1->SetMarkerColor(17);
  biplot_delta_1->Draw("Ap");
  biplot_delta_2->SetMarkerStyle(25);
  biplot_delta_2->SetMarkerColor(17);
  biplot_delta_2->Draw("psame");
  biplot_delta_1bad->SetMarkerStyle(8);
  biplot_delta_1bad->SetMarkerColor(kRed);
  biplot_delta_1bad->Draw("psame");
  biplot_delta_2bad->SetMarkerStyle(21);
  biplot_delta_2bad->SetMarkerColor(kBlue);
  biplot_delta_2bad->Draw("psame");

  TLegend *legend2 = new TLegend(0.3,0.8,0.45,0.9);
  legend2->AddEntry(biplot_delta_1, "Delta Eres 978/1048", "P");
  legend2->AddEntry(biplot_delta_1bad, "Delta Eres 978/1048 bad", "P");
  legend2->AddEntry(biplot_delta_2, "Delta Eres 978/1055", "P");
  legend2->AddEntry(biplot_delta_2bad, "Delta Eres 978/1055 bad", "P");
  legend2->Draw("same");

  TCanvas *c2 = new TCanvas;
  c2->cd();
  biplot_mean_1->GetXaxis()->SetTitle("OM");
  biplot_mean_1->GetYaxis()->SetTitle("Mean");
  biplot_mean_1->SetMarkerStyle(4);
  biplot_mean_1->SetMarkerColor(17);
  biplot_mean_1->Draw("Ap");
  biplot_mean_2->SetMarkerStyle(25);
  biplot_mean_2->SetMarkerColor(17);
  biplot_mean_2->Draw("psame");
  biplot_mean_1bad->SetMarkerStyle(8);
  biplot_mean_1bad->SetMarkerColor(kRed);
  biplot_mean_1bad->Draw("psame");
  biplot_mean_2bad->SetMarkerStyle(21);
  biplot_mean_2bad->SetMarkerColor(kBlue);
  biplot_mean_2bad->Draw("psame");

  TLegend *legend3 = new TLegend(0.3,0.8,0.45,0.9);
  legend3->AddEntry(biplot_mean_1, "Mean 1048", "P");
  legend3->AddEntry(biplot_mean_1bad, "Mean 1048 bad", "P");
  legend3->AddEntry(biplot_mean_2, "Mean 1055", "P");
  legend3->AddEntry(biplot_mean_2bad, "Mean 1055 bad", "P");
  legend3->Draw("same");

  TCanvas *c3 = new TCanvas;
  c3->cd();
  biplot_delta_new->GetXaxis()->SetTitle("OM");
  biplot_delta_new->GetYaxis()->SetTitle("Delta eres 1055-1048");
  biplot_delta_new->SetMarkerStyle(4);
  biplot_delta_new->SetMarkerColor(17);
  biplot_delta_new->Draw("Ap");
  biplot_delta_newbad->SetMarkerStyle(25);
  biplot_delta_newbad->SetMarkerColor(kRed);
  biplot_delta_newbad->Draw("psame");

  TCanvas *c4 = new TCanvas;
  c4->cd();
  biplot_mean_new->GetXaxis()->SetTitle("OM");
  biplot_mean_new->GetYaxis()->SetTitle("Delta mean 1055-1048");
  biplot_mean_new->SetMarkerStyle(4);
  biplot_mean_new->SetMarkerColor(17);
  biplot_mean_new->Draw("Ap");
  biplot_mean_newbad->SetMarkerStyle(25);
  biplot_mean_newbad->SetMarkerColor(kRed);
  biplot_mean_newbad->Draw("psame");


  TCanvas *c5 = new TCanvas;
  c5->cd();
  dist_delta_eres_974_1048->Draw();
  dist_delta_eres_974_1048->SetLineColor(27);
  dist_delta_eresbad_974_1048->Draw("histsames");
  dist_delta_eresbad_974_1048->SetLineColor(kRed);
  dist_delta_eresbad_974_1048->SetFillColor(kRed);
  dist_delta_eresbad_974_1048->SetFillStyle(3001);
  dist_delta_eresbadfr_974_1048->Draw("histsames");
  dist_delta_eresbadfr_974_1048->SetLineColor(kGreen);
  dist_delta_eresbadfr_974_1048->SetFillColor(kGreen);
  dist_delta_eresbadfr_974_1048->SetFillStyle(3001);
  dist_delta_eresbadit_974_1048->Draw("histsames");
  dist_delta_eresbadit_974_1048->SetLineColor(kBlue);
  dist_delta_eresbadit_974_1048->SetFillColor(kBlue);
  // dist_delta_eresbadit_974_1048->SetFillStyle(3001);

  TLegend *legend5 = new TLegend(0.3,0.8,0.45,0.9);
  legend5->AddEntry(dist_delta_eres_974_1048, "Delta FWHM @ 1 MeV", "F");
  legend5->AddEntry(dist_delta_eresbad_974_1048, "Delta FWHM @ 1 MeV bad", "F");
  legend5->AddEntry(dist_delta_eresbadfr_974_1048, "Delta FWHM @ 1 MeV bad fr cluster", "F");
  legend5->AddEntry(dist_delta_eresbadit_974_1048, "Delta FWHM @ 1 MeV bad it cluster", "F");
  legend5->Draw("same");

  TCanvas *c6 = new TCanvas;
  c6->cd();
  dist_delta_eres_974_1055->Draw();
  dist_delta_eres_974_1055->SetLineColor(27);
  dist_delta_eresbad_974_1055->Draw("histsames");
  dist_delta_eresbad_974_1055->SetLineColor(kRed);
  dist_delta_eresbad_974_1055->SetFillColor(kRed);
  dist_delta_eresbad_974_1055->SetFillStyle(3001);
  dist_delta_eresbadfr_974_1055->Draw("histsames");
  dist_delta_eresbadfr_974_1055->SetLineColor(kGreen);
  dist_delta_eresbadfr_974_1055->SetFillColor(kGreen);
  dist_delta_eresbadfr_974_1055->SetFillStyle(3001);
  dist_delta_eresbadit_974_1055->Draw("histsames");
  dist_delta_eresbadit_974_1055->SetLineColor(kBlue);
  dist_delta_eresbadit_974_1055->SetFillColor(kBlue);
  // dist_delta_eresbadit_974_1055->SetFillStyle(3001);

  TLegend *legend6 = new TLegend(0.3,0.8,0.45,0.9);
  legend6->AddEntry(dist_delta_eres_974_1055, "Delta FWHM @ 1 MeV", "F");
  legend6->AddEntry(dist_delta_eresbad_974_1055, "Delta FWHM @ 1 MeV bad", "F");
  legend6->AddEntry(dist_delta_eresbadfr_974_1055, "Delta FWHM @ 1 MeV bad fr cluster", "F");
  legend6->AddEntry(dist_delta_eresbadit_974_1055, "Delta FWHM @ 1 MeV bad it cluster", "F");
  legend6->Draw("same");


  TCanvas *c7 = new TCanvas;
  c7->cd();
  dist_delta_eres_1048_1055->Draw();
  dist_delta_eres_1048_1055->SetLineColor(27);
  dist_delta_eresbad_1048_1055->Draw("histsames");
  dist_delta_eresbad_1048_1055->SetLineColor(kRed);
  dist_delta_eresbad_1048_1055->SetFillColor(kRed);
  dist_delta_eresbad_1048_1055->SetFillStyle(3001);
  dist_delta_eresbadfr_1048_1055->Draw("histsames");
  dist_delta_eresbadfr_1048_1055->SetLineColor(kGreen);
  dist_delta_eresbadfr_1048_1055->SetFillColor(kGreen);
  dist_delta_eresbadfr_1048_1055->SetFillStyle(3001);
  dist_delta_eresbadit_1048_1055->Draw("histsames");
  dist_delta_eresbadit_1048_1055->SetLineColor(kBlue);
  dist_delta_eresbadit_1048_1055->SetFillColor(kBlue);
  // dist_delta_eresbadit_1048_1055->SetFillStyle(3001);

  TLegend *legend7 = new TLegend(0.3,0.8,0.45,0.9);
  legend7->AddEntry(dist_delta_eres_1048_1055, "Delta FWHM @ 1 MeV", "F");
  legend7->AddEntry(dist_delta_eresbad_1048_1055, "Delta FWHM @ 1 MeV bad", "F");
  legend7->AddEntry(dist_delta_eresbadfr_1048_1055, "Delta FWHM @ 1 MeV bad fr cluster", "F");
  legend7->AddEntry(dist_delta_eresbadit_1048_1055, "Delta FWHM @ 1 MeV bad it cluster", "F");
  legend7->Draw("same");
}

void eres_distrib ()
{

  int bad_ep[16] = {2,5,84,192,292,316,330,340,348,366,368,370,458,470,471,491};
  int bad_it[9] = {118,120,121,122,123,131,132,133,134};
  int bad_fr[9] = {495,508,510,511,512,513,514,515,516};

  float eres;
  int om;

  TFile *file = new TFile("Simu/eres_974.root", "READ");
  TTree* tree = (TTree*)file->Get("Result_tree");
  tree->SetBranchStatus("*",0);
  tree->SetBranchStatus("om",1);
  tree->SetBranchAddress("om", &om);
  tree->SetBranchStatus("eres",1);
  tree->SetBranchAddress("eres", &eres);

  TH1D* distrib_full = new TH1D("distrib_full", "distrib_full", 100,0,25);
  TH1D* distrib_it = new TH1D("distrib_it", "distrib_it", 100,0,25);
  TH1D* distrib_fr = new TH1D("distrib_fr", "distrib_fr", 100,0,25);
  TH1D* distrib_other = new TH1D("distrib_other", "distrib_other", 100,0,25);
  int compteur_it =0;
  int compteur_fr =0;
  int compteur_other =0;

  for (int omnum=0; omnum<tree->GetEntries(); ++omnum){ // MW
    tree->GetEntry(omnum);
    distrib_full->Fill(eres);

    if (bad_ep[compteur_other] == om) {
      distrib_other->Fill(eres);
      compteur_other++;
    }
    if (bad_it[compteur_it] == om) {
      cout << omnum << endl;

      distrib_it->Fill(eres);
      compteur_it++;
    }
    if (bad_fr[compteur_fr] == om) {
      distrib_fr->Fill(eres);
      compteur_fr++;
    }


  }

  distrib_full->Draw();
  distrib_fr->Draw("same");
  distrib_it->Draw("same");
  distrib_other->Draw("same");

  TLegend *legend = new TLegend(0.3,0.8,0.45,0.9);
  legend->AddEntry(distrib_full, "FWHM at 1 MeV", "F");
  legend->AddEntry(distrib_it, "FWHM at 1 MeV italian bad cluster", "F");
  legend->AddEntry(distrib_fr, "FWHM at 1 MeV french bad cluster", "F");
  legend->AddEntry(distrib_other, "FWHM at 1 MeV other bad OM", "F");
  legend->Draw("same");
}

void gain_variation() {

  gStyle->SetOptStat(0);
  gStyle->SetPalette(107);


  std::vector<int> *om_number = new std::vector<int>;
  std::vector<int> *energy = new std::vector<int>;
  std::vector<int> *flag_e_event = new std::vector<int>;
  std::vector<int> *flag_associated_nohit = new std::vector<int>;
  std::vector<double> *z_last_gg = new std::vector<double>;
  std::vector<double> *last_column = new std::vector<double>;
  std::vector<int> *calo_nohit_om_time = new std::vector<int>;
  std::vector<long> *calo_timestamp = new std::vector<long>;

  TFile tree_file("cut_1048.root", "READ");
  TTree* tree = (TTree*)tree_file.Get("Result_tree");
  tree->SetBranchStatus("*",0);
  tree->SetBranchStatus("om_number",1);
  tree->SetBranchAddress("om_number", &om_number);
  tree->SetBranchStatus("flag_e_event",1);
  tree->SetBranchAddress("flag_e_event", &flag_e_event);
  tree->SetBranchStatus("energy",1);
  tree->SetBranchAddress("energy", &energy);
  tree->SetBranchStatus("calo_timestamp",1);
  tree->SetBranchAddress("calo_timestamp", &calo_timestamp);
  gROOT->cd();

  TH3D *spectre = new TH3D("spectre","spectre", 520,0,520,100,0,10,100,0,4);
  tree->Project("spectre","om_number:calo_timestamp*6.25e-9/3600:energy","flag_e_event > 0 && flag_e_event < 9");


  TCanvas* cit = new TCanvas("cit","cit",4000,1000);
  cit->SetLogz();
  cit->Divide(13,20);

  TCanvas* cfr = new TCanvas("cfr","cfr",4000,1000);
  cfr->SetLogz();
  cfr->Divide(13,20);


  gROOT->cd();
  for (int i = 0; i < 520; i++) {
    spectre->GetXaxis()->SetRange(i+1,i+1);
    TH2D* gain_variation = (TH2D*) spectre->Project3D("x");
    gain_variation->GetXaxis()->SetTitle("Time (h)");
    gain_variation->GetYaxis()->SetTitle("Energy (MeV)");
    gain_variation->SetTitle(Form("energy variation OM %d",i));
    if (i < 260) {
      cit->cd(i);
      gain_variation->Draw("colz");
    }
    if (i > 259) {
      cfr->cd(i-260);
      gain_variation->Draw("colz");
    }
  }


}


void gamma_drawer ()
{
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
  std::vector<int> *flag_last_z = new std::vector<int>;
  std::vector<double> *z_last_gg = new std::vector<double>;
  std::vector<int> *flag = new std::vector<int>;
  std::vector<int> *calo_nohit_om_time = new std::vector<int>;

  TFile *file = new TFile("cut_bi_new.root", "READ");
  TTree* tree = (TTree*)file->Get("Result_tree");
  tree->SetBranchStatus("*",0);
  tree->SetBranchStatus("ncalo_tot",1);
  tree->SetBranchAddress("ncalo_tot", &ncalo_tot);
  tree->SetBranchStatus("om_num",1);
  tree->SetBranchAddress("om_num", &om_number);
  tree->SetBranchStatus("flag",1);
  tree->SetBranchAddress("flag", &flag);
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
  tree->SetBranchStatus("flag_last_z",1);
  tree->SetBranchAddress("flag_last_z", &flag_last_z);
  tree->SetBranchStatus("energyvis_ubc",1);
  tree->SetBranchAddress("energyvis_ubc", &energyvis_ubc);
  tree->SetBranchStatus("z_last_gg",1);
  tree->SetBranchAddress("z_last_gg", &z_last_gg);

  TH2D* gamma = new TH2D("gamma","gamma",520,0,520,1000,0,2);

  int e_compteur;
  int compteur =0;

  double rate_simu[712];
  memset (rate_simu, 0, 712*sizeof(double));
  for (int i = 0; i < tree->GetEntries()/10; i++) {
    tree->GetEntry(i);
    e_compteur = 0;
    for (int j = 0; j < flag_e_event->size(); j++)
      if (flag_e_event->at(j) == 1 && energyvis_ubc->at(j)>0.7 && ncalo_tot > 1) e_compteur++;

    for (int j = 0; j < flag_e_event->size(); j++) {
      if (flag_e_event->at(j) == 0 && e_compteur > 0 && flag_charge->at(j) == 1 && flag_MW->at(j) == 1 && flag_calo_square->at(j) == 0 && energyvis_ubc->at(j) > 0.22) {
        gamma->Fill(om_number->at(j), energyvis_ubc->at(j));
        rate_simu[om_number->at(j)]++;
        compteur++;
      }
    }
  }

    // gamma->Draw();
    // return;
  cout << "total simu event : " << compteur << endl;
  std::vector<long> *timestamp = new std::vector<long>;

  TFile *file2 = new TFile("../cut_974.root", "READ");
  TTree* tree2 = (TTree*)file2->Get("Result_tree");
  tree2->SetBranchStatus("*",0);
  tree2->SetBranchStatus("ncalo_tot",1);
  tree2->SetBranchAddress("ncalo_tot", &ncalo_tot);
  tree2->SetBranchStatus("om_number",1);
  tree2->SetBranchAddress("om_number", &om_number);
  tree2->SetBranchStatus("flag",1);
  tree2->SetBranchAddress("flag", &flag);
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
  tree2->SetBranchStatus("flag_last_z",1);
  tree2->SetBranchAddress("flag_last_z", &flag_last_z);
  tree2->SetBranchStatus("charge",1);
  tree2->SetBranchAddress("charge", &charge);
  tree2->SetBranchStatus("z_last_gg",1);
  tree2->SetBranchAddress("z_last_gg", &z_last_gg);
  tree2->SetBranchStatus("calo_nohit_om_time",1);
  tree2->SetBranchAddress("calo_nohit_om_time", &calo_nohit_om_time);
  tree2->SetBranchStatus("energy",1);
  tree2->SetBranchAddress("energy", &energyvis_ubc);
  tree2->SetBranchStatus("calo_timestamp",1);
  tree2->SetBranchAddress("calo_timestamp", &timestamp);

  int last_j;
  double rate[712];
  memset (rate, 0, 712*sizeof(double));
  int compteur2 =0;
  for (int i = 0; i < tree2->GetEntries()/10; i++) {
    tree2->GetEntry(i);
    e_compteur = 0;

    for (int j = 0; j < flag_e_event->size(); j++)
      if (flag_e_event->at(j) == 1 ) {
        e_compteur++;
        last_j = j;
      }
    if (e_compteur <1) continue;
    for (int j = 0; j < om_number->size(); j++) {

      float time = abs(timestamp->at(j)-timestamp->at(last_j));
            // cout << time << endl;/
      // if (time < 16) cout << time << " and "
      // if (ncalo_tot < 2 && flag_last_z->at(j) == 1) {
      // if ( energyvis_ubc->at(j)>0.7 && calo_nohit_om_time->at(j) < 2 && flag_charge->at(j) == 1 && flag_MW->at(j) == 1 && flag_associated_nohit->at(j)==1 && flag_last_column->at(j) == 1 && flag_calo_square->at(j) == 1 && flag_source_square->at(j) == 1 && flag_last_z->at(j) == 1)  {
      if (time < 16 && flag_e_event->at(j) == 0 && flag_charge->at(j) == 1 && flag_MW->at(j) == 1 && flag_calo_square->at(j) == 0 && energyvis_ubc->at(j) > 0.22) {
        rate[om_number->at(j)]++;
        compteur2++;
      }
      time = 0;
    }
  }
  //
  cout << "total data event : " << compteur2 << endl;
  // for (size_t i = 0; i < 520; i++) {
  //   std::cout << "rate simu = " << rate_simu[i] << " rate = " << rate[i] << '\n';
  // }



  // merge IT and FR canvas side by side using image magick (if installed)
  // gSystem->Exec("which convert > /dev/null && convert sndisplay-calorimeter-test-it.png sndisplay-calorimeter-test-fr.png +append sndisplay-calorimeter-test.png");
}


void Cut_rate ()
{

  int ncalo_tot;
  std::vector<int> *flag_e_event = new std::vector<int>;
  std::vector<int> *om_number = new std::vector<int>;
  std::vector<double> *energyvis_ubc = new std::vector<double>;
  std::vector<double> *energy = new std::vector<double>;
  std::vector<double> *charge = new std::vector<double>;
  std::vector<int> *flag_charge = new std::vector<int>;
  std::vector<int> *flag_associated_nohit = new std::vector<int>;
  std::vector<int> *flag_MW = new std::vector<int>;
  std::vector<int> *flag_calo_square = new std::vector<int>;
  std::vector<int> *flag_source_square = new std::vector<int>;
  std::vector<int> *flag_last_column = new std::vector<int>;
  std::vector<int> *flag_last_z = new std::vector<int>;
  std::vector<double> *z_last_gg = new std::vector<double>;
  std::vector<int> *flag = new std::vector<int>;
  std::vector<int> *calo_nohit_om_time = new std::vector<int>;

  TFile *file = new TFile("Simu/cut_bi_new.root", "READ");
  TTree* tree = (TTree*)file->Get("Result_tree");
  tree->SetBranchStatus("*",0);
  tree->SetBranchStatus("ncalo_tot",1);
  tree->SetBranchAddress("ncalo_tot", &ncalo_tot);
  tree->SetBranchStatus("om_num",1);
  tree->SetBranchAddress("om_num", &om_number);
  tree->SetBranchStatus("flag",1);
  tree->SetBranchAddress("flag", &flag);
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
  tree->SetBranchStatus("flag_last_z",1);
  tree->SetBranchAddress("flag_last_z", &flag_last_z);
  tree->SetBranchStatus("energyvis_ubc",1);
  tree->SetBranchAddress("energyvis_ubc", &energyvis_ubc);
  tree->SetBranchStatus("z_last_gg",1);
  tree->SetBranchAddress("z_last_gg", &z_last_gg);

  int compteur =0;
  double rate_simu_energie[712];
  double rate_simu[712];
  memset (rate_simu, 0, 712*sizeof(double));
  memset (rate_simu_energie, 0, 712*sizeof(double));
  for (int i = 0; i < tree->GetEntries()/10; i++) {
    tree->GetEntry(i);
    for (int j = 0; j < flag_e_event->size(); j++) {
      // if (energyvis_ubc->at(j) > 0.7 && ncalo_tot < 2 && flag_charge->at(j) == 1 && flag_MW->at(j) == 1 && flag_associated_nohit->at(j)==1 && flag_last_column->at(j) == 1 && flag_calo_square->at(j) == 1 && flag_source_square->at(j) == 1 && flag_last_z->at(j) == 1)  {
      if (ncalo_tot < 2 && flag_e_event->at(j) == 1) rate_simu[om_number->at(j)]++;
      if (ncalo_tot < 2 && flag_e_event->at(j) == 1 && energyvis_ubc->at(j)>0.7) {
        rate_simu_energie[om_number->at(j)]++;
        compteur++;
      }
    }
  }
  cout << "total simu event : " << compteur << endl;

  TFile *file2 = new TFile("cut_974.root", "READ");
  TTree* tree2 = (TTree*)file2->Get("Result_tree");
  tree2->SetBranchStatus("*",0);
  tree2->SetBranchStatus("ncalo_tot",1);
  tree2->SetBranchAddress("ncalo_tot", &ncalo_tot);
  tree2->SetBranchStatus("om_number",1);
  tree2->SetBranchAddress("om_number", &om_number);
  tree2->SetBranchStatus("flag",1);
  tree2->SetBranchAddress("flag", &flag);
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
  tree2->SetBranchStatus("flag_last_z",1);
  tree2->SetBranchAddress("flag_last_z", &flag_last_z);
  tree2->SetBranchStatus("charge",1);
  tree2->SetBranchAddress("charge", &charge);
  tree2->SetBranchStatus("z_last_gg",1);
  tree2->SetBranchAddress("z_last_gg", &z_last_gg);
  tree2->SetBranchStatus("calo_nohit_om_time",1);
  tree2->SetBranchAddress("calo_nohit_om_time", &calo_nohit_om_time);
  tree2->SetBranchStatus("energy",1);
  tree2->SetBranchAddress("energy", &energy);

  double rate[712];
  double rate_energie[712];
  memset (rate, 0, 712*sizeof(double));
  memset (rate_energie, 0, 712*sizeof(double));
  int compteur2 =0;
  for (int i = 0; i < tree2->GetEntries()/10; i++) {
    tree2->GetEntry(i);
    for (int j = 0; j < om_number->size(); j++) {
      // if (ncalo_tot < 2 && flag_last_z->at(j) == 1) {
      // if ( energyvis_ubc->at(j)>0.7 && calo_nohit_om_time->at(j) < 2 && flag_charge->at(j) == 1 && flag_MW->at(j) == 1 && flag_associated_nohit->at(j)==1 && flag_last_column->at(j) == 1 && flag_calo_square->at(j) == 1 && flag_source_square->at(j) == 1 && flag_last_z->at(j) == 1)  {
      if (calo_nohit_om_time ->at(j)< 2 && flag_e_event->at(j) == 1) rate[om_number->at(j)]++;
      if (calo_nohit_om_time->at(j) < 2 && flag_e_event->at(j) == 1 && energy->at(j)>0.7) {
        rate_energie[om_number->at(j)]++;
        compteur2++;
      }
    }
  }
  //
  cout << "total data event : " << compteur2 << endl;
  // for (size_t i = 0; i < 520; i++) {
  //   std::cout << "rate simu = " << rate_simu[i] << " rate = " << rate[i] << '\n';
  // }
  gROOT->cd();
  TH1D* Ratio = new TH1D("ratio","ratio",250,0,4);
  TH1D* Ratio_energy = new TH1D("ratio_energy","ratio_energy",250,0,4);

  for (int omnum=0; omnum<520; ++omnum) {// MW
    if (rate_simu[omnum] > 0 && rate[omnum] > 0 && omnum%13 != 1 && omnum%13 != 11){
      Ratio->Fill(rate_simu[omnum]/(rate[omnum]/1.78));
      Ratio_energy->Fill(rate_simu_energie[omnum]/(rate_energie[omnum]/1.78));
    }
  }

  Ratio->Draw();
  Ratio_energy->Draw("sames");
  Ratio_energy->SetLineColor(kRed);
}
