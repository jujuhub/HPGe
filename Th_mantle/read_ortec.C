// Reads HPGe spectra in txt (Ortec format) //
// - can combine files from job script      //
// Author: Vincent FISCHER
// A lot was copied from Steven GARDINER's script
// Modified by: Julie HE

#include <TFile.h>
#include <TSystem.h>
#include <TSystemDirectory.h>
#include <TSystemFile.h>
#include <TList.h>
#include <TMath.h>
#include <TH1.h>
#include <TString.h>
#include <TGraph.h>
#include <TCanvas.h>
#include <TPad.h>
#include <TLine.h>

#include <sys/stat.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>

using namespace std;
using namespace TMath;

void read_ortec(const char *dirname="/home/juhe/HPGe/Th_mantle/data/", const char *fext=".Spe") {

  // ---------------------------------------- //
  // ************ INITIALISATION ************ //
  // ---------------------------------------- //

  // create file to save plots in
//  TFile g("Th_mantle.root", "RECREATE");

  // look in data directory for all files
  string path = dirname;
  TSystemDirectory dir(dirname, dirname);
  TList *flist = dir.GetListOfFiles();

  // vectors/histograms to hold energy/channel+counts for ea file type
  TH1D *h_bg_chn;
  TH1D *h_bg_en;
  TH1D *h_mantle_chn;
  TH1D *h_mantle_en;
  TH1D *h_thox_chn;
  TH1D *h_thox_en;

  TH1D *h_mantle_chn_nobg;
  TH1D *h_mantle_en_nobg;
  TH1D *h_thox_chn_nobg;
  TH1D *h_thox_en_nobg;

  int fcounter = 0;
  Int_t line_number = 0;
  Int_t nb_chn = 0;
  Int_t nb_en_bins = 0;
  Int_t nb_energy_fit_param = 0;
  Double_t nb_bins_chn_plot, nb_bins_energy_plot;
  vector<Double_t> energy_fit_param;

  const char *pref1 = "background";
  const char *pref2 = "lantern_mantle";
  const char *pref3 = "thoxide";
  vector<Int_t> chn_bin_p1, nb_cnts_p1; // pref1
  vector<Int_t> chn_bin_p2, nb_cnts_p2; // pref2
  vector<Int_t> chn_bin_p3, nb_cnts_p3; // pref3
  vector<string> prefixes, times;
  
  Int_t live_time = 0, total_time = 0;
  
  // ---------------------------------------- //
  // ************** READ FILES ************** //
  // ---------------------------------------- //

  if (flist) {
    TSystemFile *file, *f0;
    TString fname, f0name;
    TIter next(flist);

    // might even want to define these outside if-statement
    string str;
    string sbuf;
    Double_t dbuf;

    // ONE-TIME PARAMS: set vector size, energy fit params
    size_t j = 0;
    f0 = dynamic_cast<TSystemFile *>(dir.GetListOfFiles()->At(0));
    f0name = f0->GetName();

    if (!f0name.EndsWith(fext)) {
      // if first file not .Spe file, look at next file
      while (!f0name.EndsWith(fext)) {
        ++j;
        f0 = dynamic_cast<TSystemFile *>(dir.GetListOfFiles()->At(j));
        f0name = f0->GetName();
      }
    }
    string f0name_str = f0name;
    cout << "initial Spe file: " << f0->GetName() << endl; // check

    size_t line_n = 0;
    string fullpath = path + f0name_str; 
    std::ifstream infile0(fullpath.c_str()); 
    while (getline(infile0, str)) {
      if (line_n == 11) { // $DATA
        stringstream ss(str);
        vector<string> nlines_vec;
        while (ss >> sbuf) {
          nlines_vec.push_back(sbuf);
        }
        nb_chn = atoi(nlines_vec[1].c_str());
//        break; // skip the rest of the file
      }
      if (nb_chn != 0) {
        if (line_n == nb_chn + 23) { // $MCA_CAL (2nd line)
          stringstream ss(str);
          while (ss >> dbuf) {
            energy_fit_param.push_back(dbuf);
            cout << dbuf << endl; // check
          }
          break;
        }
      }
      line_n++;
    }
    cout << "nb_chn: " << nb_chn << endl; // check

    fullpath.clear();
    chn_bin_p1.resize(nb_chn); nb_cnts_p1.resize(nb_chn); // pref1
    chn_bin_p2.resize(nb_chn); nb_cnts_p2.resize(nb_chn); // pref2
    chn_bin_p3.resize(nb_chn); nb_cnts_p3.resize(nb_chn); // pref3
    nb_en_bins = Floor(energy_fit_param[0] + nb_chn*energy_fit_param[1] + Power(nb_chn, 2)*energy_fit_param[2]);
    cout << "nb_en_bins: " << nb_en_bins << endl; // check


    // read in all files
    string fname_str;
    while ( (file = (TSystemFile*)next()) ) {
      fname = file->GetName();
      fname_str = fname;
      fullpath = path + fname_str;
 //     cout << fullpath << endl; // check
      std::ifstream infile(fullpath.c_str()); 

      if (!infile.good() || file->IsDirectory()) {
        cout << "***** File is not good. Skipping to next file... *****" << endl;
        cout << fname << endl;
        continue;
      }

      if (fname.EndsWith(fext)) {
        ++fcounter;
        while (getline(infile, str)) {
          // live and measured times
          if (line_number == 9) { // $MEAS_TIM
            stringstream ss(str);
            while (ss >> sbuf) {
              if (fname.BeginsWith(pref1)) {
                prefixes.push_back(pref1); //cout << pref1 << endl;
              } else if (fname.BeginsWith(pref2)) {
                prefixes.push_back(pref2); //cout << pref2 << endl;
              } else if (fname.BeginsWith(pref3)) {
                prefixes.push_back(pref3); //cout << pref2 << endl;
              }
              times.push_back(sbuf); //cout << sbuf << endl;
            }
          }

          // counts data
          if (line_number > 11 && line_number < nb_chn+13) {
            // background
            if (fname.BeginsWith(pref1)) {
                chn_bin_p1[line_number-12] = line_number-12;
                nb_cnts_p1[line_number-12] += atoi(str.c_str());
              }

            // lantern mantle
            if (fname.BeginsWith(pref2)) {
                chn_bin_p2[line_number-12] = line_number-12;
                nb_cnts_p2[line_number-12] += atoi(str.c_str());
            }

            // Th oxide
            if (fname.BeginsWith(pref3)) {
                chn_bin_p3[line_number-12] = line_number-12;
                nb_cnts_p3[line_number-12] += atoi(str.c_str());
            }
          }

          line_number++;
        }
      }

      // reset for next file
      line_number=0;
      fullpath.clear();
    }
  } else { cout << "!! WARNING !! No files found !!" << endl; }
  cout << "num of files processed: " << fcounter << endl;


  // ---------------------------------------- //
  // ************* HISTO DEFINE ************* //
  // ---------------------------------------- //
  h_bg_chn = new TH1D("bg_chn", "Channel spectrum: background", nb_chn, 0, nb_chn);
  h_bg_chn->Sumw2(); // not sure what this does
  h_mantle_chn = new TH1D("mantle_chn", "Channel spectrum: lantern mantle", nb_chn, 0, nb_chn);
  h_mantle_chn->Sumw2();
  h_thox_chn = new TH1D("thox_chn", "Channel spectrum: Th oxide", nb_chn, 0, nb_chn);
  h_thox_chn->Sumw2();

  h_bg_en = new TH1D("bg_en", "Energy [keV] spectrum: background", nb_en_bins, 0, nb_en_bins);
  h_bg_en->Sumw2();
  h_mantle_en = new TH1D("mantle_en", "Energy [keV] spectrum: lantern mantle", nb_en_bins, 0, nb_en_bins);
  h_mantle_en->Sumw2();
  h_thox_en = new TH1D("thox_en", "Energy [keV] spectrum: Th oxide", nb_en_bins, 0, nb_en_bins);
  h_thox_en->Sumw2();

  h_mantle_chn_nobg = new TH1D("mantle_chn_nobg", "Channel spectrum: lantern mantle (bg subtracted)", nb_chn, 0, nb_chn);
  h_mantle_en_nobg = new TH1D("mantle_en_nobg", "Energy [keV] spectrum: lantern mantle (bg subtracted", nb_en_bins, 0, nb_en_bins);
  h_thox_chn_nobg = new TH1D("thox_chn_nobg", "Channel spectrum: Th oxide (bg subtracted)", nb_chn, 0, nb_chn);
  h_thox_en_nobg = new TH1D("thox_en_nobg", "Energy [keV] spectrum: Th oxide (bg subtracted)", nb_en_bins, 0, nb_en_bins);


  // ---------------------------------------- //
  // ************** HISTO FILL ************** //
  // ---------------------------------------- //
  for(Int_t i = 0; i < nb_chn; ++i) {
    h_bg_chn->Fill(chn_bin_p1[i], nb_cnts_p1[i]);
    h_bg_en->Fill(energy_fit_param[0] + chn_bin_p1[i]*energy_fit_param[1] + chn_bin_p1[i]*chn_bin_p1[i]*energy_fit_param[2], nb_cnts_p1[i]);
    h_mantle_chn->Fill(chn_bin_p2[i], nb_cnts_p2[i]);
    h_mantle_en->Fill(energy_fit_param[0] + chn_bin_p2[i]*energy_fit_param[1] + chn_bin_p2[i]*chn_bin_p2[i]*energy_fit_param[2], nb_cnts_p2[i]);
    h_thox_chn->Fill(chn_bin_p3[i], nb_cnts_p3[i]);
    h_thox_en->Fill(energy_fit_param[0] + chn_bin_p3[i]*energy_fit_param[1] + chn_bin_p3[i]*chn_bin_p3[i]*energy_fit_param[2], nb_cnts_p3[i]);

    // background subtracted
    h_mantle_chn_nobg->Fill(chn_bin_p2[i], nb_cnts_p2[i]-nb_cnts_p1[i]);
    h_mantle_en_nobg->Fill(energy_fit_param[0] + chn_bin_p2[i]*energy_fit_param[1] + chn_bin_p2[i]*chn_bin_p2[i]*energy_fit_param[2], nb_cnts_p2[i]-nb_cnts_p1[i]);
    h_thox_chn_nobg->Fill(chn_bin_p3[i], nb_cnts_p3[i]-nb_cnts_p1[i]);
    h_thox_en_nobg->Fill(energy_fit_param[0] + chn_bin_p3[i]*energy_fit_param[1] + chn_bin_p3[i]*chn_bin_p3[i]*energy_fit_param[2], nb_cnts_p3[i]-nb_cnts_p1[i]);

//    Spectrum_energy->Fill(energy_fit_param[0] + chn_bin[i]*energy_fit_param[1] + chn_bin[i]*chn_bin[i]*energy_fit_param[2], height[i]);
  }


  // ---------------------------------------- //
  // ************** HISTO PLOT ************** //
  // ---------------------------------------- //

  /*
  TCanvas *c1 = new TCanvas("c1", "background_Chn", 0, 0, 1200, 1000);
  c1->cd();
  h_bg_chn->Draw("HIST");
  c1->Update();

  TCanvas *c2 = new TCanvas("c2", "mantle_Chn", 0, 0, 1200, 1000);
  c2->cd();
  h_mantle_chn->Draw("HIST");
  c2->Update();

  TCanvas *c3 = new TCanvas("c3", "thoxide_Chn", 0, 0, 1200, 1000);
  c3->cd();
  h_thox_chn->Draw("HIST");
  c3->Update();
  */

  TCanvas *c1 = new TCanvas("c1", "Th mantle analysis", 0, 0, 1200, 1000);
  c1->Divide(2,3);
  c1->cd(1);
  h_bg_chn->Draw("HIST"); // background
  c1->cd(2);
  h_bg_en->Draw("HIST");

  c1->cd(3);
  h_mantle_chn->Draw("HIST"); // lantern mantle
  c1->cd(4);
  h_mantle_en->Draw("HIST");

  c1->cd(5);
  h_thox_chn->Draw("HIST"); // Th oxide
  c1->cd(6);
  h_thox_en->Draw("HIST");
  c1->Update();

  // background subtracted plots
  TLine *l_238 = new TLine(238.6, 0., 238.6, 18000.); // Pb-212
  l_238->SetLineColor(kRed);
  l_238->SetLineStyle(2);
  TLine *l_338 = new TLine(338.3, 0., 338.3, 18000.); 
  l_338->SetLineColor(kRed);
  l_338->SetLineStyle(2);
  TLine *l_583 = new TLine(583.2, 0., 583.2, 18000.);
  l_583->SetLineColor(kRed);
  l_583->SetLineStyle(2);
  TLine *l_911 = new TLine(911.2, 0., 911.2, 18000.);
  l_911->SetLineColor(kRed);
  l_911->SetLineStyle(2);
  TLine *l_2614 = new TLine(2614.5, 0., 2614.5, 18000.); // Tl-208
  l_2614->SetLineColor(kRed);
  l_2614->SetLineStyle(2);

  TCanvas *c2 = new TCanvas("c2", "Th mantle analysis (background subtracted)", 0, 0, 1200, 1000);
  c2->Divide(2,2);
  c2->cd(1);
  h_mantle_chn_nobg->Draw("HIST"); // lantern mantle
  c2->cd(2);
  h_mantle_en_nobg->Draw("HIST");
  l_238->Draw("same");
  l_338->Draw("same");
  l_583->Draw("same");
  l_911->Draw("same");
  l_2614->Draw("same");

  c2->cd(3);
  h_thox_chn_nobg->Draw("HIST"); // Th oxide
  c2->cd(4);
  h_thox_en_nobg->Draw("HIST");
  l_238->Draw("same");
  l_338->Draw("same");
  l_583->Draw("same");
  l_911->Draw("same");
  l_2614->Draw("same");
  c2->Update();

}
