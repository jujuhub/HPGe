//////// Reads HPGE spectra in TKA (Canberra ASCII format) or txt (Ortec format) ////////
//////// Author: Vincent FISCHER ////////
//////// A lot was copied from Steven Gardiner's script ////////
//////// Modified by: Julie He ////////

#include <TFile.h>
#include <TTree.h>
#include <TH1.h>
#include <TMath.h>
#include <TCanvas.h>
#include <TLeaf.h>
#include <Rtypes.h>
#include <TROOT.h>
#include <TRandom3.h>
#include <TH2.h>
#include <TPad.h>
#include <TVector3.h>
#include <TString.h>
#include <TPRegexp.h>
#include <TGraph.h>

#include <sys/stat.h>
#include <iostream>
#include <fstream>
#include <vector>
#include <sstream>
#include <string>
#include <vector>

using namespace std;
using namespace TMath;

void read_spectra(const char* filename, TString HPGE_type="ortec") {
  
  // -----------------------------------------------------------------//
  // ************************ INITIALISATION ************************ //
  // -----------------------------------------------------------------//
  std::ifstream infile(filename);
  
  if (!infile.good()) {
    cout << "***** File is not good. Aborting... *****" << endl; return;
  }
  if (HPGE_type != "ortec" && HPGE_type != "canberra99" && HPGE_type != "canberra50" && HPGE_type != "canberra25" && HPGE_type != "canberra08") {
    cout << "***** HPGE type (brand) not recognized, only ortec, canberra99, canberra50, canberra25, canberra08 are accepted. Aborting... *****" << endl; return;
  }
  
  Int_t line_number = 0;
  Int_t number_channels = 0;
  Int_t energy_fit_param_nb = 0;
  Double_t Nb_bins_channel_plot, Nb_bins_energy_plot;
  vector<Double_t> energy_fit_param, canberra_energy_calib;
  energy_fit_param.clear(); canberra_energy_calib.clear();
  string str;
  string spectrum_name;
  vector<Int_t> channel_bin, height;
  vector<Double_t> energy_bin;
  channel_bin.clear();height.clear();energy_bin.clear();
  TH1D *Spectrum_channel, *Spectrum_energy;
  
  Int_t live_time = 0, total_time = 0;
  
  if (HPGE_type == "canberra99") { canberra_energy_calib.push_back(-0.360);canberra_energy_calib.push_back(0.501);canberra_energy_calib.push_back(2.2460e-8);}
  if (HPGE_type == "canberra50") { canberra_energy_calib.push_back(-0.250);canberra_energy_calib.push_back(0.533);canberra_energy_calib.push_back(3e-9);}
  if (HPGE_type == "canberra25") { canberra_energy_calib.push_back(-0.626);canberra_energy_calib.push_back(0.500);canberra_energy_calib.push_back(7.08e-8);}
  if (HPGE_type == "canberra08") { canberra_energy_calib.push_back(-0.343);canberra_energy_calib.push_back(0.500);canberra_energy_calib.push_back(3.5737e-8);}
  
  // -----------------------------------------------------------------//
  // ************************ READ THE FILES ************************ //
  // -----------------------------------------------------------------//
  if (HPGE_type == "ortec"){
    while(getline(infile, str)) {
      
      if (line_number == 1) {
        spectrum_name = str;
      }
      
      if (line_number == 9) { // $MEAS_TIM
        string buf; // Have a buffer string
        stringstream ss(str);
        vector<string> times;
        while (ss >> buf) {
          times.push_back(buf);
        } 
        live_time = atoi(times[0].c_str()); total_time = atoi(times[1].c_str());
      }
      
      if (line_number == 11) { // $DATA
        string sbuf; // Have a buffer string
        stringstream ss(str);
        vector<string> lines_vec;
        while (ss >> sbuf) {
          lines_vec.push_back(sbuf);
        }
        number_channels = atoi(lines_vec[1].c_str());
//        cout << number_channels << endl;
      }
      
      if (line_number > 11 && line_number < number_channels + 13) {
//        cout << str << endl;
        channel_bin.push_back(line_number-12); // start channel_bin at 0
        height.push_back(atoi(str.c_str()));
      }
      
      if (line_number == number_channels + 22) { // includes first 11 rows; $MCA_CAL
        energy_fit_param_nb = atoi(str.c_str());
      }
      
      if (line_number == number_channels + 23) { // $MCA_CAL (2nd line)
        Double_t dbuf; // Have a buffer double
        stringstream ss(str);
        while (ss >> dbuf) {
          energy_fit_param.push_back(dbuf);
        }
      }      
      
      line_number++;
    }
  } else if (HPGE_type == "canberra99" || HPGE_type == "canberra50" || HPGE_type == "canberra25" || HPGE_type == "canberra08") {
    while(getline(infile, str)) {
      if (line_number == 0) {
        live_time = atoi(str.c_str());
      }
      if (line_number == 1) {
        total_time = atoi(str.c_str());
      }
      
      if (line_number > 1 && str != "") {
        channel_bin.push_back(line_number-2); height.push_back(atoi(str.c_str()));
      }
      
      line_number++;
    }
  } 
  
  // -----------------------------------------------------------------//
  // ************************ HISTO DECLARATION ********************* //
  // -----------------------------------------------------------------//
  if (HPGE_type == "ortec") {
    Nb_bins_channel_plot = 20000, Nb_bins_energy_plot = Floor(energy_fit_param[0] + Nb_bins_channel_plot*energy_fit_param[1] + Power(Nb_bins_channel_plot, 2)*energy_fit_param[2]);

    Spectrum_channel = new TH1D("Spectrum_channel", "Spectrum in channels; Channel number; Entries", Nb_bins_channel_plot, 0, Nb_bins_channel_plot);
    Spectrum_channel->Sumw2();
    
    Spectrum_energy = new TH1D("Spectrum_energy", "Spectrum in energy; Energy [keV]; Entries", 20000, 0, Nb_bins_energy_plot);
    Spectrum_energy->Sumw2();

  } else if (HPGE_type == "canberra99" || HPGE_type == "canberra50" || HPGE_type == "canberra25" || HPGE_type == "canberra08") {
    Nb_bins_channel_plot = 5000, Nb_bins_energy_plot = canberra_energy_calib[0] + (Nb_bins_channel_plot*canberra_energy_calib[1]) + (canberra_energy_calib[2]*Power(Nb_bins_channel_plot, 2));

    Spectrum_channel = new TH1D("Spectrum_channel", "Spectrum in channels; Channel number; Entries", Nb_bins_channel_plot, 0, Nb_bins_channel_plot);
    Spectrum_channel->Sumw2();
    Spectrum_energy = new TH1D("Spectrum_energy", "Spectrum in energy; Energy [keV]; Entries", 5000, 0, Nb_bins_energy_plot);
    Spectrum_energy->Sumw2();
  }
  
  // -------------------------------------------------------------//
  // ************************ HISTO FILLING ********************* //
  // -------------------------------------------------------------//
  if (HPGE_type == "ortec") {
    for(Int_t i = 0; i < number_channels; ++i) {
      Spectrum_channel->Fill(channel_bin[i],height[i]);
      Spectrum_energy->Fill(energy_fit_param[0] + channel_bin[i]*energy_fit_param[1] + channel_bin[i]*channel_bin[i]*energy_fit_param[2], height[i]);
    }
  }
  if (HPGE_type == "canberra99" || HPGE_type == "canberra50" || HPGE_type == "canberra25" || HPGE_type == "canberra08") {
    for(UInt_t i = 0; i < channel_bin.size(); ++i) {
      Spectrum_channel->Fill(channel_bin[i], height[i]);
      Spectrum_energy->Fill(canberra_energy_calib[0] + (channel_bin[i]*canberra_energy_calib[1]) + (canberra_energy_calib[2]*Power(channel_bin[i], 2)), height[i]);

      // Energy-channel of 50 % HPGe: energy (keV)=-0.210+(0.532*Channel number)+(0.0000000203*channel number^2)
      // Energy-channel of 99 % HPGe: energy (keV)=-0.365+(0.501*channel number)-(0.0000000414*channel number^2)
    }
  }
  
  // --------------------------------------------------------------//
  // ************************ HISTO PLOTTING ********************* //
  // --------------------------------------------------------------//
  if (HPGE_type == "ortec") {
    TCanvas* c1 = new TCanvas("c1", "ortec_Chn", 0, 0, 1200, 1000);
    c1->cd();
    Spectrum_channel->Draw("HIST");
    c1->Update();
    TCanvas* c2 = new TCanvas("c2", "ortec_Spe", 0, 0, 1200, 1000);
    c2->cd();
    Spectrum_energy->Draw("HIST");
    c2->Update();
  }
  
  if (HPGE_type == "canberra99" || HPGE_type == "canberra50" || HPGE_type == "canberra25" || HPGE_type == "canberra08") {
    TCanvas* c1 = new TCanvas("c1", "canberra_Chn", 0,0, 1200, 1000);
    Spectrum_channel->DrawCopy("HIST");
    TCanvas* c2 = new TCanvas("c2", "canberra_Spe", 0,0, 1200, 1000);
    Spectrum_energy->DrawCopy("HIST");
  }

}
