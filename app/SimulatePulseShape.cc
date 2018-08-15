#include <iostream>
#include <string>
//ROOT
#include <TCanvas.h>
#include <TFile.h>
#include <TTree.h>
//LOCAL
#include <PulseShape.hh>
#include <Configuration.hh>

bool PulseShape::_info    = false;
bool PulseShape::_debug   = false;
bool PulseShape::_warning = false;

int main ( int argc, char** argv )
{



  TH1F* h = new TH1F("pulse_time", "pulse_time", 2000, -10,10);
  TTree* pulse = new TTree("pulse", "Digitized waveforms");

  Configuration* config = new Configuration(true);
  config->GetCommandLineArgs(argc, argv);
  config->ParseConfigurationFile(config->config_file);

  const double n_threshold = config->n_threshold;
  const int n_experiments = config->n_experiments;
  const double DCR = config->DCR;
  const double Npe = config->Npe;
  const double scintillation_decay_constant = config->scintillation_decay_constant;
  const double single_photon_risetime_response = config->single_photon_risetime_response;
  const double single_photon_decaytime_response = config->single_photon_decaytime_response;
  const double high_pass_filter_RC = config->high_pass_filter_RC;

  std::cout << "number of experiments is: " << n_experiments << std::endl;
  std::cout << "Npe: " << Npe << std::endl;
  std::cout << "n_threshold: " << n_threshold << std::endl;
  std::cout << "scintillation decay constant: " << scintillation_decay_constant << " [ns]" << std::endl;
  std::cout << "scintillation_risetime: " << config->scintillation_risetime << " [ns]" << std::endl;
  std::cout << "single_photon_risetime_response: " << single_photon_risetime_response << " [ns]" << std::endl;
  std::cout << "single_photon_decaytime_response:" << single_photon_decaytime_response << " [ns]" << std::endl;
  std::cout << "high_pass_filter_RC: " << high_pass_filter_RC << std::endl;
  std::cout << "DCR: " << DCR << " [GHz]" << std::endl;


  PulseShape* ps;
  TGraph* total_pulse;
  TGraph* scintillation_pulse;
  TGraph* dark_noise;


  double step = 0.01;
  double x_low  = 0;
  double x_high = 50;

  const int npoints  = (x_high-x_low)/step;
  std::cout << "[INFO] number of points per pulse: " << npoints << std::endl;
  std::cout << "[INFO] sampling rate is: " << step  << " ns" << std::endl;
  double x[npoints];
  double y[npoints], y_sc[npoints], y_dc[npoints];
  int i_evt;
  pulse->Branch("i_evt", &i_evt, "i_evt/i");
  pulse->Branch("channel", y, Form("channel[%d]/D", npoints));
  pulse->Branch("channel_sc", y_sc, Form("channel_sc[%d]/D", npoints));
  pulse->Branch("channel_dc", y_dc, Form("channel_dc[%d]/D", npoints));
  pulse->Branch("time", x, Form("time[%d]/D", npoints));
  for ( int j = 0; j < n_experiments; j++ )
  {
    if ( j % 100 == 0 )std::cout << "experiment #" << j << std::endl;
    //reset variables and objects
    ps = new PulseShape("gauss");
    ps->SetNpe( Npe );
    ps->SetDCR( DCR );
    //ps->SetSinglePhotonResponse( single_pe_risetime );//units in ns
    ps->SetSinglePhotonRisetimeResponse( single_photon_risetime_response );
    ps->SetSinglePhotonDecaytimeResponse( single_photon_decaytime_response );
    ps->SetScintillationDecay( scintillation_decay_constant );//units in ns
    ps->SetHighPassFilterRC(high_pass_filter_RC);
    ps->NormalizeSinglePhotonResponse();
    //std::cout << ps->GetSinglePhotonResponseNormalization() << std::endl;
    ps->NormalizeSinglePhotonResponseHighPassFilter();
    //std::cout << ps->GetSinglePhotonResponseNormalization() << std::endl;
    for( int i = 0; i < npoints; i++ ) y[i] = x[i] = 0.0;
    double y_max = 0;

    //normalize pulse shape to area = 1;
    double normalization = 0;
    for( int i = 0; i < int(100 / 0.01); i++ ) normalization += ps->LGADShapedPulse( i * 0.01) * 0.01;

    for( int i = 0; i < npoints; i++ )
    {
      x[i] = x_low + double(i)*step;
      //if ( i % 1000 == 0 ) std::cout << "iteration #" << i << std::endl;
      //y[i]  = ps->Convolution(x[i], "Gauss", "RandomExp");
      y_sc[i]  = ps->LGADShapedPulse(x[i]) / normalization;
      //cout << i << " : " << y_sc[i] << "\n";
      // y_dc[i]  = ps->DarkNoise(x[i], x_low, x_high);
      y[i]     = y_sc[i];// + y_dc[i];
      if( y[i] > y_max ) y_max = y[i];
    }
    delete ps;//release memory of pulseshape object.

    double t_stamp = -999;
    for( int i = 0; i < npoints; i++ )
    {
      if ( y[i] > n_threshold )
      {
        t_stamp = (x[i]+x[i+1])/2.;
        break;
      }
    }

    h->Fill(t_stamp);
    i_evt = j;
    pulse->Fill();
  }

  ps = NULL;

  TFile* f = new TFile(config->output_file.c_str(), "recreate");
  total_pulse = new TGraph( npoints, x, y);
  scintillation_pulse = new TGraph( npoints, x, y_sc);
  dark_noise = new TGraph( npoints, x, y_dc);
  TCanvas *cv = new TCanvas("Cv","Cv", 800,800);
  cv->SetLeftMargin(0.13);
  cv->SetBottomMargin(0.12);
  cv->SetRightMargin(0.05);
  //h->GetXaxis()->SetRangeUser(-1e5, 0);
  total_pulse->Draw("AC*");
  total_pulse->GetXaxis()->SetTitle("time [ns]");
  total_pulse->GetYaxis()->SetTitle("Amplitude [normalized]");
  total_pulse->GetYaxis()->SetTitleOffset(1.6);
 
  cv->SaveAs("Convolution1.pdf");
  cv->SetLogy();
  cv->SaveAs("Convolution2.pdf");
  pulse->Write();
  total_pulse->Write("total_pulse");
  scintillation_pulse->Write("scintillation_pulse");
  dark_noise->Write("dark_noise");
  h->Write("h");
  f->Close();

  return 0;
};
