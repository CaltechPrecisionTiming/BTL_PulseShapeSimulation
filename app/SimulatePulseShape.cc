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

  const int NFilter = config->NFilter;
  const int n_experiments = config->n_experiments;
  const double ShapingTime = config->ShapingTime;
  const double SNR = config->SNR;
  const int randomSeed = config->randomSeed;

  std::cout << "number of experiments is: " << n_experiments << std::endl;
  std::cout << "NFilter: " << NFilter << std::endl;
  std::cout << "ShapingTime: " << ShapingTime << " ns" << std::endl;
  std::cout << "Signal-to-Noise Ratio: " << SNR << std::endl;


  PulseShape* ps;
  TGraph* total_pulse;
  TGraph* scintillation_pulse;
  TGraph* dark_noise;
  TGraph* white_noise;
  TGraph* white_noise_pulse;

  double step = 0.01;
  double x_low  = 0;
  double x_high = 50;

  const int npoints  = (x_high-x_low)/step;
  std::cout << "[INFO] number of points per pulse: " << npoints << std::endl;
  std::cout << "[INFO] sampling rate is: " << step  << " ns" << std::endl;
  float x[npoints];
  float y[npoints], y_signal[npoints], y_dc[npoints], y_wn[npoints], y_wnps[npoints];
  int i_evt;
  pulse->Branch("i_evt", &i_evt, "i_evt/i");
  pulse->Branch("channel", y, Form("channel[%d]/F", npoints));
  pulse->Branch("shapednoise", y_wnps, Form("shapednoise[%d]/F", npoints));
  pulse->Branch("noise", y_wn, Form("noise[%d]/F", npoints));
  pulse->Branch("time", x, Form("time[%d]/F", npoints));
  for ( int j = 0; j < n_experiments; j++ )
  {
    if ( j % 1 == 0 )std::cout << "experiment #" << j << std::endl;
    //reset variables and objects
    for( int i = 0; i < npoints; i++ ) y[i] = x[i] = 0.0;
    double y_max = 0;

    //create pulse shape object
    ps = new PulseShape( ShapingTime, NFilter , SNR, randomSeed);


    //normalize pulse shape to have pulse height at 1.0
    double normalization = 0;
    for( int i = 0; i < int(100 / 0.01); i++ ) {
      double tmp = ps->LGADShapedPulse( i * 0.01);
      if ( tmp > normalization ) normalization = tmp;
    }


    //populate the pulse shape
    for( int i = 0; i < npoints; i++ )
    {
      x[i] = x_low + double(i)*step;
      //if ( i % 1000 == 0 ) std::cout << "iteration #" << i << std::endl;
      //y[i]  = ps->Convolution(x[i], "Gauss", "RandomExp");
      y_signal[i]   = ps->LGADShapedPulse(x[i]) / normalization;
      y_wn[i]   = ps->WhiteNoise(0, 1./30.);
      y_wnps[i] = ps->WhiteNoiseShapedPulse(x[i], 0, 1./30.);
      //y_signal[i]  = ps->NormalizedImpulseResponse(x[i]);

      //cout << i << " : " << y_signal[i] << "\n";
      // y_dc[i]  = ps->DarkNoise(x[i], x_low, x_high);
      y[i]     = y_signal[i] + y_wnps[i];// + y_dc[i];
      if( y[i] > y_max ) y_max = y[i];
    }
    delete ps;//release memory of pulseshape object.


    double t_stamp = -999;
    for( int i = 0; i < npoints; i++ )
    {
      if ( y[i] > 0.15 )
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
  scintillation_pulse = new TGraph( npoints, x, y_signal);
  dark_noise = new TGraph( npoints, x, y_dc);
  white_noise = new TGraph( npoints, x, y_wn);
  white_noise_pulse = new TGraph( npoints, x, y_wnps);
  /*
  plot total pulse
  */
  TCanvas *cv = new TCanvas("Cv","Cv", 800,800);
  cv->SetLeftMargin(0.13);
  cv->SetBottomMargin(0.12);
  cv->SetRightMargin(0.05);
  //h->GetXaxis()->SetRangeUser(-1e5, 0);
  total_pulse->SetTitle("");
  total_pulse->Draw("AC*");
  total_pulse->GetXaxis()->SetTitle("time [ns]");
  total_pulse->GetYaxis()->SetTitle("Amplitude [normalized]");
  total_pulse->GetYaxis()->SetTitleOffset(1.6);
  cv->SaveAs("Convolution1.pdf");
  cv->SetLogy();
  cv->SaveAs("Convolution2.pdf");

  /*
  plot white noise
  */
  white_noise->SetTitle("");
  white_noise->SetMarkerStyle(20);
  white_noise->SetMarkerSize(0.3);
  white_noise->SetMarkerColor(kBlue);
  white_noise->SetLineColor(kBlue);
  white_noise->Draw("APL");
  white_noise->GetXaxis()->SetTitle("time [ns]");
  white_noise->GetYaxis()->SetTitle("Amplitude [normalized]");
  white_noise->GetYaxis()->SetTitleOffset(1.6);
  white_noise->GetYaxis()->SetRangeUser(-0.3,0.3);
  white_noise->GetXaxis()->SetRangeUser(0,10);
  cv->SetLogy(0);
  cv->SaveAs("WhiteNoise1.pdf");

  /*
  plot white noise shaped
  */
  white_noise_pulse->SetTitle("");
  white_noise_pulse->SetMarkerStyle(20);
  white_noise_pulse->SetMarkerSize(0.3);
  white_noise_pulse->SetMarkerColor(kBlue);
  white_noise_pulse->SetLineColor(kBlue);
  total_pulse->SetMarkerStyle(20);
  total_pulse->SetMarkerSize(0.3);
  total_pulse->SetMarkerColor(kRed);
  total_pulse->SetLineColor(kRed);
  white_noise_pulse->Draw("APL");
  total_pulse->Draw("PL");
  white_noise_pulse->GetXaxis()->SetTitle("time [ns]");
  white_noise_pulse->GetYaxis()->SetTitle("Amplitude [normalized]");
  white_noise_pulse->GetYaxis()->SetTitleOffset(1.6);
  white_noise_pulse->GetYaxis()->SetRangeUser(-0.1,1.1);
  white_noise_pulse->GetXaxis()->SetRangeUser(0,50);
  cv->SetLogy(0);
  cv->SaveAs("WhiteNoiseShaped1.pdf");

  /*
  write objects to TTree
  */
  pulse->Write();
  // total_pulse->Write("total_pulse");
  // scintillation_pulse->Write("scintillation_pulse");
  // dark_noise->Write("dark_noise");
  // h->Write("h");
  f->Close();

  return 0;
};
