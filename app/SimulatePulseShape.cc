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
  

  std::cout << "number of experiments is: " << n_experiments << std::endl;
  std::cout << "NFilter: " << NFilter << std::endl;
  std::cout << "ShapingTime: " << ShapingTime << " ns" << std::endl;
  std::cout << "Signal-to-Noise Ratio: " << SNR << std::endl;


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
    for( int i = 0; i < npoints; i++ ) y[i] = x[i] = 0.0;
    double y_max = 0;

    //create pulse shape object
    ps = new PulseShape( ShapingTime, NFilter );
    
    
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
      y_sc[i]  = ps->LGADShapedPulse(x[i]) / normalization;
      //y_sc[i]  = ps->NormalizedImpulseResponse(x[i]);

      //cout << i << " : " << y_sc[i] << "\n";
      // y_dc[i]  = ps->DarkNoise(x[i], x_low, x_high);
      y[i]     = y_sc[i];// + y_dc[i];
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
