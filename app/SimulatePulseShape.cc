#include <iostream>
#include <string>
#include <vector>
//ROOT
#include <TCanvas.h>
#include <TFile.h>
#include <TTree.h>
#include <TH1F.h>
//LOCAL
#include <PulseShape.hh>
#include <Configuration.hh>
#include <assert.h>

bool PulseShape::_info    = false;
bool PulseShape::_debug   = false;
bool PulseShape::_warning = false;

int main ( int argc, char** argv )
{

  const float gain = 1e7;
  TH1F* h = new TH1F("pulse_time", "pulse_time", 2000, -10,10);
  TH1F* input_voltage = new TH1F("voltage", "voltage", 100, 0,1-3);
  TH1F* voltage_rms = new TH1F("voltage_rms", "voltage_rms", 100, 0,1-3);

  TTree* pulse = new TTree("pulse", "Digitized waveforms");

  Configuration* config = new Configuration(true);
  config->GetCommandLineArgs(argc, argv);
  config->ParseConfigurationFile(config->config_file);

  const int NFilter = config->NFilter;
  const int n_experiments = config->n_experiments;
  const double ShapingTime = config->ShapingTime;
  const double SNR = config->SNR;
  double SignalAmplitudeMean = config->SignalAmplitudeMean;
  double NoiseRMS = SignalAmplitudeMean / SNR;
  const int randomSeed = config->randomSeed;
  string LGADSignalFilename = config->LGADSignalFilename;
  TRandom3 *random = new TRandom3(randomSeed);

  bool useLGADPulseLibrary = false;
  TFile *LGADSignalFile = 0;
  TTree *LGADSignalTree = 0;
  if (LGADSignalFilename != "") {
    LGADSignalFile = new TFile(LGADSignalFilename.c_str(),"READ");
    if (!LGADSignalFile) { cout << "Error: LGAD Signal File " << LGADSignalFilename << " was not opened successfully. \n"; assert(false);}
    LGADSignalTree = (TTree*)LGADSignalFile->Get("pulse");
    assert(LGADSignalTree);
    useLGADPulseLibrary = true;
  }

  //******************************************************************
  //Load the LGAD Signal Pulses into memory
  //******************************************************************
  std::vector<std::vector< std::pair<double,double > > > LGADSignalLibrary;
  if (useLGADPulseLibrary) {
    std::cout << "Loading LGAD Signal Pulses into memory from : " << LGADSignalFilename << "\n";
    const int npointsSignal = 1500;
    const double impedance = 50; //use 50 Ohm impedance
    float tmpTime[npointsSignal];
    float tmpAmp[npointsSignal];
    LGADSignalTree->SetBranchAddress("time",&tmpTime);
    LGADSignalTree->SetBranchAddress("channel",&tmpAmp);
    for( int i=0; i < LGADSignalTree->GetEntries();i++) {
      //cout << "reading signal event: " << i << "\n";
      std::vector< std::pair <double,double > > pulse;
      LGADSignalTree->GetEntry(i);
      double max_voltage = 0;
      double signal_rms = 0;
      for (int j=0; j < npointsSignal; j++) {
        double current_voltage = impedance * tmpAmp[j]*gain;
        pulse.push_back( std::pair<double,double>( tmpTime[j] , current_voltage ));
        if ( max_voltage < impedance * tmpAmp[j] ) max_voltage = current_voltage;
        if ( j < 1100 ) signal_rms = pow(current_voltage, 2.0);
      }
      LGADSignalLibrary.push_back(pulse);
      input_voltage->Fill(max_voltage);
      voltage_rms->Fill(sqrt(signal_rms/1100.));
      //std::cout << max_voltage << std::endl;
    }
    LGADSignalFile->Close();
  }

  TF1* landau = new TF1("landau", "landau", 0.4*input_voltage->GetMean(), 1.5*input_voltage->GetMean());
  input_voltage->Fit("landau", "LRW");
  double mpv_max_voltage = landau->GetParameter(1);
  //NoiseRMS = SignalAmplitudeMean / SNR;
  //****************
  //RMS
  //****************
  TF1* landau_rms = new TF1("landau_rms", "landau", 0.2*voltage_rms->GetMean(), 1.6*voltage_rms->GetMean());
  voltage_rms->Fit("landau_rms", "LRW");
  double  mpv_voltage_rms = landau_rms->GetParameter(1);

  std::cout << "mpv max voltage: " << mpv_max_voltage << " , mpv rms voltage: " <<  mpv_voltage_rms << std::endl;

  SignalAmplitudeMean = mpv_voltage_rms;
  NoiseRMS = SignalAmplitudeMean / SNR;
  //std::cout << "integral: " << input_voltage->GetMean() << std::endl;

  /*
  auto sort_lgad_pulses = []( std::pair<double,double > a, std::pair<double,double > b){ return a.second() > b.second() ?  true : false; };
  std::sort( LGADSignalLibrary.at(0).begin(), LGADSignalLibrary.at(0).end(), sort_lgad_pulse );
  */
  std::cout << "number of experiments is: " << n_experiments << std::endl;
  std::cout << "NFilter: " << NFilter << std::endl;
  std::cout << "ShapingTime: " << ShapingTime << " ns" << std::endl;
  std::cout << "Signal-to-Noise Ratio: " << SNR << std::endl;
  std::cout << "SignalAmplitudeMean: " << SignalAmplitudeMean << std::endl;
  std::cout << "NoiseRMS: " << NoiseRMS << std::endl;


  //******************************************************************
  //Set up output tree
  //******************************************************************
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
  float y[npoints], y_signal[npoints], y_dc[npoints], y_wnps[npoints];
  double* noise_v2 = new double[npoints];
  double noise[npoints];
  int i_evt;
  pulse->Branch("i_evt", &i_evt, "i_evt/i");
  pulse->Branch("channel", y, Form("channel[%d]/F", npoints));
  pulse->Branch("noise", &noise[0], Form("noise[%d]/D", npoints));
  pulse->Branch("shapednoise", y_wnps, Form("shapednoise[%d]/F", npoints));
  pulse->Branch("time", x, Form("time[%d]/F", npoints));
  for ( int j = 0; j < n_experiments; j++ )
  {
    if ( j % 1 == 0 )std::cout << "Generating Event #" << j << std::endl;
    //reset variables and objects
    for( int i = 0; i < npoints; i++ ) y[i] = x[i] = 0.0;
    double y_max = 0;

    //create pulse shape object
    double normalization = 0;
    if (useLGADPulseLibrary)
    {
      ps = new PulseShape( ShapingTime, NFilter , NoiseRMS, randomSeed+j, LGADSignalLibrary );
      normalization = 1; //don't do normalization for this case
    }
    else
    {
      ps = new PulseShape( ShapingTime, NFilter , NoiseRMS, randomSeed+j );

      //normalize pulse shape to have pulse height at 1.0V
      for( int i = 0; i < int(100 / 0.01); i++ )
      {
        //std::cout << "norm " << i << " ";
        double tmp = ps->LGADShapedPulse( i * 0.01);
        //std::cout << tmp << "\n";
        if ( tmp > normalization ) normalization = tmp;
      }
      //introduce a 20% fluctuation on the normalization
      normalization = normalization / (1.0 + random->Gaus(0,0.2));
    }

    //populate the pulse shape
    for( int i = 0; i < npoints; i++ )
    {
      x[i] = x_low + double(i)*step;
      ////if ( i % 1000 == 0 ) std::cout << "iteration #" << i << std::endl;
      y_signal[i]   = ps->LGADShapedPulse(x[i]) / normalization;
      y_wnps[i] = ps->WhiteNoiseShapedPulse(x[i]);
      //y_wnps[i] = ps->WhiteNoiseShapedPulse_ZOH(x[i]);
      //y_wnps[i] = ps->DiscriteWhiteNoiseShapedPulse(i);
      //cout << i << " : " << y_signal[i] << "\n";
      y[i]     = y_signal[i] + y_wnps[i];// + y_dc[i];
      if( y[i] > y_max ) y_max = y[i];
    }



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
    noise_v2 = ps->GetNoiseArray();
    double input_noise_rms_simulation = 0;
    double output_noise_rms_simulation = 0;
    for( int i = 0; i < npoints; i++)
    {
      noise[i] = noise_v2[i];
      //std::cout << "here: " << noise[i] << std::endl;
      input_noise_rms_simulation += noise[i]*noise[i];
      output_noise_rms_simulation += y_wnps[i]*y_wnps[i];
    }
    input_noise_rms_simulation = sqrt(input_noise_rms_simulation/npoints);
    output_noise_rms_simulation = sqrt(output_noise_rms_simulation/npoints);
    std::cout << "input noise rms: " << input_noise_rms_simulation << std::endl;
    std::cout << "output noise rms: " << output_noise_rms_simulation << std::endl;
    pulse->Fill();
    delete ps;//release memory of pulseshape object.
  }

  ps = NULL;

  TFile* f = new TFile(config->output_file.c_str(), "recreate");
  total_pulse = new TGraph( npoints, x, y);
  scintillation_pulse = new TGraph( npoints, x, y_signal);
  //dark_noise = new TGraph( npoints, x, noise);
  //white_noise = new TGraph( npoints, x, noise);
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
  /*
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
  */

  /*
  plot white noise shaped
  */
  /*
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
  */
  input_voltage->Draw("HISTO");
  landau->Draw("same");
  cv->SaveAs("voltage.pdf");

  voltage_rms->Draw("HISTO");
  landau_rms->Draw("same");
  cv->SaveAs("voltage_rms.pdf");

  /*
  write objects to TTree
  */
  pulse->Write();
  //total_pulse->Write("total_pulse");
  //scintillation_pulse->Write("scintillation_pulse");
  //dark_noise->Write("dark_noise");
  h->Write("h");
  input_voltage->Write("voltage");
  voltage_rms->Write("voltage_rms");
  f->Close();

  return 0;
};
