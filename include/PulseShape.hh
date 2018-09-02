#ifndef PulseShape_HH
#define PulseShape_HH

#include <iostream>
#include <string>
//ROOT
#include <TF1.h>
#include <TROOT.h>
#include <TRandom3.h>
#include <TH1F.h>
#include <TMath.h>
#include <TGraph.h>
#include <TString.h>
#include <TTree.h>
#include <vector>

class PulseShape
{
public:

  static bool _debug;
  static bool _info;
  static bool _warning;

  PulseShape( double tau , int nf, float SNR, int seed);
  PulseShape( double tau , int nf, float SNR, int seed, std::vector<std::vector<std::pair<double,double> > > &LGADPulseLibrary);
  PulseShape( std::string function_name );
  PulseShape( std::string function_name, std::string integration_method );
  ~PulseShape();

  double Gauss( double x, double mean, double sigma, bool norm = true );
  double Exp( double x, double exponent );
  double RandomExp( double x, double exponent );

  double Convolution( double x, std::string function_name1, std::string function_name2 );

  //ETL_ASIC PART
  double LGADPulse( double x );
  double LGADShapedPulse( double x );
  //NOISE
  double WhiteNoise(double mean, double rms);
  double WhiteNoiseShapedPulse( double x, double mean, double rms );
  //
  double ImpulseResponse( double x );
  double NormalizedImpulseResponse( double x );
  bool SetSinglePhotonResponse( std::string function_name );
  bool SetIntegrationMethod(std::string integration_method );
  void SetNpe( int npe ){ Npe = npe;};
  void SetDCR( double dcr ){ DCR = dcr;};//in GHz
  void SetSinglePhotonResponse( double sigma ){ single_photon_response_sigma = sigma;};//units in ns
  void SetSinglePhotonRisetimeResponse( double risetime ){ single_photon_risetime_response = risetime;};//units in ns
  void SetSinglePhotonDecaytimeResponse( double decaytime ){ single_photon_decaytime_response = decaytime;};//units in ns
  void SetScintillationDecay( double tau_s ){ scintillation_decay_constant = tau_s;};//units in ns
  void SetHighPassFilterRC( double rc ){ high_pass_filter_RC = rc;};
  double GetSinglePhotonResponseNormalization(){return single_photon_response_normalization;};

protected:
  std::string function_name;
  std::string integration_method;
  int Npe;// Number of photo-electrons (dE/dx*thickness*LightCollectionEfficiency*SiPM_PDE), about 4500 in LYSO
  double scintillation_decay_constant;//decay constant of the scintillator (LYSO is 40 ns )
  double scintillation_risetime;//rise time of the scintillator (LYSO is 60 ps)
  double single_photon_response_sigma;//sigma of the gaussian used to model the single photon response
  double single_photon_risetime_response;//tau1 in of A*t/tau1*exp(-t/tau1) - B*t/tau2*exp(-t/tau2) used to model the single photon response
  double single_photon_decaytime_response;//tau2 in of A*t/tau1*exp(-t/tau1) - B*t/tau2*exp(-t/tau2) used to model the single photon response
  double high_pass_filter_RC;//tau of the RC HighPassFilter (R*C) in ns
  double DCR;//dark count rate in GHz
  double single_photon_response_normalization;
  double *noise;
  double *LGADSignal;

  //LGAD parameters
  bool useLGADLibrary_ = false;
  double shapingTime_;
  int NFilter_;
  double ImpulseNormalization_;

  //internal parameters 
  float SNR_;
  float integrationWindowLow_;
  float integrationWindowHigh_;
  int NIntegrationPoints_;
  
  int NIntegrationPointsLGADSignal_;
  float IntegralTimeStepSignal_;

  TRandom3 *random_;
  int randomSeed_;
  int randomSignalEvent_;
  std::vector<double> t_sc_random;
  std::vector<double> t_dc_random;
  const double A = 3.0;
  const double B = 0.5;

};

#endif
