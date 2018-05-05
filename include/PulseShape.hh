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

class PulseShape
{
public:

  PulseShape();
  PulseShape( std::string function_name );
  PulseShape( std::string function_name, std::string integration_method );
  ~PulseShape();

  double Gauss( double x, double mean, double sigma, bool norm = true );
  double Exp( double x, double exponent );
  double RandomExp( double x, double exponent );

  double Convolution( double x, std::string function_name1, std::string function_name2 );
  bool SetSinglePhotonResponse( std::string function_name );
  bool SetIntegrationMethod(std::string integration_method );

protected:
  std::string function_name;
  std::string integration_method;
  double Npe;
  double scintillation_decay_constant;
  double scintillation_risetime;
  
};

#endif
