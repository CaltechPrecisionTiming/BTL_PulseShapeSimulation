#include <string>
#include <PulseShape.hh>

PulseShape::PulseShape( double tau, int nf)
{
  //t_sc_random = NULL;
  //t_dc_random = NULL;
  shapingTime_ = tau;
  NFilter_ = nf;

  //Normalize the pulse height to 1.0
  double tmp = 0;
  for (int i=0; i < (100 / 0.01) ; i++ ) {
    double r = ImpulseResponse (i*0.01);
    if (r > tmp) tmp = r;
  }
  ImpulseNormalization_ = tmp;

};

PulseShape::PulseShape( std::string function_name )
{
  //this->function_name = function_name;
  //t_sc_random = NULL;
  //t_dc_random = NULL;
};

PulseShape::PulseShape( std::string function_name, std::string integration_method )
{
  //t_sc_random = NULL;
  //t_dc_random = NULL;
};

PulseShape::~PulseShape()
{
  /*
 if ( t_sc_random != NULL )
 {
   delete [] t_sc_random;
   //t_sc_random = NULL;
   if ( _debug )std::cout <<  "[DEBUG]: deleting memory allocated for SC array" << std::endl;
 }
 if ( t_dc_random != NULL )
 {
   delete [] t_dc_random;
   //t_dc_random = NULL;
   if ( _debug )std::cout <<  "[DEBUG]: deleting memory allocated for DC array" << std::endl;
 }
 */
};

double PulseShape::Gauss( double x, double mean, double sigma, bool norm )
{
  return 5e-3*TMath::Gaus( x, mean, sigma, norm);
};

double PulseShape::Exp( double x, double exponent )
{
  if( x < 0 ) return 0.0;
  return exp(-1.0*exponent*x);
};

double PulseShape::RandomExp( double x, double exponent )
{
  TRandom3 r(0);
  return r.Poisson( 4.5e3*Exp( x, exponent ) );
};

double PulseShape::Convolution( double x, std::string function_name1, std::string function_name2 )
{
  double value = .0;
  double step_size = 1; //nano seconds units of the whole thing
  double x_low  = -1e1;//-1 micro second
  double x_high = 1e3;// +1 micro second
  int steps = int( (x_high-x_low)/step_size );
  if (function_name1 == "Gauss" && function_name2 == "RandomExp")
  {
    //Simpson's rule 1/3
    double h  = step_size/2.0;
    for ( int i = 0; i < steps; i++ )
    {
      double x0 = x_low + i*step_size;
      double x2 = x_low + (i+1)*step_size;
      double x1 = (x0+x2)/2.;
      value += (h/3.)*( RandomExp(x0, 1./40.)*Gauss(x - x0, 10,1) +4.*RandomExp(x1, 1./40.)*Gauss(x - x1, 10,1) + RandomExp(x2, 1./40.)*Gauss(x - x2, 10,1));
    }
  }
  return value;
};

bool SetSinglePhotonResponse( std::string function_name )
{
  //this->function_name = function_name;
  return true;
};

bool SetIntegrationMethod(std::string integration_method )
{
  //this->integration_method = integration_method
  return true;
};

double PulseShape::ScintillationPulse( double x )
{
  if ( Npe <= 0 )
  {
    std::cerr << "[Error] Npe is zero or negative, Npe = " << Npe << std::endl;
    exit(0);
  }
  if ( scintillation_decay_constant <= 0 )
  {
    std::cerr << "[Error] scintillation_decay_constant is zero or negative, scintillation_decay_constant = "
              << scintillation_decay_constant << std::endl;
    exit(0);
  }
  if ( single_photon_risetime_response <= 0 )
  {
    std::cerr << "[Error] single_photon_risetime_response is zero or negative, single_photon_risetime_response = "
              << single_photon_risetime_response << std::endl;
    exit(0);
  }
  if ( single_photon_decaytime_response <= 0 )
  {
    std::cerr << "[Error] single_photon_decaytime_response is zero or negative, single_photon_decaytime_response = "
              << single_photon_decaytime_response << std::endl;
    exit(0);
  }
  //if scintillation times are not yet been drawn then we draw them
  if ( t_sc_random.size() == 0 )
  {
    t_sc_random.clear();
    TRandom3 r(0);//define random variable
    if ( _debug ) std::cout << "[DEBUG] filling vector with containing random times for SC" << std::endl;
    for ( int i = 0; i < Npe; i++ )
    {
      t_sc_random.push_back( r.Exp(scintillation_decay_constant) );
    }
  }
  double eval = 0;
  //t_sc_random.at(0) = 2;
  for ( int i = 0; i < Npe; i++ )
  {
    //eval += TMath::Gaus( x-t_sc_random.at(i), 0, single_photon_response_sigma);
    //eval += 0.5*(x-t_sc_random.at(i))/1.5*exp( -(x-t_sc_random.at(i))/1.5 ) - 0.5*(x-t_sc_random.at(i))/3.*exp( -(x-t_sc_random.at(i))/3.0 );
    if ( x-t_sc_random.at(i) >= 0 )
    {
      /*eval += A*((x-t_sc_random.at(i))/single_photon_risetime_response)*exp( -(x-t_sc_random.at(i))/single_photon_risetime_response )
              - B*(x-t_sc_random.at(i))/single_photon_decaytime_response*exp( -(x-t_sc_random.at(i))/single_photon_decaytime_response );*/
      eval += HighPassFilterResponse(x-t_sc_random.at(i));
    }
  }

  //return eval/single_photon_response_normalization;
  return eval;
};

double PulseShape::DarkNoise( double x, double x_low, double x_high )//Dark Noise in the [x_low, x_high] region, units in ns
{
  int DC = int( DCR*(x_high-x_low) );//number of dark counts in the time window
  //if scintillation times are not yet been drawn then we draw them
  if ( _warning && DC == 0 )
  {
    std::cerr << "[WARNING] DC is zero, are you sure about this? DCR =  " << DCR << std::endl;
  }
  if ( single_photon_risetime_response <= 0 )
  {
    std::cerr << "[Error] single_photon_risetime_response is zero or negative, single_photon_risetime_response = "
              << single_photon_risetime_response << std::endl;
    exit(0);
  }
  if ( single_photon_decaytime_response <= 0 )
  {
    std::cerr << "[Error] single_photon_decaytime_response is zero or negative, single_photon_decaytime_response = "
              << single_photon_decaytime_response << std::endl;
    exit(0);
  }
  if ( t_dc_random.size() == 0 )
  {
    t_dc_random.clear();
    TRandom3 r(0);//define random variable
    if ( _debug ) std::cout << "[DEBUG] filling vector with random times for DC" << std::endl;
    for ( int i = 0; i < DC; i++ )
    {
      t_dc_random.push_back( r.Uniform(x_low,x_high) );
    }
  }
  double eval = 0;
  for ( int i = 0; i < DC; i++ )
  {
    if ( x-t_dc_random.at(i) >= 0 )
    {
      /*eval += A*((x-t_dc_random.at(i))/single_photon_risetime_response)*exp( -(x-t_dc_random.at(i))/single_photon_risetime_response )
            - B*((x-t_dc_random.at(i))/single_photon_decaytime_response)*exp( -(x-t_dc_random.at(i))/single_photon_decaytime_response );*/

      eval += HighPassFilterResponse(x-t_dc_random.at(i));
    }
  }
  return eval/single_photon_response_normalization;

};


/*
4-point LGAD response.
Need to promote this to Nicolo's code
*/
double PulseShape::LGADPulse( double x )
{
  //x is assumed to be in units of ns
  double eval = 0;

  //4-point signal from Gregory Deptuch
  if (x >= 0 && x<0.2) eval = (0.8 / 0.2) * x ;
  else if (x >= 0.2 && x<0.7) eval = 0.8 - (0.1 / 0.5)*(x - 0.2);
  else if (x>=0.7 && x < 1.5) eval = 0.7 - (0.7 / 0.8)*(x - 0.7);
  else eval = 0;

  //Delta Function
  //if (x==0) eval = 1;

  return eval;
};


/*
This is the convolution using the Simpson rule, its now taking the semi-gaussian impulse ImpulseResponse
and then convoluting now with the 4-point LGAD response.
*/
double PulseShape::LGADShapedPulse( double x )
{
  double eval = 0;

  const double integrationStep = 1.0; //in ns
  const double integrationLow = -100;
  const double integrationHigh = 100;
  const int NIntegrationPoints = (integrationHigh - integrationLow) / integrationStep;
  // for (int i=0; i < NIntegrationPoints; i++) {
  //   double s = integrationLow + i * integrationStep;
  //   eval += LGADPulse(s) * NormalizedImpulseResponse(x-s) * integrationStep;
  // }
  //eval = eval / (integrationHigh - integrationLow);

  //Simpson's rule 1/3
  double h  = integrationStep/2.0;
  for ( int i = 0; i < NIntegrationPoints; i++ ) {
      double x0 = integrationLow + i*integrationStep;
      double x2 = integrationLow + (i+1)*integrationStep;
      double x1 = (x0+x2)/2.;
      eval += (h/3.)*( LGADPulse(x0)*NormalizedImpulseResponse(x-x0)
			+ 4.0 * LGADPulse(x1)*NormalizedImpulseResponse(x-x1)
			+ LGADPulse(x2)*NormalizedImpulseResponse(x-x2)
			);
  }




  //return eval/single_photon_response_normalization;
  return eval;
};


void PulseShape::NormalizeSinglePhotonResponse()
{
  double x_low  = .0;//ns
  double x_high = 200;//ns
  double step = 1e-3; //1ps
  const int n_iterations = (x_high-x_low)/step;
  double max_val = 0;
  for ( int i = 0; i < n_iterations; i++ )
  {
    double x = x_low + double(i)*step;
    double f_x =  A*(x/single_photon_risetime_response)*exp( -x/single_photon_risetime_response )
                - B*(x/single_photon_decaytime_response)*exp( -x/single_photon_decaytime_response );
    if ( f_x > max_val ) max_val = f_x;
  }

  if ( max_val > 0 )
  {
    single_photon_response_normalization = max_val;
    return;
  }
  else
  {
    std::cerr << "[ERROR] single_photon_response_normalization zero or negative\nEXIT" << std::endl;
    exit(0);
  }

};

double PulseShape::ImpulseResponse( double x )
{
  double eval = 0;
  double omegashaper = NFilter_ / shapingTime_;

  if (x>=0) {
    eval = exp(-omegashaper*x) * pow(x,NFilter_);
  } else {
    eval = 0;
  }

  //First attempt, hard coded from Mathematica notebook
  //NFilter = 2; shaping time = 4ns ; x is in units of ns
  // if (x >= 0) eval = 4.61816e-1 * exp(-0.5*x)*x*x;
  // else eval = 0;

  return eval;
};

double PulseShape::NormalizedImpulseResponse( double x )
{
  return ImpulseResponse(x) / ImpulseNormalization_;
};

double PulseShape::HighPassFilterResponse( double x )
{
  double eval = high_pass_filter_RC*(
                A*single_photon_risetime_response*exp(-x/single_photon_risetime_response)/pow(single_photon_risetime_response-high_pass_filter_RC,2)
                +A*x*exp(-x/single_photon_risetime_response)/(high_pass_filter_RC*single_photon_risetime_response-pow(single_photon_risetime_response,2))
                -B*single_photon_decaytime_response*exp(-x/single_photon_decaytime_response)/pow(single_photon_decaytime_response-high_pass_filter_RC,2)
                -B*x*exp(-x/single_photon_decaytime_response)/(high_pass_filter_RC*single_photon_decaytime_response-pow(single_photon_decaytime_response,2))
                +(2*(A-B)*high_pass_filter_RC*single_photon_risetime_response*single_photon_decaytime_response
                +single_photon_risetime_response*single_photon_decaytime_response*(B*single_photon_risetime_response-A*single_photon_decaytime_response)
                +pow(high_pass_filter_RC,2)*(B*single_photon_decaytime_response-A*single_photon_risetime_response))*
                exp(-x/high_pass_filter_RC)/(pow(single_photon_risetime_response-high_pass_filter_RC,2)*pow(single_photon_decaytime_response-high_pass_filter_RC,2))
              )/single_photon_response_normalization;
  return eval;
};

void PulseShape::NormalizeSinglePhotonResponseHighPassFilter()
{
  double x_low  = .0;//ns
  double x_high = 200;//ns
  double step = 1e-3; //1ps
  const int n_iterations = (x_high-x_low)/step;
  double max_val = 0;
  for ( int i = 0; i < n_iterations; i++ )
  {
    double x = x_low + double(i)*step;
    double f_x =  high_pass_filter_RC*(
                  A*single_photon_risetime_response*exp(-x/single_photon_risetime_response)/pow(single_photon_risetime_response-high_pass_filter_RC,2)
                  +A*x*exp(-x/single_photon_risetime_response)/(high_pass_filter_RC*single_photon_risetime_response-pow(single_photon_risetime_response,2))
                  -B*single_photon_decaytime_response*exp(-x/single_photon_decaytime_response)/pow(single_photon_decaytime_response-high_pass_filter_RC,2)
                  -B*x*exp(-x/single_photon_decaytime_response)/(high_pass_filter_RC*single_photon_decaytime_response-pow(single_photon_decaytime_response,2))
                  +(2*(A-B)*high_pass_filter_RC*single_photon_risetime_response*single_photon_decaytime_response
                  +single_photon_risetime_response*single_photon_decaytime_response*(B*single_photon_risetime_response-A*single_photon_decaytime_response)
                  +pow(high_pass_filter_RC,2)*(B*single_photon_decaytime_response-A*single_photon_risetime_response))*
                  exp(-x/high_pass_filter_RC)/(pow(single_photon_risetime_response-high_pass_filter_RC,2)*pow(single_photon_decaytime_response-high_pass_filter_RC,2))
                );
    if ( f_x > max_val ) max_val = f_x;
  }

  if ( max_val > 0 )
  {
    single_photon_response_normalization = max_val;
    return;
  }
  else
  {
    std::cerr << "[ERROR] HighPassFilter zerp or negative\nEXIT" << std::endl;
    exit(0);
  }
};
