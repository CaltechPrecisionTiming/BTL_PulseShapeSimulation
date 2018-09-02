#include <string>
#include <chrono>
#include <random>
#include <assert.h>
#include <PulseShape.hh>

PulseShape::PulseShape( double tau, int nf, float SNR, int seed, std::vector<std::vector<std::pair<double,double> > > &LGADPulseLibrary)
{
  //t_sc_random = NULL;
  //t_dc_random = NULL;
  shapingTime_ = tau;
  NFilter_ = nf;
  integrationWindowLow_ = -100;
  integrationWindowHigh_ = 100;
  double tmpNSteps = (integrationWindowHigh_ - integrationWindowLow_) / (shapingTime_ / 10);
  NIntegrationPoints_ = int(tmpNSteps);
  if (NIntegrationPoints_ != tmpNSteps) {
    std::cout << "Error: integration window from " << integrationWindowLow_ << " to " << integrationWindowHigh_ 
  	 << " is not divisible by shapingTime_/10 = " << (shapingTime_ / 10) 
  	 << "\n";
    assert(false);
  } 



  //Normalize the pulse height to 1.0
  double tmp = 0;
  for (int i=0; i < (100 / 0.01) ; i++ ) {
    double r = ImpulseResponse (i*0.01);
    if (r > tmp) tmp = r;
  }
  ImpulseNormalization_ = tmp;

  //Create TRandom3 object
  randomSeed_ = seed;
  random_ = new TRandom3(randomSeed_);

  //**********************************************
  //Randomly pick a signal pulse from the library
  //**********************************************
  useLGADLibrary_ = true;
  randomSignalEvent_ = random_->Integer(LGADPulseLibrary.size());
  const int npoints = 1500;
  IntegralTimeStepSignal_ = shapingTime_ / 100;

  double tmpNStepsLGADSignal = (integrationWindowHigh_ - integrationWindowLow_) / IntegralTimeStepSignal_;
  NIntegrationPointsLGADSignal_ = int(tmpNStepsLGADSignal);
  if (NIntegrationPointsLGADSignal_ != tmpNStepsLGADSignal) {
    std::cout << "Error: integration window from " << integrationWindowLow_ << " to " << integrationWindowHigh_ 
	      << " is not divisible by shapingTime_/10 = " << IntegralTimeStepSignal_
	      << "\n";
    assert(false);
  } 
  LGADSignal = new double[NIntegrationPointsLGADSignal_];

  //std::cout << "randomSignalEvent_ : " << randomSignalEvent_ << " \n";
  for ( int i  = 0; i < NIntegrationPointsLGADSignal_; i++ ) {
    int tmpSignalPulseIndex =  i*(1000 * IntegralTimeStepSignal_);
    if (tmpSignalPulseIndex < 1500) {
      LGADSignal[i] = LGADPulseLibrary[randomSignalEvent_][tmpSignalPulseIndex].second;
    } else {
      LGADSignal[i] = 0;
    }
    //std::cout << i << " : " << tmpSignalPulseIndex << " -> " << LGADSignal[i] << "\n";
  }
  

  //****************************
  //Simulate the Noise
  //****************************
  SNR_ = SNR;
  noise = new double[NIntegrationPoints_];
  for ( int i  = 0; i < NIntegrationPoints_; i++ )
  {
    noise[i] = WhiteNoise(0,1./SNR_);
  }

};

PulseShape::PulseShape( double tau, int nf, float SNR, int seed)
{
  //t_sc_random = NULL;
  //t_dc_random = NULL;
  useLGADLibrary_ = false;
  shapingTime_ = tau;
  NFilter_ = nf;
  integrationWindowLow_ = -100;
  integrationWindowHigh_ = 100;
  double tmpNSteps = (integrationWindowHigh_ - integrationWindowLow_) / (shapingTime_ / 10);
  NIntegrationPoints_ = int(tmpNSteps);
  if (NIntegrationPoints_ != tmpNSteps) {
    std::cout << "Error: integration window from " << integrationWindowLow_ << " to " << integrationWindowHigh_ 
  	 << " is not divisible by shapingTime_/10 = " << (shapingTime_ / 10) 
  	 << "\n";
    assert(false);
  } 

  //Normalize the pulse height to 1.0
  double tmp = 0;
  for (int i=0; i < (100 / 0.01) ; i++ ) {
    double r = ImpulseResponse (i*0.01);
    if (r > tmp) tmp = r;
  }
  ImpulseNormalization_ = tmp;

  //Create TRandom3 object
  randomSeed_ = seed;
  random_ = new TRandom3(randomSeed_);

  SNR_ = SNR;
  noise = new double[NIntegrationPoints_];
  for ( int i  = 0; i < NIntegrationPoints_; i++ )
  {
    noise[i] = WhiteNoise(0,1./SNR_);
  }

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


/*
4-point LGAD response.
Need to promote this to Nicolo's code
*/
double PulseShape::LGADPulse( double x )
{
  const double timeShift = 20;
  double t = x - timeShift;
  //t is assumed to be in units of ns
  double eval = 0;  
  //std::cout << "LGADPulse " << useLGADLibrary_ << "\n";
  if (useLGADLibrary_) {
    //std::cout << "Pulse " << t ;
    if (t < 0 || t > 1.5) {
      eval = 0;
    }
    else {
      int tmpSignalIndex = int ( std::round ( t / IntegralTimeStepSignal_ ));
      //std::cout << t / IntegralTimeStepSignal_ << " " << int ( std::round ( t / IntegralTimeStepSignal_ )) ;
      eval = LGADSignal[tmpSignalIndex];
    }
    //std::cout << " ; \n";
  } else {    
    //4-point signal from Gregory Deptuch
    if (t >= 0 && t<0.2) eval = (0.8 / 0.2) * t ;
    else if (t >= 0.2 && t<0.7) eval = 0.8 - (0.1 / 0.5)*(t - 0.2);
    else if (t>=0.7 && t < 1.5) eval = 0.7 - (0.7 / 0.8)*(t - 0.7);
    else eval = 0;
  }

  //Delta Function
  //if (t==0) eval = 1;

  return eval;
};

double PulseShape::WhiteNoise( double mean, double rms )
{
  //x is assumed to be in units of ns

  // construct a trivial random generator engine from a time-based seed:
  // unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
  // std::default_random_engine generator (seed);  
  // std::normal_distribution<double> distribution (mean, 
  // return distribution(generator);

  //Use the TRandom3 version of this
  return random_->Gaus(mean, rms);


};

/*
This is the convolution using the Simpson rule, its now taking the semi-gaussian impulse ImpulseResponse
and then convoluting now with the 4-point LGAD response.
*/
double PulseShape::LGADShapedPulse( double x )
{
  double eval = 0;

  const double integrationStep = (shapingTime_ / 10); //in ns

  //Use Simpson's rule 
  double h  = integrationStep/2.0;
  //std::cout << "LGADShapedPulse " << NIntegrationPoints_ << "\n";
  for ( int i = 0; i < NIntegrationPoints_; i++ ) {    
      double x0 = integrationWindowLow_ + i*integrationStep;
      double x2 = integrationWindowLow_ + (i+1)*integrationStep;
      double x1 = (x0+x2)/2.;
      //std::cout << "integral " << i << " : " << x0 << " " << x1 << " " << x2 << " : " << LGADPulse(x0) << "\n";
      eval += (h/3.)*( LGADPulse(x0)*NormalizedImpulseResponse(x-x0)
			+ 4.0 * LGADPulse(x1)*NormalizedImpulseResponse(x-x1)
			+ LGADPulse(x2)*NormalizedImpulseResponse(x-x2)
			);
  }

  return eval;
};

double PulseShape::WhiteNoiseShapedPulse( double x, double mean, double rms )
{
  double eval = 0;

  const double integrationStep = (shapingTime_ / 10);
  for (int i=0; i < NIntegrationPoints_; i++) {
     double s = integrationWindowLow_ + i * integrationStep;
     eval += noise[i] * NormalizedImpulseResponse(x-s) * integrationStep;
   }

  return eval;
};


double PulseShape::ImpulseResponse( double x )
{
  const double timeShift = 0.0;
  double eval = 0;
  double omegashaper = NFilter_ / shapingTime_;

  if (x>=0) {
    eval = exp(-omegashaper*(x-timeShift)) * pow(x-timeShift,NFilter_);
  } else {
    eval = 0;
  }


  return eval;
};

double PulseShape::NormalizedImpulseResponse( double x )
{
  return ImpulseResponse(x) / ImpulseNormalization_;
};


