#include <string>
#include <chrono>
#include <random>
#include <assert.h>
#include <PulseShape.hh>

PulseShape::PulseShape( double tau, int nf, float NoiseRMS, int seed, std::vector<std::vector<std::pair<double,double> > > &LGADPulseLibrary)
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
  //ImpulseNormalization_ = tmp;
  ImpulseNormalization_ = 1.0;

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
  npoints_noise_ = int((integrationWindowHigh_ - integrationWindowLow_)/integration_noise_step_);
  noise = new double[npoints_noise_];
  for ( int i  = 0; i < npoints_noise_; i++ )
  {
    noise[i] = WhiteNoise(0,NoiseRMS);
  }


};

PulseShape::PulseShape( double tau, int nf, float NoiseRMS, int seed)
{
  //t_sc_random = NULL;
  //t_dc_random = NULL;
  useLGADLibrary_ = false;
  shapingTime_ = tau;
  NFilter_ = nf;
  integrationWindowLow_ = -100;
  integrationWindowHigh_ = 100;
  double tmpNSteps = (integrationWindowHigh_ - integrationWindowLow_) / (shapingTime_ / 10.);
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
  //ImpulseNormalization_ = tmp;
  ImpulseNormalization_ = 1.0;

  //Create TRandom3 object
  randomSeed_ = seed;
  random_ = new TRandom3(randomSeed_);

  //****************************
  //Simulate the Noise
  //****************************
  //const int npoints_noise = int((integrationWindowHigh_ - integrationWindowLow_)/integration_noise_step_);
  noise = new double[npoints_noise_];
  for ( int i  = 0; i < npoints_noise_; i++ )
  {
    noise[i] = WhiteNoise(0,NoiseRMS);
  }

  //****************************
  //Simulate the Noise
  //****************************
  random_p1_time = random_->Gaus(0, 0.1);
  random_p2_time = random_->Gaus(0, 0.1);
  random_p3_time = random_->Gaus(0, 0.1);
  random_p4_time = random_->Gaus(0, 0.1);
  random_p1_amp = random_->Gaus(0, 0.1);
  random_p2_amp = random_->Gaus(0, 0.1);
  random_p3_amp = random_->Gaus(0, 0.1);
  std::cout << "random: " << random_p1_time << "\n";

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
    double p1_time = 0.05*(1+random_p1_time);
    double p2_time = 0.2*(1+random_p2_time);
    double p3_time = 0.7*(1+random_p3_time);
    double p4_time = 1.5*(1+random_p4_time);
    double p1_amp = 0.2*(1+random_p1_amp);
    double p2_amp = 0.8*(1+random_p1_amp);
    double p3_amp = 0.7*(1+random_p1_amp);
    const double p4_amp = 0.0;


    if (t >= 0 && t<p1_time) eval = (p1_amp/p1_time) * t ;
    else if (t >= p1_time && t<p2_time) eval = p1_amp + ((p2_amp-p1_amp) / (p2_time-p1_time)) * (t - p1_time) ;
    else if (t >= p2_time && t<p3_time) eval = p2_amp + ((p3_amp-p2_amp) / (p3_time-p2_time)) * (t - p2_time) ;
    else if (t >= p3_time && t<p4_time) eval = p3_amp + ((p4_amp-p3_amp) / (p4_time-p3_time)) * (t - p3_time) ;
    else eval = 0;

    //Fixed Pulse from Deptuch
    // if (t >= 0 && t<0.05) eval = (0.8 / 0.2) * t ;
    // else if (t >= 0 && t<0.2) eval = (0.8 / 0.2) * t ;
    // else if (t >= 0.2 && t<0.7) eval = 0.8 - (0.1 / 0.5)*(t - 0.2);
    // else if (t>=0.7 && t < 1.5) eval = 0.7 - (0.7 / 0.8)*(t - 0.7);
    // else eval = 0;
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
      eval += (h/3.)*(
        LGADPulse(x0)*NormalizedImpulseResponse(x-x0) +
        4.0 * LGADPulse(x1)*NormalizedImpulseResponse(x-x1) +
        LGADPulse(x2)*NormalizedImpulseResponse(x-x2)
			);
  }

  return eval;
};

double PulseShape::WhiteNoiseShapedPulse( double x )
{
  double eval = 0;
  double h = integration_noise_step_/2.0;
  for (int i=0; i < npoints_noise_; i++)
  {
    double x0 = integrationWindowLow_ + i * integration_noise_step_;
    double x2 = integrationWindowLow_ + (i+1)*integration_noise_step_;
    double x1 = (x0+x2)/2.;
     //double s_1 = integrationWindowLow_ + (i+1) * integration_noise_step_;
     //eval += noise[i] * NormalizedImpulseResponse(x-s) * integrationStep;
     //eval += noise[i] * NormalizedImpulseResponse(x-s);//white noise does not need integration step, cancels out when taking the delta function integral
    eval += (h/3.)*(
      noise[i] * NormalizedImpulseResponse(x-x0) +
      (4./2.)*(noise[i]+ noise[i+1])* NormalizedImpulseResponse(x-x1) +
      noise[i+1] * NormalizedImpulseResponse(x-x2)
    );
   }

  return eval;
};

double PulseShape::WhiteNoiseShapedPulse_ZOH( double x )
{
  double eval = 0;
  for (int i=0; i < npoints_noise_; i++)
  {
    eval += noise[i] * ZOH_Response(x,integration_noise_step_, i);
  }

  return eval;
};

double PulseShape::DiscriteWhiteNoiseShapedPulse( int i )
{
  double eval = 0;
  double samples_impulse_response[npoints_noise_];
  for (int j=0; j < npoints_noise_; j++) {
     double s = integrationWindowLow_ + j * integration_noise_step_;
     samples_impulse_response[i] = NormalizedImpulseResponse(s);
     //eval += noise[i] * NormalizedImpulseResponse(x-s) * integrationStep;
     //eval += noise[i] * NormalizedImpulseResponse(x-s);//white noise does not need integration step, cancels out when taking the delta function integral
   }
   for ( int j = 0; j < npoints_noise_; j++)
   {
    eval += noise[j] * samples_impulse_response[i-j];
   }

  return eval;
};


double PulseShape::ImpulseResponse( double x )
{
  //const double timeShift = 0.0;
  double eval = 0;
  double omegashaper = NFilter_ / shapingTime_;

  if (x>=0) {
    eval = exp(-omegashaper*x) * pow(x,NFilter_)*(1/pow(shapingTime_,NFilter_+1));
  } else {
    eval = 0;
  }
  return eval;
};

double PulseShape::ZOH_Response(double x, double T, int k){
  const double timeShift = 0.0;
  double eval = 0;
  double omegashaper = 1. / shapingTime_;

  if (x>=0)
  {
    if ( x >= k*T && x < (k+1)*T )
    {
      eval = -(1/pow(shapingTime_,2.))*shapingTime_*(-shapingTime_ + exp(-x/shapingTime_)*exp(k*T/shapingTime_)*(x - k*T + shapingTime_)) ;
    }
    else if ( x >= (k+1)*T )
    {
      eval = -(1/pow(shapingTime_,2.))*shapingTime_*(
        ( -shapingTime_ + exp(-x/shapingTime_)*exp(k*T/shapingTime_)*(x - k*T + shapingTime_))+
        ( exp(-x/shapingTime_)*exp((k+1)*T/shapingTime_)*((k+1)*T-x-shapingTime_) +shapingTime_ ));
    }
    else
    {
      eval = 0;
    }
  }
  else
  {
    eval = 0;
  }

  return eval;
};

double PulseShape::NormalizedImpulseResponse( double x )
{
  return ImpulseResponse(x) / ImpulseNormalization_;
};

double* PulseShape::GetNoiseArray()
{
  double* return_noise = new double[npoints_noise_];
  for ( int i = 0; i < npoints_noise_; i++ ) return_noise[i] = noise[i];
  return return_noise;
};

float PulseShape::FrequencySpectrum(double freq, double tMin, double tMax, unsigned int n_samples, float* my_channel, float* my_time)
{
  const int range = 0; // extension of samples to be used beyond [tMin, tMax]
	double deltaT = (my_time[n_samples - 1] - my_time[0])/(double)(n_samples); // sampling time interval
	double fCut = 0.5/deltaT; // cut frequency = 0.5 * sampling frequency from WST
	int n_min = floor(tMin/deltaT) - range; // first sample to use
	int n_max = ceil(tMax/deltaT) + range; // last sample to use
	n_min = std::max(n_min,0); // check low limit
	n_max = std::min(n_max, (int)(n_samples)); // check high limit
	int n_0 = (n_min + n_max)/2;

	TComplex s(0.,0.); // Fourier transform at freq
	TComplex I(0.,1.); // i

  //std::cout << "deltaT: " << deltaT << n_min << " " << n_max << std::endl;
  //deltaT = 0.01;
	for(int n = n_min; n <= n_max; n++)
	{
    double time = double(n*deltaT);
    //double time = double(deltaT*(n-n_0));
    //s += deltaT*(double)my_channel[n]*TComplex::Exp(-I*(2.*TMath::Pi()*freq*(n-n_0)*deltaT));//maybe don't need n_0 here, I think it will just add a phase to the fourier transform
    s += deltaT*(double)my_channel[n]*TComplex::Exp(-I*(2.*TMath::Pi()*freq*time));
    //std::cout << "s: " << s << " " << time << std::endl;
	}
  return s.Rho();
};
