#include <string>
#include <PulseShape.hh>

PulseShape::PulseShape()
{

};

PulseShape::PulseShape( std::string function_name )
{
  //this->function_name = function_name;
};

PulseShape::PulseShape( std::string function_name, std::string integration_method )
{

};

PulseShape::~PulseShape()
{

};

double PulseShape::Gauss( double x, double mean, double sigma, bool norm )
{
  return TMath::Gaus( x, mean, sigma, norm);
};

double PulseShape::Exp( double x, double exponent )
{
  if( x < 0 ) return 0.0;
  return exp(-1.0*exponent*x);
};

double PulseShape::RandomExp( double x, double exponent )
{
  TRandom3 r(0);
  return r.Poisson( 1e2*Exp( x, exponent ) );
};

double PulseShape::Convolution( double x, std::string function_name1, std::string function_name2 )
{
  double value = .0;
  double step_size = 0.5; //pico seconds units of the whole thing
  double x_low  = -1e3;
  double x_high = 1e3;
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
