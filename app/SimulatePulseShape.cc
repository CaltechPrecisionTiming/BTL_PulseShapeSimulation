#include <iostream>
#include <string>
//ROOT
#include <TCanvas.h>
#include <TFile.h>
//LOCAL
#include <PulseShape.hh>

bool PulseShape::_info    = false;
bool PulseShape::_debug   = false;
bool PulseShape::_warning = true;

int main ( int argc, char** argv )
{
  const double pe_threshold = 20;
  const int n_experiments = 1000;
  PulseShape* ps;
  TGraph* total_pulse;
  TGraph* scintillation_pulse;
  TGraph* dark_noise;

  TH1F* h = new TH1F("pulse_time", "pulse_time", 2000, -10,10);

  double step = 0.01;
  double x_low  = -1e1;
  double x_high = 3e1;

  const int npoints  = (x_high-x_low)/step;
  std::cout << "[INFO] number of points per pulse: " << npoints << std::endl;
  std::cout << "[INFO] sampling rate is: " << step  << " ns" << std::endl;
  double x[npoints];
  double y[npoints], y_sc[npoints], y_dc[npoints];
  for ( int j = 0; j < n_experiments; j++ )
  {
    if ( j % 100 == 0 )std::cout << "experiment #" << j << std::endl;
    //reset variables and objects
    ps = new PulseShape("gauss");
    ps->SetNpe( 4.0e3 );
    ps->SetDCR( 30. );
    ps->SetSinglePhotonResponse( .1 );//units in ns
    ps->SetScintillationDecay( 40. );//units in ns
    for( int i = 0; i < npoints; i++ ) y[i] = x[i] = 0.0;
    double y_max = 0;
    for( int i = 0; i < npoints; i++ )
    {
      x[i] = x_low + double(i)*step;
      //if ( i % 1000 == 0 ) std::cout << "iteration #" << i << std::endl;
      //y[i]  = ps->Convolution(x[i], "Gauss", "RandomExp");
      y_sc[i]  = ps->ScintillationPulse(x[i]);
      y_dc[i]  = ps->DarkNoise(x[i], x_low, x_high);
      y[i]     = y_sc[i] + y_dc[i];
      if( y[i] > y_max ) y_max = y[i];
    }
    delete ps;//release memory of pulseshape object.

    double t_stamp = -999;
    for( int i = 0; i < npoints; i++ )
    {
      if ( y[i] > pe_threshold )
      {
        t_stamp = (x[i]+x[i+1])/2.;
        break;
      }
    }

    h->Fill(t_stamp);
  }

  ps = NULL;

  TFile* f = new TFile("out.root", "recreate");
  total_pulse = new TGraph( npoints, x, y);
  scintillation_pulse = new TGraph( npoints, x, y_sc);
  dark_noise = new TGraph( npoints, x, y_dc);
  TCanvas *cv = new TCanvas("Cv","Cv", 800,800);
  cv->SetLeftMargin(0.13);
  cv->SetBottomMargin(0.12);
  cv->SetRightMargin(0.05);
  //h->GetXaxis()->SetRangeUser(-1e5, 0);
  total_pulse->Draw("AC*");
  cv->SaveAs("Convolution1.pdf");
  cv->SetLogy();
  cv->SaveAs("Convolution2.pdf");
  total_pulse->Write("total_pulse");
  scintillation_pulse->Write("scintillation_pulse");
  dark_noise->Write("dark_noise");
  h->Write("h");
  f->Close();

  return 0;
};
