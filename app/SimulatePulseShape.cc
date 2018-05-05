#include <iostream>
#include <string>
//ROOT
#include <TCanvas.h>
#include <TFile.h>
//LOCAL
#include <PulseShape.hh>

int main ( int argc, char** argv )
{
  PulseShape* ps = new PulseShape("gauss");
  TGraph* gr;
  TH1F* h = new TH1F("pulse_time", "pulse_time", 100, 8,12);

  double step = 0.01;
  double x_low  = -1e1;
  double x_high = 3e1;

  const int npoints  = (x_high-x_low)/step;
  std::cout << "[INFO] number of points per pulse: " << npoints << std::endl;
  std::cout << "[INFO] sampling rate is: " << step  << " ns" << std::endl;
  double x[npoints];
  double y[npoints];
  for ( int j = 0; j < 100; j++ )
  {
    for( int i = 0; i < npoints; i++ ) y[i] = x[i] = 0.0;
    double y_max = 0;
    for( int i = 0; i < npoints; i++ )
    {
      x[i] = x_low + double(i)*step;
      if ( i % 1000 == 0 ) std::cout << "iteration #" << i << std::endl;
      y[i]  = ps->Convolution(x[i], "Gauss", "RandomExp");
      if( y[i] > y_max ) y_max = y[i];
    }

    double t_stamp = -999;
    for( int i = 0; i < npoints; i++ )
    {
      if ( y[i] > 20.e-3 )
      {
        t_stamp = (x[i]+x[i+1])/2.;
        break;
      }
    }

    h->Fill(t_stamp);
  }
  TFile* f = new TFile("out.root", "recreate");
  gr = new TGraph( npoints, x, y);
  TCanvas *cv = new TCanvas("Cv","Cv", 800,800);
  cv->SetLeftMargin(0.13);
  cv->SetBottomMargin(0.12);
  cv->SetRightMargin(0.05);
  //h->GetXaxis()->SetRangeUser(-1e5, 0);
  gr->Draw("AC*");
  cv->SaveAs("Convolution1.pdf");
  cv->SetLogy();
  cv->SaveAs("Convolution2.pdf");
  gr->Write("graph_ps");
  h->Write("h");
  f->Close();

  return 0;
};
