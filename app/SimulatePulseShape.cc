#include <iostream>
#include <string>
//ROOT
#include <TCanvas.h>
//LOCAL
#include <PulseShape.hh>

int main ( int argc, char** argv )
{
  PulseShape* ps = new PulseShape("gauss");
  TGraph* gr;

  double step = 0.01;
  double x_low  = -1e1;
  double x_high = 1.5e2;

  const int npoints  = (x_high-x_low)/step;
  double x[npoints];
  double y[npoints];
  TH1F* h = new TH1F("h", "h", 2e3,-1e3, 1e3);
  for( int i = 0; i < npoints; i++ )
  {
    x[i] = x_low + double(i)*step;
    if ( i % 1000 == 0 ) std::cout << "iteration #" << i << std::endl;
    y[i] = ps->Convolution(x[i], "Gauss", "RandomExp");
    //h->SetBinContent(i+1, conv);
    //if ( conv != 0 )  std::cout << i+1 << " " <<  conv << std::endl;
  }

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

  return 0;
};
