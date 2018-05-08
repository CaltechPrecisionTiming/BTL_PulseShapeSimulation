#include <iostream>
#include "TGraph.h"
#include "TFile.h"
#include "TH1F.h"

void time_resolution_no_dcr()
{
  double sigma_t[] = {15.4, 19.7, 21.3, 23.5};
  double n_threshold[] = {20., 40., 60., 80.};
  int n_points = 4;
  TGraph* g = new TGraph(n_points,n_threshold,sigma_t);
  TCanvas *cv = new TCanvas("Cv","Cv", 800,800);
  cv->SetLeftMargin(0.13);
  cv->SetBottomMargin(0.12);
  cv->SetRightMargin(0.05);
  //h->GetXaxis()->SetRangeUser(-1e5, 0);
  g->Draw("AC*");
  cv->SaveAs("deltaT_vs_npe.pdf");
  cv->SetLogy();
  cv->SaveAs("deltaT_vs_npe_log.pdf");
};
