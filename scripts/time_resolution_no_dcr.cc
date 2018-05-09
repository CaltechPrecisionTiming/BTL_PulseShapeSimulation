#include <iostream>
#include "TGraph.h"
#include "TFile.h"
#include "TH1F.h"
#include "TF1.h"

void time_resolution_no_dcr()
{
  double sigma_t20[] = {47.1, 27.1, 15.5, 9.3};
  double sigma_t40[] = {57.8, 32.3, 19.7, 11.5};
  double sigma_t60[] = {71.8,37.4,21.3,12.1};
  double sigma_t80[] = {78.8, 42.5, 23.5, 13.7};
  TF1* f_20 = new TF1("f_20", "sqrt(20.)*40000./x", 2000, 50000);
  TF1* f_40 = new TF1("f_40", "sqrt(40.)*40000./x", 2000, 50000);
  TF1* f_60 = new TF1("f_60", "sqrt(60.)*40000./x", 2000, 50000);
  TF1* f_80 = new TF1("f_80", "sqrt(80.)*40000./x", 2000, 50000);
  double Npe[] = {4000., 8000., 16000., 32000.};
  int n_points = 4;
  TGraph* g1 = new TGraph(n_points, Npe,sigma_t20);
  TGraph* g2 = new TGraph(n_points, Npe,sigma_t40);
  TGraph* g3 = new TGraph(n_points, Npe,sigma_t60);
  TGraph* g4 = new TGraph(n_points, Npe,sigma_t80);
  g1->SetMarkerColor(2);
  g1->SetMarkerSize(2);
  g1->SetMarkerStyle(20);

  g2->SetMarkerColor(4);
  g2->SetMarkerSize(2);
  g2->SetMarkerStyle(20);

  g3->SetMarkerColor(kPink-4);
  g3->SetMarkerSize(2);
  g3->SetMarkerStyle(20);

  g4->SetMarkerColor(kGreen);
  g4->SetMarkerSize(2);
  g4->SetMarkerStyle(20);
  TCanvas *cv = new TCanvas("Cv","Cv", 800,800);
  cv->SetLeftMargin(0.13);
  cv->SetBottomMargin(0.12);
  cv->SetRightMargin(0.05);
  //h->GetXaxis()->SetRangeUser(-1e5, 0);

  g4->Draw("AP");
  g2->Draw("P+same");
  g3->Draw("P+same");
  g1->Draw("P+same");
  f_20->SetLineStyle(2);
  f_20->Draw("same");
  f_40->SetLineStyle(2);
  f_40->SetLineColor(4);
  f_40->Draw("same");
  f_60->SetLineStyle(2);
  f_60->SetLineColor(kPink-4);
  f_60->Draw("same");
  f_80->SetLineStyle(2);
  f_80->SetLineColor(kGreen);
  f_80->Draw("same");


  cv->SaveAs("deltaT_vs_npe.pdf");
  cv->SetLogy();
  cv->SetLogx();
  cv->SaveAs("deltaT_vs_npe_log.pdf");
};
