#include <iostream>
#include <fstream>
#include "TTree.h"
#include "TFile.h"

void ConvertLGADSimulationToRoot(string inputFilenamePattern, string outputFilename ) {
  
  
  TTree* pulse = new TTree("pulse", "Digitized waveforms");

  const int npoints = 1500;
  const float stepSize = 0.001; //in ns
  float time[npoints];
  float channel[npoints];
  uint i_evt = 1;

  pulse->Branch("i_evt", &i_evt, "i_evt/i");
  pulse->Branch("channel", channel, Form("channel[%d]/F", npoints));
  pulse->Branch("time", time, Form("time[%d]/F", npoints));

  ifstream inputfile;
  inputfile.open(Form("%s_%d.txt",inputFilenamePattern.c_str(),i_evt));
  

  while (inputfile.is_open()) {
    cout << "Loading Event " << i_evt << "\n";

    //clear pulse shapes
    for (int i=0; i<npoints;i++) {
      channel[i] = 0;
    }


    string tmp;
    //skip first 2 lines, seems to be some header information
    std::getline(inputfile,tmp);
    std::getline(inputfile,tmp);

    int tmpIndex = 0;
    const double timeZero = 5e-9;
    while (!inputfile.eof() && tmpIndex < npoints) {
      double tmpAmp;
      double tmpTime;
      inputfile >> tmpTime >> tmpAmp;     
 
      //save points every 1ps      
      if ( int(std::round((tmpTime - timeZero)*1e13)) % 10 == 0 ) {
	//cout << tmpIndex << " : " << tmpTime << " " << tmpAmp << "\n";
	time[tmpIndex] = std::round((tmpTime - timeZero)*1e12);
	channel[tmpIndex] = -1*tmpAmp;
	tmpIndex++;
      }
    }

    pulse->Fill();

    //go to next event
    i_evt++;
    inputfile.close();    
    inputfile.open(Form("%s_%d.txt",inputFilenamePattern.c_str(),i_evt));   
  }


  TFile *outputfile = new TFile(outputFilename.c_str(),"RECREATE");
  pulse->Write();
  outputfile->Close();
  

}
