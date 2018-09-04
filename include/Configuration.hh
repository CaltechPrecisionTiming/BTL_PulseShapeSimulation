#ifndef Configuration_HH
#define Configuration_HH

// Std include
#include <iostream>
#include <fstream>
#include <sstream>
#include <map>

#include "TString.h"
#include "TRegexp.h"
using namespace std;

// This class reads a configuration file and extracts polarity, attenuation,
// and amplification factors for each digitizer channel.
//
// Configuration file lines should be of the following form (lines beginning with a '#' are ignored):
// baseline i_start n_samples
// CH  POLARITY  AMPLIFICATION  ATTENUATION (dB)  ALGORITHM  FILTER_WIDTH
// with:
// CH - channel number (integer)
// POLARITY - sign of the pulse ('+' or '-'). Pulses are supposed to have the peak below the baseline
// AMPLIFICATION - in dB, amount of amplification that was applied to the input (float)
// ATTENUATION - in dB, amount of attenuation that was applied to the input (float)
// ALGORITHM - indicates the algorithm to run to extract pulse times (string).Usa a '+' to separate them
//      G: gaussian fit, G##
//      Re: linear fit at the rising edge, Re##-## possible in %
//      LP#: Local polinomial fit. For the moment min 1, max 3. e.g LP2
//      else --> No action
// FILTER_WIDTH - gaussian kernel width for Weierstrass transform (gaussian filter).
//      If 0, no Weierstrass transform will be applied.

class Configuration {
    public:
      struct Channel {
        unsigned int N = 0;

        int polarity = +1;
        unsigned int counter_auto_pol_switch = 0;

        float amplification = 0;
        float attenuation = 0;
        TString algorithm = "";
        float gaus_fraction = 0.4;
        float re_bounds[2] = {0.15, 0.75};
        vector<int> PL_deg;

        float weierstrass_filter_width = 0;
      };

      // Read and initilize from the configuration file
      Configuration(){};
      Configuration(bool verb){ verbose = verb;};
      Configuration(std::string config_file, bool verb);
      void ParseConfigurationFile(std::string config_file);

      TString ParseCommandLine( int argc, char** argv, TString opt );
      void GetCommandLineArgs(int argc, char **argv);
      //command line arguments
      std::string output_file;
      std::string config_file;
      int n_experiments;

      map<unsigned int,Channel> channels;

      bool isValid();
      bool verbose = false;
      bool _warning = true;
      // get overall multiplier including polarity, amplification, and attenuation
      float getChannelMultiplicationFactor(unsigned int ch);

      // returns true if the specified channel is present in the config
      bool hasChannel(unsigned int ch);

      unsigned int baseline[2] = {20, 150};
      // must be between 0.1 and 0.9
      vector<float> constant_fraction = {0.15, 0.3, 0.45};
      vector<float> constant_threshold = {};

      vector<float> z_DUT = {-50., 50.};
      //simulation specifics
      int NFilter = 0; 
      double ShapingTime = 4.0;
      double SNR = 30.0;
      double SignalAmplitudeMean = 0;
      int randomSeed = 0;
      string LGADSignalFilename = "";

    private:
      // process one line of the config file
      void parseConfigurationLine(std::string line);
      // get next non-space token from config
      int nextConfigurationElement(std::stringstream &ss, std::string &item);
};

#endif
