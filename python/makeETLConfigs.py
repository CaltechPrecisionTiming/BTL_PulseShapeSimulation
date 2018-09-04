import sys
import argparse


if __name__ == '__main__':
    
    parser = argparse.ArgumentParser()
    #parser.add_argument("--tag", default="Razor2016_MoriondRereco",
    #        help="Analysis tag, e.g. Razor2015")
    #args = parser.parse_args()
    
    shapingTimeOptions = [0.5,1.0,1.5,2.0,2.5,3.0,4.0]
    SNROptions = [5,10,15,20,30,40,50,100,1000]
    LGADPulseLibraryOptions = ['',"LGADPulseLibrary_35micron_Gain15_Prerad.root","LGADPulseLibrary_55micron_Gain15_Prerad.root","LGADPulseLibrary_75micron_Gain15_Prerad.root","LGADPulseLibrary_W6_55micron_Gain16_5E14_480V.root","LGADPulseLibrary_W6_55micron_Gain10_1E15_550V.root","LGADPulseLibrary_W6_55micron_Gain4_3E15_700V.root"]

    shapingTimeName = {
	0.5 : "0p5",
        1.0 : "1p0",
        1.5 : "1p5",
        2.0 : "2p0",
        2.5 : "2p5",
        3.0 : "3p0",
        4.0 : "4p0"
        }  

    LGADPulseLibraryName = {
		"" : "DefaultFixedSignal",
		"LGADPulseLibrary_35micron_Gain15_Prerad.root" : "35MicronGain15Prerad",
		"LGADPulseLibrary_55micron_Gain15_Prerad.root" : "55MicronGain15Prerad",
		"LGADPulseLibrary_75micron_Gain15_Prerad.root" : "75MicronGain15Prerad",
		"LGADPulseLibrary_W6_55micron_Gain16_5E14_480V.root" : "W6_55MicronGain16_5E14_480V",
		"LGADPulseLibrary_W6_55micron_Gain10_1E15_550V.root" : "W6_55MicronGain10_1E15_550V",
		"LGADPulseLibrary_W6_55micron_Gain4_3E15_700V.root" : "W6_55MicronGain4_3E15_700V",
                }  
    LGADPulseLibrarySignalAmplitudeMean = {
        "" : 1.0,
        "LGADPulseLibrary_35micron_Gain15_Prerad.root" : 0.484e-3,
        "LGADPulseLibrary_55micron_Gain15_Prerad.root" :  0.564e-3,
        "LGADPulseLibrary_75micron_Gain15_Prerad.root" :  0.563e-3,
        "LGADPulseLibrary_W6_55micron_Gain16_5E14_480V.root" :  0.701e-3,
        "LGADPulseLibrary_W6_55micron_Gain10_1E15_550V.root" :  0.379e-3,
        "LGADPulseLibrary_W6_55micron_Gain4_3E15_700V.root" :  0.148e-3,
        }  
  

    for shapingTime in shapingTimeOptions:
        for SNR in SNROptions:
            for LGADPulseFile in LGADPulseLibraryOptions:
                filename = "config/config_ShapingTime"+shapingTimeName[shapingTime]+"_SNR"+str(SNR)+"_"+LGADPulseLibraryName[LGADPulseFile]+".cfg"
                print filename
                f = open(filename, 'w')
                f.write("NFilter 2 \n")
                f.write("ShapingTime "+str(shapingTime)+"\n")
                f.write("SNR "+str(SNR)+"\n")
                f.write("SignalAmplitudeMean "+str(LGADPulseLibrarySignalAmplitudeMean[LGADPulseFile])+"\n")
                f.write("randomSeed 1\n")
                if (LGADPulseFile != ""):
                    f.write("LGADSignalFilename /afs/cern.ch/work/s/sixie/public/releases/run2/Timing/CMSSW_9_0_2/src/BTL_PulseShapeSimulation/LGADPulseLibrary/" + LGADPulseFile + "\n")
                f.close()





  
