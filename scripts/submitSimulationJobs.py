import os,sys
import argparse


if __name__ == '__main__':
    
    parser = argparse.ArgumentParser()
    parser.add_argument("--outputDir", default="/eos/cms/store/user/sixie/Timing/ETLSimulation/",
                        help="Output Directory for storing the ntuples")
    args = parser.parse_args()    
    outputDir = args.outputDir

    shapingTimeOptions = [
        #0.5,
        #1.0,
        #1.5,
        2.0,
        #2.5,
        #3.0,
        4.0
        ]
    SNROptions = [
        #5,
        10,15,
        20,
        #30,40,50,100,
        #1000
        ]
    LGADPulseLibraryOptions = [#'',
                               #"LGADPulseLibrary_35micron_Gain15_Prerad.root",
                               "LGADPulseLibrary_55micron_Gain15_Prerad.root",
                               #"LGADPulseLibrary_75micron_Gain15_Prerad.root",
                               #"LGADPulseLibrary_W6_55micron_Gain16_5E14_480V.root",
                               #"LGADPulseLibrary_W6_55micron_Gain10_1E15_550V.root",
                               #"LGADPulseLibrary_W6_55micron_Gain4_3E15_700V.root"
                               ]

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
    
    Jobs = [
            1,
            2,
            3,
            4,
            5,
             6,
             7,
             8,
             9,
             10,
            # 11,
            # 12,
            # 13,
            # 14,
            # 15,
            # 16,
            # 17,
            # 18,
            # 19,
            # 20,
            # 21,
            # 22,
            # 23,
            # 24,
            # 25,                 
            ]
    NEvents = 100

    for shapingTime in shapingTimeOptions:
        for SNR in SNROptions:
            for LGADPulseFile in LGADPulseLibraryOptions:
                configFile = "/afs/cern.ch/work/s/sixie/public/releases/run2/Timing/CMSSW_9_0_2/src/BTL_PulseShapeSimulation/config/config_ShapingTime"+shapingTimeName[shapingTime]+"_SNR"+str(SNR)+"_"+LGADPulseLibraryName[LGADPulseFile]+".cfg"
                #loop over jobs
                for jobnumber in Jobs:
                    
                    label = "ShapingTime"+shapingTimeName[shapingTime]+"_SNR"+str(SNR)+"_"+LGADPulseLibraryName[LGADPulseFile]
                    outputFilename = "ETLSimulation_"+label+"_Job"+str(jobnumber)+".root"
                    #print outputFilename
            
                    command = "bsub -q 1nd -o /afs/cern.ch/user/s/sixie/work/private/condor/res/ETLSimulation/ETLSimulation_"+label+"_Job"+str(jobnumber)+".out"+ " -J ETLSimulation_"+label+"_Job"+str(jobnumber)+" /afs/cern.ch/work/s/sixie/public/releases/run2/Timing/CMSSW_9_0_2/src/BTL_PulseShapeSimulation/scripts/runJob.csh "+configFile+" "+str(NEvents)+ " " + outputFilename + " " + outputDir
                    print command
                    os.system(command)



  
