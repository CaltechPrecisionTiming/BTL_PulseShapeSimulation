#!/Users/cmorgoth/anaconda2/bin/python
import os
Npes = [4000, 5000, 6000, 8000, 10000]
nThresholds = [5, 10, 20, 40, 80, 100]
ScintillationDecayConstants = [40]#ns
ScintillationRisetimes = [0]
SinglePhotonRisetimeResponses = [1, 1.5, 3, 6]
SinglePhotonDecaytimeResponses = [2, 4, 6, 8, 15, 30, 60]
DarkCountRates = [0, 1, 10, 30, 50, 100]
HighPassFilterRCs = [0.1, 0.25, 0.5, 0.75, 1, 3, 6, 100000]


for Npe in Npes:
    for nThreshold in nThresholds:
        for ScintillationDecayConstant in ScintillationDecayConstants:
            for ScintillationRisetime in ScintillationRisetimes:
                for SinglePhotonRisetimeResponse in SinglePhotonRisetimeResponses:
                    for SinglePhotonDecaytimeResponse in SinglePhotonDecaytimeResponses:
                        for DarkCountRate in DarkCountRates:
                            for HighPassFilterRC in HighPassFilterRCs:
                                #print HighPassFilterRC
                                os.system("mkdir -p config/config_Npe_"+str(Npe))
                                filename = "config/config_Npe_"+str(Npe)+"/Npe_"+str(Npe)+"_nTh_"+str(nThreshold)+"_ScDecay_"+str(ScintillationDecayConstant)+"_ScRise_"+str(ScintillationRisetime)+"_SPRT_"+str(SinglePhotonRisetimeResponse)+"_SPDT_"+str(SinglePhotonDecaytimeResponse)+"_DCR_"+str(DarkCountRate)+"_HighPassFilterRC_"+str(HighPassFilterRC)+".cfg"
                                f = open(filename, 'w')
                                f.write("Npe "+str(Npe))
                                f.write("\nnThreshold "+str(nThreshold))
                                f.write("\nScintillationDecayConstant "+str(ScintillationDecayConstant))
                                f.write("\nScintillationRisetime "+str(ScintillationRisetime))
                                f.write("\nSinglePhotonRisetimeResponse "+str(SinglePhotonRisetimeResponse))
                                f.write("\nSinglePhotonDecaytimeResponse " + str(SinglePhotonDecaytimeResponse))
                                f.write("\nDarkCountRate "+str(DarkCountRate))
                                f.write("\nHighPassFilterRC "+str(HighPassFilterRC)+"\n")
                                f.close()
