#!/Users/cmorgoth/anaconda2/bin/python

Npes = [2000, 4000, 6000, 8000, 10000, 20000]
nThresholds = [5, 10, 20, 30, 40, 60, 80, 100]
ScintillationDecayConstants = [40000]#ns
ScintillationRisetimes = [0]
SinglePhotonRisetimeResponses = [1, 1.5, 3,6]
SinglePhotonDecaytimeResponses = [2, 4, 6, 8, 15, 30, 60]
DarkCountRates = [0, 0.1, 0.5, 1, 5, 10, 30, 50, 100]
HighPassFilterRCs = [0.1, 0.25, 0.5, 0.75, 1., 3, 6, 100000]


for Npe in Npes:
    for nThreshold in nThresholds:
        for ScintillationDecayConstant in ScintillationDecayConstants:
            for ScintillationRisetime in ScintillationRisetimes:
                for SinglePhotonRisetimeResponse in SinglePhotonRisetimeResponses:
                    for SinglePhotonDecaytimeResponse in SinglePhotonDecaytimeResponses:
                        for DarkCountRate in DarkCountRates:
                            for HighPassFilterRC in HighPassFilterRCs:
                                #print HighPassFilterRC
                                filename = "Npe_"+str(Npe)+"_nTh_"+str(nThreshold)+"_ScDecay_"+str(ScintillationDecayConstant)+"_ScRise_"+str(ScintillationRisetime)+"_SPRT_"+str(SinglePhotonRisetimeResponse)+"_SPDT_"+str(SinglePhotonDecaytimeResponse)+"_DCR_"+str(DarkCountRate)+"_HighPassFilterRC_"+str(HighPassFilterRC)+".cfg"
                                print filename
