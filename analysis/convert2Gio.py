import ROOT
import array
import math

## NEEDS TO BE FIXED
# to correctly account for all different error correlation of all contributions
# across sources, across HT,NJ,NB regions, and across MT2 bins


#inFile = ROOT.TFile("EventYields_data_Run2015D_25nsGolden_qcd_2p3ifb/qcdControlRegion/macro/qcdEstimateData.root","READ")
#inFile = ROOT.TFile("EventYields_data_Run2016B_Golden/qcdControlRegion/macro/qcdEstimateData.root","READ")
inFile = ROOT.TFile("EventYields_data_Run2016_13ifb/qcdControlRegion//macro/nonQCDunc20/qcdEstimateData.root","READ")

inFile.cd("qcdEstimate")

regions = [subdir.GetName() for subdir in ROOT.gDirectory.GetListOfKeys()]

bins = [200, 300, 400, 600, 800, 1000, 1500]
nBins = len(bins)-1
h_macroRegions = {}
for x in range(1,3,1):
    for y in range(1,3,1):
        for z in range(1,4,1):
            if z==3 and (x==2 or y==2): continue
            macroRegion = "srMacroHT%dNJ%dNB%d" % (x,y,z)
            h_macroRegions[macroRegion] = ROOT.TH1F("h_mt2rebin_"+macroRegion,"",nBins,array.array('d',bins) )


def INT(str):
    return -1 if str == 'Inf' else int(str)

def addhistos(h1, h2): # add errors linearly
    for b in range(1,h1.GetNbinsX()+1):
        h1.SetBinContent(b, h1.GetBinContent(b)+h2.GetBinContent(b))
        h1.SetBinError  (b, h1.GetBinError  (b)+h2.GetBinError  (b))


for region in regions:
   HTregion = region.split("_")[0]
   NJregion = region.split("_")[1]
   NBregion = region.split("_")[2]

   
   HTmin, HTmax = INT(HTregion[2:HTregion.find("to")]), INT(HTregion[HTregion.find("to")+2:len(HTregion)])
   NJmin, NJmax = INT(NJregion[1:NJregion.find("to")]), INT(NJregion[NJregion.find("to")+2:len(NJregion)])
   if NBregion.find("to")==-1 : NBregion += 'to'+NBregion[-1]
   NBmin, NBmax = INT(NBregion[1:NBregion.find("to")]), INT(NBregion[NBregion.find("to")+2:len(NBregion)])
   
   this_hist = inFile.Get("qcdEstimate/"+region+"/yield_qcdEstimate_"+region)

   iHT = 1 if HTmin<1000 else 2
   iNJ = 1 if NJmin <4   else 2
   iNB = 1 if NBmin==0   else 2 if NBmin<3 else 3

   if iNB==3: iHT, iNJ = 1, 1

   # 3b's are going to srMacroHT1NJ1NB3 for all ht and all nj
   macroRegion = "srMacroHT%dNJ%dNB%d" % (iHT,iNJ,iNB)

   #if NBmin==3 and NJmin==2:
   #    continue              # skip 2-6j, 3b, as agreed with Gio

   print region, macroRegion

   if macroRegion == "srMacroHT2NJ2NB2":
       b = h_macroRegions[macroRegion].FindBin(250)
       print h_macroRegions[macroRegion].GetBinContent(b), h_macroRegions[macroRegion].GetBinError(b),this_hist.GetBinContent(b),this_hist.GetBinError(b)

   #h_macroRegions[macroRegion].Add(h_macroRegions[macroRegion],this_hist)
   addhistos(h_macroRegions[macroRegion],this_hist) # add errors linearly, more correct due to strong correlation
   if macroRegion == "srMacroHT2NJ2NB2":
       print h_macroRegions[macroRegion].GetBinContent(b), h_macroRegions[macroRegion].GetBinError(b)


outFile = ROOT.TFile("MacroQcdEstimateData.root","RECREATE")

hh_macroRegions = {} # final binning

def ReBin(h1,h2):
    for b in range(1,h1.GetNbinsX()+1):
        x1,y1,e1 = h1.GetBinCenter(b),h1.GetBinContent(b),h1.GetBinError(b)

        b2 = h2.FindBin(x1)
        if b2 == h2.GetNbinsX()+1: b2 -= 1 # add overflow to last bin

        y2,e2 = h2.GetBinContent(b2),h2.GetBinError(b2)
        y2 += y1
        #e2 = math.sqrt(e2*e2+e1*e1)
        e2 = math.sqrt(e2*e2+e1*e1) # linearly is more correct given strong correlation

        h2.SetBinContent(b2,y2)
        h2.SetBinError  (b2,e2)


for mr in h_macroRegions.keys():
    
    dir_mr = outFile.mkdir(mr)
    dir_mr.cd()

    if mr == "srMacroHT2NJ1NB1" or  mr == "srMacroHT2NJ2NB1" :
        bins = [200, 300, 400, 600, 800, 1000, 1500]
        hh_macroRegions[mr] = ROOT.TH1F("h_mt2rebin","",len(bins)-1,array.array('d',bins) )
    elif mr == "srMacroHT1NJ1NB1" or mr == "srMacroHT1NJ1NB2" or mr == "srMacroHT1NJ2NB1" or mr == "srMacroHT2NJ1NB2" or mr == "srMacroHT2NJ2NB2":
        bins = [200, 300, 400, 600, 800, 1000]
        hh_macroRegions[mr] = ROOT.TH1F("h_mt2rebin","",len(bins)-1,array.array('d',bins) )
    elif mr == "srMacroHT1NJ2NB2" or mr == "srMacroHT1NJ1NB3":
        bins = [200, 300, 400, 600, 1000]
        hh_macroRegions[mr] = ROOT.TH1F("h_mt2rebin","",len(bins)-1,array.array('d',bins) )
    else:
        print "macro region not defined ?", mr

    ReBin(h_macroRegions[mr],hh_macroRegions[mr])

outFile.Write()
