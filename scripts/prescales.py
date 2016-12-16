import ROOT
from array import array
import os
from math import sqrt

# this was created when MT2DrawTools could be compiled standalone
ROOT.gSystem.Load("/mnt/t3nfs01/data01/shome/casal/ana743/src/MT2DrawTools_cc")
from ROOT import MT2DrawTools

lumi = 27.7

MT2DrawTools.setStyle()
lumilabel = MT2DrawTools.getLabelTop(lumi)


OVERWRITE = False # use cached files if they exist, if true recreate cached files

base_dir = "dcap://t3se01.psi.ch:22125//pnfs/psi.ch/cms/trivcat/store/user/casal/MT2production/80X/PostProcessed/prodNov08_allRuns_forQCD/"

# can contain subdirs, and wildcards
# if separated in chunks, they will be cached for posible quick access later
### re-reco up to runG 27.7/fb ###
files = [ "JetHT_Run2016B_23Sep2016_v3*.root", "JetHT_Run2016C_23Sep2016_v1*.root", "JetHT_Run2016D_23Sep2016_v1*.root", "JetHT_Run2016E_23Sep2016_v1*.root", "JetHT_Run2016F_23Sep2016_v1*.root", "JetHT_Run2016G_23Sep2016_v1*.root" ]

### full 2016 re-reco dataset ###
#files = [ "JetHT_Run2016*.root" ]


def getHistos(base, afile, ow): # will return [ht900, ht475, ht350, ht125]

    dirName = base.rstrip("/").split("/")[-1]
    if ( not os.path.isdir(dirName) ):
        os.mkdir(dirName)

    theFile = dirName + "/" + afile.replace("*","_a_")
    if ( not theFile.endswith(".root") ):
        theFile += ".root"
    subname = theFile.rstrip(".root").split("/")[-1]
    
    if ( not ow and os.path.isfile(theFile) ): #get histos from cached file
        print "reading from cached file", theFile
        rfile = ROOT.TFile(theFile, "READ")
        ROOT.gROOT.cd()
        return [rfile.Get("ht900"+"_"+subname).Clone(), 
                rfile.Get("ht475"+"_"+subname).Clone(), 
                rfile.Get("ht350"+"_"+subname).Clone(), 
                rfile.Get("ht125"+"_"+subname).Clone() ]

    else:

        print "creating histos from", base+"/"+afile

        t = ROOT.TChain("mt2")
        t.Add(base+"/"+afile)
   
        ht900 = ROOT.TH1F("ht900"+"_"+subname,"ht900",68,300,2000)
        ht475 = ROOT.TH1F("ht475"+"_"+subname,"ht475",68,300,2000)
        ht350 = ROOT.TH1F("ht350"+"_"+subname,"ht350",68,300,2000)
        ht125 = ROOT.TH1F("ht125"+"_"+subname,"ht125",68,300,2000)

        t.Draw("ht>>ht900"+"_"+subname,"isGolden && HLT_PFHT900"         , "goff");  ht900.Sumw2()
        t.Draw("ht>>ht475"+"_"+subname,"isGolden && HLT_PFHT475_Prescale", "goff");  ht475.Sumw2() 
        t.Draw("ht>>ht350"+"_"+subname,"isGolden && HLT_PFHT350_Prescale", "goff");  ht350.Sumw2() 
        t.Draw("ht>>ht125"+"_"+subname,"isGolden && HLT_PFHT125_Prescale", "goff");  ht125.Sumw2() 

        # cache histos in file
        rfile = ROOT.TFile(theFile, "RECREATE")
        rfile.cd()
        ht900.Write()
        ht475.Write()
        ht350.Write()
        ht125.Write()
        rfile.Close()

        return [ht900, ht475, ht350, ht125]
        

for i,afile in enumerate(files):
    [h1, h2, h3, h4] = getHistos(base_dir, afile, OVERWRITE)
    print "got histos for file", afile
    if i==0:
        [ht900, ht475, ht350, ht125] = [h1, h2, h3, h4]
    else:
        ht900.Add(h1)
        ht475.Add(h2)
        ht350.Add(h3)
        ht125.Add(h4)


r900over475 = ht900.Clone("r900over475");  r900over475.Divide(ht475)
r900over350 = ht900.Clone("r900over350");  r900over350.Divide(ht350)
r475over350 = ht475.Clone("r475over350");  r475over350.Divide(ht350)
r475over125 = ht475.Clone("r475over125");  r475over125.Divide(ht125)
# double binning for low stats 900over125
ht900.Rebin(2); ht125.Rebin(2)
r900over125 = ht900.Clone("r900over125");  r900over125.Divide(ht125)
r900over475.SetLineWidth(2); r900over475.SetLineColor(1)
r900over350.SetLineWidth(2); r900over350.SetLineColor(2)
r900over125.SetLineWidth(2); r900over125.SetLineColor(4)
r475over350.SetLineWidth(2); r475over350.SetLineColor(8)
r475over125.SetLineWidth(2); r475over125.SetLineColor(51)

r900over350.GetXaxis().SetTitle("H_{T} [GeV]")
r900over350.GetYaxis().SetTitle("ratio of events")
r900over350.GetXaxis().SetTitleSize(0.06); r900over350.GetYaxis().SetTitleSize(0.06); 
r900over350.GetYaxis().SetRangeUser(0.1,12000)
can = ROOT.TCanvas()
can.SetLogy()

r900over350.Draw("e")
r900over475.Draw("esame")
r900over125.Draw("esame")
r475over350.Draw("esame")
r475over125.Draw("esame")

f_ps1  = ROOT.TF1("ps1" , "pol0(0)", 1000, 2000); f_ps1 .SetLineColor(1) ; f_ps1 .SetLineWidth(2)
f_ps2  = ROOT.TF1("ps2" , "pol0(0)", 1000, 2000); f_ps2 .SetLineColor(2) ; f_ps2 .SetLineWidth(2)
f_ps3  = ROOT.TF1("ps3" , "pol0(0)", 1000, 2000); f_ps3 .SetLineColor(4) ; f_ps3 .SetLineWidth(2)
f_ps12 = ROOT.TF1("ps12", "pol0(0)", 575, 2000); f_ps12.SetLineColor(8) ; f_ps12.SetLineWidth(2); f_ps12.SetLineStyle(7)
f_ps13 = ROOT.TF1("ps13", "pol0(0)", 575, 2000); f_ps13.SetLineColor(51); f_ps13.SetLineWidth(2); f_ps13.SetLineStyle(7)

r900over475.Fit(f_ps1 ,"R")
r900over350.Fit(f_ps2 ,"R")
r900over125.Fit(f_ps3 ,"R")
r475over350.Fit(f_ps12,"R")
r475over125.Fit(f_ps13,"R")

ps1 , eps1  = f_ps1 .GetParameter(0), f_ps1 .GetParError(0)
ps2 , eps2  = f_ps2 .GetParameter(0), f_ps2 .GetParError(0)
ps3 , eps3  = f_ps3 .GetParameter(0), f_ps3 .GetParError(0)
ps12, eps12 = f_ps12.GetParameter(0), f_ps12.GetParError(0)
ps13, eps13 = f_ps13.GetParameter(0), f_ps13.GetParError(0)

ps2b, eps2b = ps1*ps12, sqrt( (eps1*ps12)**2 + (ps1*eps12)**2 )
ps3b, eps3b = ps1*ps13, sqrt( (eps1*ps13)**2 + (ps1*eps13)**2 )

leg = ROOT.TLegend(0.356, 0.18, 0.5, 0.386)
leg.AddEntry(f_ps1 , "PFHT900/475: %.1f #pm %.1f" % (ps1 , eps1 ), "L")
leg.AddEntry(f_ps2 , "PFHT900/350: %.1f #pm %.1f" % (ps2 , eps2 ), "L")
leg.AddEntry(f_ps3 , "PFHT900/125: %.0f #pm %.0f" % (ps3 , eps3 ), "L")
leg.AddEntry(f_ps12, "PFHT475/350: %.2f #pm %.2f #rightarrow ps350 = %.1f #pm %.1f" % (ps12, eps12, ps2b, eps2b), "L")
leg.AddEntry(f_ps13, "PFHT475/125: %.1f #pm %.1f #rightarrow ps125 = %.0f #pm %.0f" % (ps13, eps13, ps3b, eps3b), "L")
leg.SetTextSize(0.026)

#weighted averages 
def weightedAverage(data, errors):
    weights = [1/x/x for x in errors]
    mean = sum([x*y for x,y in zip(data,weights)])/sum(weights)
    std  = sqrt(1/sum(weights))
    return mean,std

data, error = [ps2,ps2b], [eps2, eps2b]
mean,std = weightedAverage(data, error)
print "averaged ps350 = %.1f +/- %.1f" % (mean, std)
data, error = [ps3,ps3b], [eps3, eps3b]
mean,std = weightedAverage(data, error)
print "averaged ps125 = %.0f +/- %.0f" % (mean, std)


leg.Draw("same")

lumilabel.Draw("same")

can.SaveAs(base_dir.rstrip("/").split("/")[-1]+"/prescales.pdf")
can.SaveAs(base_dir.rstrip("/").split("/")[-1]+"/prescales.png")
