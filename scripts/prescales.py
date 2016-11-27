import ROOT
from array import array
import os
from math import sqrt

ROOT.gSystem.Load("/mnt/t3nfs01/data01/shome/casal/ana743/src/MT2DrawTools_cc")
from ROOT import MT2DrawTools

lumi = 24.5

MT2DrawTools.setStyle()
lumilabel = MT2DrawTools.getLabelTop(lumi)

# 211/pb
#f = ROOT.TFile.Open("dcap://t3se01.psi.ch:22125//pnfs/psi.ch/cms/trivcat/store/user/pandolf/MT2production/80X/PostProcessed/data2016_v3/JetHT_Run2016B_PromptReco_v2_post.root")
# 589/pb
#f = ROOT.TFile.Open("dcap://t3se01.psi.ch:22125//pnfs/psi.ch/cms/trivcat/store/user/mangano/babies/80X/MT2/PostProcessed/data27May_v1/JetHT_Run2016B_PromptReco_v2_po0st.root")
# 804/pb
#f = ROOT.TFile.Open("dcap://t3se01.psi.ch:22125//pnfs/psi.ch/cms/trivcat/store/user/mangano/babies/80X/MT2/PostProcessed/data07June_JetHT_v1/JetHT_Run2016B_PromptReco_v2_post.root")

#t = f.Get("mt2")


t = ROOT.TChain("mt2")

# 2.07/fb
#t.Add("dcap://t3se01.psi.ch:22125//pnfs/psi.ch/cms/trivcat/store/user/mangano/MT2production/80X/PostProcessed/data13June_271036-274421_newCFG_v1/JetHT_Run2016B*.root")

# 2.07 + 1.89 /fb
#t.Add("dcap://t3se01.psi.ch:22125//pnfs/psi.ch/cms/trivcat/store/user/mangano/MT2production/80X/PostProcessed/data13June_271036-274421_newCFG_v1/JetHT_Run2016B*.root")
#t.Add("dcap://t3se01.psi.ch:22125//pnfs/psi.ch/cms/trivcat/store/user/mangano/MT2production/80X/PostProcessed/prodJune21_274422-275125_forQCD_v1_v2/JetHT_Run2016B*.root")

# all runB 5.4/fb
#t.Add("dcap://t3se01.psi.ch:22125//pnfs/psi.ch/cms/trivcat/store/user/mangano/MT2production/80X/PostProcessed/prodJune30_allRunB_forQCD_v1//JetHT_Run2016B*.root")

# JECv6, eta2.4,  5.89 (B) + 1.76 (C) /fb
#t.Add("dcap://t3se01.psi.ch:22125//pnfs/psi.ch/cms/trivcat/store/user/mangano/MT2production/80X/PostProcessed/prodJuly08_runB_forQCD_v1//JetHT_Run2016B*.root")
#t.Add("dcap://t3se01.psi.ch:22125//pnfs/psi.ch/cms/trivcat/store/user/mangano/MT2production/80X/PostProcessed/prodJuly08_runC_forQCD_v1//JetHT_Run2016C*.root")

# 5.94 (B) + 2.646 (C) + 0.649 (D) = 9.24/fb
#t.Add("dcap://t3se01.psi.ch:22125//pnfs/psi.ch/cms/trivcat/store/user/mangano/MT2production/80X/PostProcessed/prodJuly08_runB_forQCD_v1/JetHT_Run2016B*.root")
#t.Add("dcap://t3se01.psi.ch:22125//pnfs/psi.ch/cms/trivcat/store/user/mangano/MT2production/80X/PostProcessed/prodJuly15_runC_all_forQCD_v1//JetHT_Run2016C*.root")
#t.Add("dcap://t3se01.psi.ch:22125//pnfs/psi.ch/cms/trivcat/store/user/mangano/MT2production/80X/PostProcessed/prodJuly15_runD_276311-276384_forQCD_v1/JetHT_Run2016D*.root")

# 5.940 (B) + 2.646 (C) + 4.330 (D) = 12.92/fb
#t.Add("dcap://t3se01.psi.ch:22125//pnfs/psi.ch/cms/trivcat/store/user/mangano/MT2production/80X/PostProcessed/prodJuly08_runB_forQCD_v1/JetHT_Run2016B*.root")
#t.Add("dcap://t3se01.psi.ch:22125//pnfs/psi.ch/cms/trivcat/store/user/mangano/MT2production/80X/PostProcessed/prodJuly15_runC_all_forQCD_v1//JetHT_Run2016C*.root")
#t.Add("dcap://t3se01.psi.ch:22125//pnfs/psi.ch/cms/trivcat/store/user/mangano/MT2production/80X/PostProcessed/prodJuly19_runD_276311-276811_forQCD_v1/JetHT_Run2016D*.root")


# 24.5/fb run G up to 279931
t.Add("dcap://t3se01.psi.ch:22125//pnfs/psi.ch/cms/trivcat/store/user/casal/MT2production/80X/PostProcessed/prodSept13_all_forQCD/JetHT_Run2016*.root")


ht900 = ROOT.TH1F("ht900","ht900",68,300,2000)
ht475 = ROOT.TH1F("ht475","ht475",68,300,2000)
ht350 = ROOT.TH1F("ht350","ht350",68,300,2000)
ht125 = ROOT.TH1F("ht125","ht125",68,300,2000)

t.Draw("ht>>ht900","isGolden && HLT_PFHT900"         , "q");  ht900.Sumw2()
t.Draw("ht>>ht475","isGolden && HLT_PFHT475_Prescale", "q");  ht475.Sumw2() 
t.Draw("ht>>ht350","isGolden && HLT_PFHT350_Prescale", "q");  ht350.Sumw2() 
t.Draw("ht>>ht125","isGolden && HLT_PFHT125_Prescale", "q");  ht125.Sumw2() 
#t.Draw("ht>>ht125","isGolden && (HLT_PFHT125_Prescale||HLT_PFHT200_Prescale)", "q");  ht125.Sumw2() 


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

#can.SaveAs("prescales.pdf")
#can.SaveAs("prescales.png")
