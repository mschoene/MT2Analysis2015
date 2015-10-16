import ROOT
from array import array
from math import sqrt
import os


# before running compile first MT2DrawTools.cc:
# .L src/MT2DrawTools.cc+
# one needs (either in rootlogon or by command line) the following before compiling:
# gSystem->SetIncludePath("-I/swshare/cms/slc6_amd64_gcc481/lcg/roofit/5.34.09-cms11/include/ -I./interface/");
ROOT.gSystem.Load("../src/MT2DrawTools_cc")
from ROOT import MT2DrawTools


lumi = 0.454 #fb

doMC = True  # use pseudo-data or real data
plotMC = 'mc'  # plot a MC band ('mc': mc stats, 'lumi': lumi stats, otherwise: don't plot)
#relumi = lumi   # if doMC, then mc uncertainties are re-scaled for the given number
relumi = 3.0   # if doMC, then mc uncertainties are re-scaled for the given number

doLogx = True

MT2DrawTools.setStyle()
pandolfiniWantsThis = MT2DrawTools.getLabelTop(lumi if not doMC else relumi)

doDiffMetMht = True
doHFjetVeto = False

doRand = False

qcdbin = ""
#qcdbin = "_from300"  # high HT

regionSet = "zurich_onlyHT"
#regionSet = "zurich_HTtriggers"
#regionSet = "zurich_HTtriggers2"
#regionSet = "zurich_bMerged"
#regionSet = "zurich"

#samples = "PHYS14_v6"
#samples = "Spring15qcd"
samples = "74X_jecV4_MET30_QCD"
samples = "Spring15_25ns_qcdSkim"

rootfile = ROOT.TFile( "../analysis/QCDcontrolRegion_%s_%s_%.3ffb%s%s%s/qcdCRmod%s.root" % (samples,regionSet, lumi, "" if doDiffMetMht else "_noDiffMetMht", "" if not doHFjetVeto else "_HFjetVeto", qcdbin, "_rand" if doRand else "") )

#regions = ["HT450to575_j2to3_b0", "HT450to575_j4to6_b0", "HT575to1000_j2to3_b0", "HT575to1000_j4to6_b0"]
#regions = ["HT1000to1500_j2to3_b0", "HT1000to1500_j4to6_b0", "HT1500toInf_j4to6_b0", "HT1500toInf_j2to3_b0"]
#regions = ["HT450to575_j2to3_b0", "HT450to575_j4to6_b0" ]#, 
#regions = ["HT575to1000_j2to3_b0", "HT575to1000_j4to6_b0"]
regions = ["HT450to575_j2to3_b0", "HT450to575_j4to6_b0", "HT450to575_j7toInf_b0", "HT450to575_j2to3_b1", "HT450to575_j4to6_b1", "HT450to575_j7toInf_b1", "HT575to1000_j2to3_b0", "HT575to1000_j4to6_b0", "HT575to1000_j7toInf_b0", "HT575to1000_j2to3_b1", "HT575to1000_j4to6_b1", "HT575to1000_j7toInf_b1", "HT1000to1500_j2to3_b0", "HT1000to1500_j4to6_b0" , "HT1000to1500_j7toInf_b0", "HT1000to1500_j2to3_b1", "HT1000to1500_j4to6_b1" , "HT1000to1500_j7toInf_b1", "HT1000to1500_j2to3_b2", "HT1000to1500_j4to6_b2" , "HT1000to1500_j7toInf_b2", "HT1500toInf_j2to3_b0", "HT1500toInf_j4to6_b0", "HT1500toInf_j7toInf_b0", "HT1500toInf_j2to3_b1", "HT1500toInf_j4to6_b1", "HT1500toInf_j7toInf_b1", "HT1500toInf_j2to3_b2", "HT1500toInf_j4to6_b2", "HT1500toInf_j7toInf_b2" ]#, 

if "bMerged" in regionSet:
    regions = ["HT450to575_j2to3_b0toInf", "HT450to575_j4to6_b0toInf", "HT450to575_j7toInf_b0toInf", "HT575to1000_j2to3_b0toInf", "HT575to1000_j4to6_b0toInf", "HT575to1000_j7toInf_b0toInf", "HT1000to1500_j2to3_b0toInf", "HT1000to1500_j4to6_b0toInf", "HT1000to1500_j7toInf_b0toInf", "HT1500toInf_j2to3_b0toInf", "HT1500toInf_j4to6_b0toInf", "HT1500toInf_j7toInf_b0toInf"]

if "onlyHT" in regionSet:
    regions = ["HT450to575_j2toInf_b0toInf", "HT575to1000_j2toInf_b0toInf", "HT1000to1500_j2toInf_b0toInf", "HT1500toInf_j2toInf_b0toInf"]
if "HTtriggers" in regionSet:
    regions = ["HT450to575_j2toInf_b0toInf", "HT575to1000_j2toInf_b0toInf", "HT1000toInf_j2toInf_b0toInf"]
    #regions = ["HT1000toInf_j2toInf_b0toInf"]
if "HTtriggers2" in regionSet:
    regions = ["HT900toInf_j2toInf_b0toInf"]

qcdDir     = "qcdCR"     if plotMC=="mc" else "pdataSubCR"
dataDir    = "dataCR"    if not doMC else "pdataCR"
dataSubDir = "dataSubCR" if not doMC else "pdataSubCR"
#dataDir    = "pdataCR"     # pseudo-data
#dataSubDir = "pdataSubCR"  # pseudo-data

#prescale = {"": 1, "ps175": 175, "ps700": 700} #, "ps2400": 2400}
prescale = {"": 1, "ps60": 60, "ps180": 180} #, "ps2400": 2400}

h_qcd, h_data, h_dataSub = {}, {},{}
gr_pow   = {}
f_pow    = {}
band_pow = {}
gr_exp   = {}
f_exp    = {}
band_exp = {}

#fitW_lo, fitW_hi = 60, 100
fitW_lo, fitW_hi = 70, 100

y_lo, y_hi = 0.02, 20
xmax = 250 if not doLogx else 800

for region in regions:
    for ps in prescale.keys():
        if plotMC=="lumi" or ps=="":
            h_qcd[region+ps] = rootfile.Get("{name}/{reg}/ratio_{name}_{reg}".format( name=(qcdDir+ps), reg=region ) )
            h_qcd[region+ps].SetFillColor  (38);           h_qcd[region+ps].SetFillStyle(3001); 
            h_qcd[region+ps].SetMarkerColor(38);           h_qcd[region+ps].SetMarkerStyle(5); 
            h_qcd[region+ps].SetXTitle("M_{T2} [GeV]");    h_qcd[region+ps].SetYTitle("r_{#phi}(M_{T2})"); 
            h_qcd[region+ps].GetXaxis().SetTitleSize(0.05)
            h_qcd[region+ps].GetYaxis().SetTitleSize(0.05)
            #h_qcd[region+ps].GetXaxis().SetRangeUser(50,100)
            h_qcd[region+ps].GetXaxis().SetRangeUser(40,xmax)
     
        #if "HT450" in region and prescale[ps]<500:
        #    continue
        if "HT450" in region and prescale[ps]!=180:
            continue
        if "HT575" in region and prescale[ps]!=60:
            continue
        if "HT450" not in region and "HT575" not in region and ps != "":
            continue  # only plots prescales for low and medium HT regions
        h_data   [region+ps] = rootfile.Get("{name}/{reg}/ratio_{name}_{reg}".format( name=(dataDir)      , reg=region ) )
        h_dataSub[region+ps] = rootfile.Get("{name}/{reg}/ratio_{name}_{reg}".format( name=(dataSubDir+ps), reg=region ) )
        gr_exp   [region+ps] = rootfile.Get("{name}/{reg}/exp_{name}_{reg}"  .format( name=(dataSubDir+ps), reg=region ) )
        gr_pow   [region+ps] = rootfile.Get("{name}/{reg}/power_{name}_{reg}".format( name=(dataSubDir+ps), reg=region ) )

        if doMC and relumi!=lumi:  # re-scale uncertainties to given lumi
            for b in range(h_dataSub[region+ps].GetNbinsX()):
                h_dataSub[region+ps].SetBinError(b+1, h_dataSub[region+ps].GetBinError(b+1)*sqrt(lumi/relumi)) # yes, i did the math
                h_data   [region+ps].SetBinError(b+1, h_data   [region+ps].GetBinError(b+1)*sqrt(lumi/relumi)) # yes, i did the math
            

        # redoing fit so TVirtualFitter gets the info. need to figure out how to get the info from the already existing fit
        f_exp[region+ps] = gr_exp[region+ps].Clone("exp_fit_"+region+ps)
        f_exp[region+ps].SetLineColor(414)
        #h_dataSub[region+ps].Fit(f_exp[region+ps],"Q0","",50,100)
        h_dataSub[region+ps].Fit(f_exp[region+ps],"Q0","",fitW_lo,fitW_hi)
        band_exp[region+ps] = ROOT.TH1D("band_exp"+region+ps,"band",200,30,500)
        band_exp[region+ps].SetFillColor(30);   band_exp[region+ps].SetFillStyle(3001)
        band_exp[region+ps].SetLineColor(414);  band_exp[region+ps].SetLineWidth(2)
        (ROOT.TVirtualFitter.GetFitter()).GetConfidenceIntervals(band_exp[region+ps],0.68)

        f_pow[region+ps] = gr_pow[region+ps].Clone("pow_fit_"+region+ps)
        f_pow[region+ps].SetLineColor(ROOT.kYellow-3)
        #h_dataSub[region+ps].Fit(f_pow[region+ps],"Q0","",50,100)
        h_dataSub[region+ps].Fit(f_pow[region+ps],"Q0","",fitW_lo,fitW_hi)
        band_pow[region+ps] = ROOT.TH1D("band_pow"+region+ps,"band",200,30,500)
        band_pow[region+ps].SetFillColor(ROOT.kYellow-6);  band_pow[region+ps].SetFillStyle(3002)
        band_pow[region+ps].SetLineColor(ROOT.kYellow-3);  band_pow[region+ps].SetLineWidth(2)
        (ROOT.TVirtualFitter.GetFitter()).GetConfidenceIntervals(band_pow[region+ps],0.68)

        h_data   [region+ps].SetLineColor(50)          ; h_data   [region+ps].SetLineWidth(2);  
        h_dataSub[region+ps].SetLineColor(ROOT.kBlue-3); h_dataSub[region+ps].SetLineWidth(2); h_dataSub[region+ps].SetLineStyle(2); 
        h_data   [region+ps].SetMarkerStyle(20); h_data   [region+ps].SetMarkerSize(0.8); h_data   [region+ps].SetMarkerColor(50); 
        h_dataSub[region+ps].SetMarkerStyle(25); h_dataSub[region+ps].SetMarkerSize(0.8); h_dataSub[region+ps].SetMarkerColor(ROOT.kBlue-3); 
      


c = {}
l = {}
legg = {}

def getRegionText(x, text):
    text = text[text.find(x)+len(x):]
    x1 = text[0:text.find("to") if text.find("to")!=-1 else text.find("_") if text.find("_")!=-1 else len(text)]
    x2 = text[text.find("to")+2:text.find("_") if text.find("_")!=-1 else len(text)]
    if (text.find("to")==-1) or (text.find("_")!=-1 and text.find("_") <= text.find("to")):
        x2=x1
    return ( "%s = %s"%("N_{j}" if x=="j" else "N_{b}" if x=="b" else "H_{T}", x1 + (" GeV" if x=="HT" else "")) if x2==x1 else "%s #leq %s #leq %s"%(x1, "N_{j}" if x=="j" else "N_{b}" if x=="b" else "H_{T}", x2+ (" GeV" if x=="HT" else "")) ) if x2!="Inf" else ( "%s #geq %s"%("N_{j}" if x=="j" else "N_{b}" if x=="b" else "H_{T}", x1+ (" GeV" if x=="HT" else "")) )

for region in regions:
    for ps in prescale.keys():
        if "HT450" in region and prescale[ps]!=180:
            continue
        if "HT575" in region and prescale[ps]!=60:
            continue
        if "HT450" not in region and "HT575" not in region and ps != "":
            continue  # only plots prescales for low and medium HT regions
        #if "HT450" in region and prescale[ps]<500:
        #    continue
        #if "HT575" in region and (prescale[ps]>500 and prescale[ps]!=700):
        #    continue

        c[region+ps] = ROOT.TCanvas("can_"+region+ps,"can_"+region+ps)
        c[region+ps].SetLogy()
        if doLogx:
            c[region+ps].SetLogx()

        h_qcd    [region].GetYaxis().SetRangeUser(y_lo,y_hi);
        if plotMC=='mc':
            h_qcd[region].Draw("e2")
        elif plotMC=='lumi':
            h_qcd[region+ps].Draw("e2")
        else:
            h_qcd[region].Draw("axis")
        band_exp   [region+ps].Draw("e3same")
        band_pow   [region+ps].Draw("e3same")
        if plotMC=='mc':
            h_qcd[region].Draw("samee2")
        elif plotMC=='lumi':
            h_qcd[region+ps].Draw("samee2")
        f_exp    [region+ps].Draw("same")  #fit copy
        f_pow    [region+ps].Draw("same")  #fit copy
        #gr_exp   [region+ps].Draw("same") # original fit
        #gr_pow   [region+ps].Draw("same") # original fit
        h_data   [region+ps].Draw("esame")
        h_dataSub[region+ps].Draw("esame")
        text = "%s, %s, %s" % (getRegionText("HT", region),getRegionText("j", region),getRegionText("b", region))
        #text = "%s, %s, %s, L = 4 fb^{-1}%s" % (getRegionText("HT", region),getRegionText("j", region),getRegionText("b", region), "" if "HT1" in region else ", prescale = "+str(prescale[ps]))
        #text = "medium HT, HLT_PFHT475" if "HT575" in region else "low HT, HLT_PFHT350" if "HT450" in region else "high HT, HLT_PFHT900"
        #l[region+ps] = ROOT.TLatex(0.18,0.9, "%s, prescale = %d, rate = %.1f Hz" % (text, prescale[ps], (758. if "HT575" in region else 2800.)/prescale[ps]))
        l[region+ps] = ROOT.TLatex(0.52,0.88, text)
        l[region+ps].SetNDC()
        l[region+ps].SetTextSize(0.03); l[region+ps].SetTextFont(42)
        l[region+ps].Draw()
        pandolfiniWantsThis.Draw()

        legg[region+ps] = ROOT.TLegend(0.5,0.7,0.93,0.86)
        legg[region+ps].SetFillColor(0)
        legg[region+ps].SetBorderSize(0)
        legg[region+ps].AddEntry(h_data   [region+ps], "Data","pel")
        legg[region+ps].AddEntry(h_dataSub[region+ps], "Data (non-QCD bkgs subtracted)","pel")
        legg[region+ps].AddEntry(band_exp [region+ps], "Fit to subtracted data (exp)","fl")
        legg[region+ps].AddEntry(band_pow [region+ps], "Fit to subtracted data (pow)","fl")
        legg[region+ps].Draw()

        line1, line2 = ROOT.TLine(fitW_lo,y_lo,fitW_lo,y_hi), ROOT.TLine(fitW_hi,y_lo,fitW_hi,y_hi)
        line1.SetLineStyle(2);  line1.Draw()
        line2.SetLineStyle(2);  line2.Draw()
        arrow = ROOT.TArrow(fitW_lo,0.1,fitW_hi,0.1, 0.015, "<|>")
        arrow.Draw()
        txt = ROOT.TLatex(65,0.075, "fit window" ); txt.SetTextFont(42); txt.SetTextSize(0.029)
        txt.Draw()


        plotsDir = "plotsQCD_%s%s%s"% (regionSet, "" if doDiffMetMht else "_noDiffMetMht",qcdbin)
        if not os.path.isdir(plotsDir):
            os.mkdir(plotsDir)
        figure = "%s/ratio%s_%s_%.3fifb%s" % (plotsDir, "_rand" if doRand else "", region, lumi, "_"+ps if ps!="" else "")
        c[region+ps].SaveAs(figure + ("_mc" if doMC else "") +  ".eps")
        c[region+ps].SaveAs(figure + ("_mc" if doMC else "") +  ".png")
        c[region+ps].SaveAs(figure + ("_mc" if doMC else "") +  ".pdf")


#ps_low, ps_high = [2400, 700], [700, 175]
#ROOT.gStyle.SetEndErrorSize(6)
#can = {}
#latex = {}
#leg = {}
#for region in regions:
#    text = "medium HT, HLT_PFHT475" if "HT575" in region else "low HT, HLT_PFHT350" if "HT450" in region else "high HT, HLT_PFHT900"
#    latex[region] = ROOT.TLatex(0.18,0.9, text)
#    latex[region].SetNDC()
#    latex[region].SetTextSize(0.03); latex[region].SetTextFont(42)
#    leg[region] = ROOT.TLegend(0.165,0.165,0.65,0.35)
#    leg[region].SetFillColor(0)
#    leg[region].SetBorderSize(0)
#    for pskey, psval in prescale.iteritems():
#        ps = ps_low if "HT450" in region else ps_high
#        if psval not in ps:
#            continue
#        h_dataSub[region+pskey].SetLineColor(9 if psval==ps[0] else 1)
#        h_dataSub[region+pskey].SetLineStyle(1)
#        h_dataSub[region+pskey].SetLineWidth(1 if psval==ps[0] else 2)
#        h_dataSub[region+pskey].GetXaxis().SetRangeUser(50,100)
#        h_dataSub[region+pskey].SetXTitle("M_{T2}")
#        h_dataSub[region+pskey].SetYTitle("ratio")
#    can[region] = ROOT.TCanvas(region,region)
#    can[region].SetLogy()
#    h_dataSub[region+prescale.keys()[prescale.values().index(ps[0])]].Draw("e1")
#    h_dataSub[region+prescale.keys()[prescale.values().index(ps[1])]].Draw("e1same")
#    gr_exp   [region+prescale.keys()[prescale.values().index(ps[1])]].Draw("same")
#    leg[region].AddEntry(h_dataSub[region+prescale.keys()[prescale.values().index(ps[0])]], "prescale = %d, rate = %.1f Hz" % (ps[0], (758. if "HT575" in region else 2800.)/ps[0]), "l")
#    leg[region].AddEntry(h_dataSub[region+prescale.keys()[prescale.values().index(ps[1])]], "prescale = %d, rate = %.1f Hz" % (ps[1], (758. if "HT575" in region else 2800.)/ps[1]), "l")
#    latex[region].Draw()
#    leg  [region].Draw()
    


