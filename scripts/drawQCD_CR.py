import ROOT
from array import array


# before running compile first MT2DrawTools.cc:
# .L src/MT2DrawTools.cc+
# one needs (either in rootlogon or by command line) the following before compiling:
# gSystem->SetIncludePath("-I/swshare/cms/slc6_amd64_gcc481/lcg/roofit/5.34.09-cms11/include/ -I./interface/");
ROOT.gSystem.Load("../src/MT2DrawTools_cc")
from ROOT import MT2DrawTools

MT2DrawTools.setStyle()
#pandolfiniWantsThis = MT2DrawTools.getLabelTop()

rootfile = ROOT.TFile("../analysis/QCDcontrolRegion_PHYS14_qcdCR_zurich_4fb/pseudodata.root")

regions = ["HT450to575_j2to3_b0", "HT450to575_j4to6_b0", "HT575to1000_j2to3_b0", "HT575to1000_j4to6_b0", "HT1000to1500_j4to6_b0"]

qcdDir     = "qcdCR"
dataDir    = "dataCR"
dataSubDir = "dataSubCR"
prescale = {"": 1, "ps125": 125, "ps1200": 1200}

h_qcd, h_data, h_dataSub = {}, {},{}
gr_exp = {}

for region in regions:
    h_qcd[region] = rootfile.Get("{name}/{reg}/ratio_{name}_{reg}".format( name=qcdDir, reg=region ) )
    h_qcd[region].SetFillColor  (38);     h_qcd[region].SetFillStyle(3001); 
    h_qcd[region].SetMarkerColor(38);     h_qcd[region].SetMarkerStyle(5); 
    h_qcd[region].SetXTitle("M_{T2} [GeV]");    h_qcd[region].SetYTitle("ratio"); 
    h_qcd[region].GetXaxis().SetRangeUser(40,180)
    for ps in prescale.keys():
        h_data[region+ps]    = rootfile.Get("{name}/{reg}/ratio_{name}_{reg}".format( name=(dataDir+ps)   , reg=region ) )
        h_dataSub[region+ps] = rootfile.Get("{name}/{reg}/ratio_{name}_{reg}".format( name=(dataSubDir+ps), reg=region ) )
        gr_exp[region+ps]    = rootfile.Get("{name}/{reg}/exp_{name}_{reg}"  .format( name=(dataSubDir+ps), reg=region ) )

        h_data   [region+ps].SetLineColor(1); h_data   [region+ps].SetLineWidth(2); 
        h_dataSub[region+ps].SetLineColor(1); h_dataSub[region+ps].SetLineStyle(2); 


c = {}
l = {}

for region in regions:
    for ps in prescale.keys():

        c[region+ps] = ROOT.TCanvas()
        c[region+ps].SetLogy()

        h_qcd    [region].Draw("e2")
        h_data   [region+ps].Draw("esame")
        h_dataSub[region+ps].Draw("esame")
        gr_exp   [region+ps].Draw("same")
        l[region+ps] = ROOT.TLatex(0.25,0.9, "%s, prescale = %d" % (region, prescale[ps]))
        l[region+ps].SetNDC()
        l[region+ps].SetTextSize(0.03); l[region+ps].SetTextFont(42)
        l[region+ps].Draw()

        c[region+ps].SaveAs("plotsQCD/ratio_%s_%s.eps" % (region, ps))
        c[region+ps].SaveAs("plotsQCD/ratio_%s_%s.png" % (region, ps))
