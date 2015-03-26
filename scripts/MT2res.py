import ROOT
from array import array

# what can you configure?
# 1. samples to run on (setup works with up to 4 samples)
fsig = ["/scratch/casal/PHYS14_Production_Feb20/SMS_T1qqqq_2J_mGl1400_mLSP100_post.root", "/scratch/casal/PHYS14_Production_Feb20/SMS_T1qqqq_2J_mGl1000_mLSP800_post.root", "/scratch/casal/PHYS14_Production_Feb20/SMS_T1tttt_2J_mGl1500_mLSP100_post.root" ]#,"/scratch/casal/PHYS14_Production_Feb20/SMS_T1tttt_2J_mGl1200_mLSP800_post.root"]

# 2. resolution as a function of (so far binning only adapted to ht and mt2*):
vs = "ht"
#vs = "mt2"
#vs = "mt2_gen"

# before running compile first MT2DrawTools.cc:
# .L src/MT2DrawTools.cc+
# one needs (either in rootlogon or by command line) the following before compiling:
# gSystem->SetIncludePath("-I/swshare/cms/slc6_amd64_gcc481/lcg/roofit/5.34.09-cms11/include/ -I./interface/");
ROOT.gSystem.Load("../src/MT2DrawTools_cc")
from ROOT import MT2DrawTools

MT2DrawTools.setStyle()
pandolfiniWantsThis = MT2DrawTools.getLabelTop()

def graphStyle(graph, title, xtitle, ytitle, style, color):
    graph.SetTitle(title)
    graph.GetXaxis().SetTitle(xtitle)
    graph.GetYaxis().SetTitle(ytitle)
    #graph.GetXaxis().SetTitleSize(0.04)
    #graph.GetYaxis().SetTitleSize(0.04)
    #graph.GetYaxis().SetTitleOffset(1.2)
    graph.SetMinimum(0)
    graph.SetMaximum(100)
    graph.SetMarkerStyle(style)
    graph.SetMarkerSize(0.6)
    graph.SetMarkerColor(color)
    graph.SetLineColor(color)

tfile = {}
ttree = {}
signals  = []
signals2 = []
hist2d = {}
x,dx = {},{}
rms,drms = {},{}
sig,dsig = {},{}
histProj = {}
graphRMS, graphSig = {}, {}

colorlist = [46, 9, 8, 41]
colors = {}

for n,f in enumerate(fsig):
    for s in f.split("_"):
        if s[0]=="T":
            model=s
        if s[0:3] == "mGl":
            m1=s[3:]
        if s[0:5] == "mStop":
            m1=s[5:]
        if s[0:4] == "mLSP":
            m2=s[4:]
    signal = model+"_"+m1+"_"+m2
    print signal
    signals .append(signal)
    signals2.append("%s(%s,%s)"%(model,m1,m2))

    colors[signal] = colorlist[n]

    tfile[signal] = ROOT.TFile(f)
    ttree[signal] = tfile[signal].Get("mt2")
    
    if "mt2" in vs:
        hist2d[signal] = ROOT.TH2F("hist2d_"+signal,"hist2d",52,200,1500,400,-1000,1000) # 25 GeV in genMT2, 5 GeV in MT2res
    else:
        hist2d[signal] = ROOT.TH2F("hist2d_"+signal,"hist2d",160,450,4450,400,-1000,1000) # 25 GeV in HT, 5 GeV in MT2res
    
    noGenLep = True
    
    sel = "&& ngenLep==0 && ngenTau==0" if noGenLep else ""
    sel += " && mt2>200 && ht > 450 && nJet40>=2"
    ttree[signal].Draw("mt2-mt2_gen:%s>>hist2d_%s" % (vs,signal),"mt2_gen>-1"+sel,"colz")
    
    if signal == "T1qqqq_1000_800":
        bins = [450,500,575, 650, 750, 850, 1000, 1200, 1500, 2000]
    elif signal == "T1tttt_1200_800":
        bins = [450,550 ,650, 750, 1000, 1200, 1500, 2000]
    elif signal == "T1qqqq_1400_100":
        bins = [450,700, 900, 1050] + range(1200,2300,100) + [2350, 2500, 2700, 3000, 3500]
    elif signal == "T1tttt_1500_100":
        bins = [450, 1000] + range(1200,2300,100) + [2350, 2500, 2700, 3000, 3500]

    if "mt2" in vs:
        if signal == "T1qqqq_1000_800" or signal == "T1tttt_1200_800":
            bins = range(200,451,25)+[500, 550, 650, 800]
        if signal == "T1qqqq_1400_100" or signal == "T1tttt_1500_100":
            bins = range(200,401,25)+range(450,501,50)+range(600,1100,100)
    
    x[signal],dx[signal] = [],[]
    rms[signal], drms[signal] = [],[]
    sig[signal], dsig[signal] = [],[]
    histProj[signal] = []
    canfit = ROOT.TCanvas ()
    for i,ht in enumerate(bins):
        htname = vs+str(bins[i])+"to"+ (str(bins[i+1]) if i<len(bins)-1 else "Inf")
        htname+= "_"+signal
        htcut  = ("%d < X < %d" % (bins[i],bins[i+1])) if i<len(bins)-1 else "X > %d" % (bins[i])
        histProj[signal].append(hist2d[signal].ProjectionY( htname, (bins[i]-bins[0])/25+1, (bins[i+1]-bins[0])/25 if i<len(bins)-1 else hist2d[signal].GetNbinsX()+1))
        histProj[signal][i].SetTitle(htcut)
        if i<len(bins)-1:
            x [signal].append((bins[i]+bins[i+1])/2.0/(1000 if vs=="ht" else 1))
            dx[signal].append((bins[i+1]-bins[i])/2.0/(1000 if vs=="ht" else 1))
        else:
            x [signal].append((bins[i]/(1000.0 if vs=="ht" else 1)+(0.5 if vs=="ht" else 200)))
            dx[signal].append(0.5 if vs=="ht" else 200)
        rms [signal].append( histProj[signal][i].GetRMS())
        drms[signal].append( histProj[signal][i].GetRMSError())
        fit = histProj[signal][i].Fit("gaus","QS","",-200,200)
        sig [signal].append(fit.Parameter(2))
        dsig[signal].append(fit.ParError (2))
    
    graphRMS[signal] = ROOT.TGraphErrors(len(x[signal]),array('f',x[signal]),array('f',rms[signal]),array('f',dx[signal]),array('f',drms[signal]))
    graphSig[signal] = ROOT.TGraphErrors(len(x[signal]),array('f',x[signal]),array('f',sig[signal]),array('f',dx[signal]),array('f',dsig[signal]))
    graphStyle(graphRMS[signal], "", "M_{T2}^{gen} [GeV]" if "mt2" in vs else "H_{T} [TeV]", 
               "M_{T2} resolution [GeV]", 21, colors[signal])
    graphStyle(graphSig[signal], "", "M_{T2}^{gen} [GeV]" if "mt2" in vs else "H_{T} [TeV]", 
               "M_{T2} resolution [GeV]", 25, colors[signal])

#can = ROOT.TCanvas("can","can",700,700)
can = ROOT.TCanvas("can","can")
#can.SetTopMargin(0.08)
#can.SetRightMargin(0.06)
#can.SetLeftMargin(0.12)
drawn = False
for signal in signals:
    if vs=="mt2":
        graphRMS[signal].SetMaximum(120)
    graphRMS[signal].Draw("AP" if not drawn else "Psame")
    drawn = True
    graphSig[signal].Draw("Psame")

mark, lat, line = {},{},{}
mark[1] = ROOT.TMarker(0.43,0.25, 25);  lat[1] = ROOT.TLatex(0.45,0.25, "#sigma_{G}(M_{T2} - M_{T2}^{gen})")
mark[2] = ROOT.TMarker(0.43,0.20, 21);  lat[2] = ROOT.TLatex(0.45,0.20, "RMS(M_{T2} - M_{T2}^{gen})")
for i,s in enumerate(signals):
    line[s] = ROOT.TLine(0.67, 0.20+0.05*i, 0.69, 0.20+0.05*i); lat[s] = ROOT.TLatex(0.7, 0.2+0.05*i, signals2[i])
    line[s].SetLineColor(colors[s]); line[s].SetLineWidth(3)

def textStyle(text):
    text.SetTextSize(0.03)
    text.SetTextFont(42)
    text.SetTextAlign(12)

for m in mark.values():
    m.SetNDC()
    m.Draw()
for l in line.values():
    l.DrawLineNDC(l.GetX1(),l.GetY1(),l.GetX2(),l.GetY2())
for l in lat.values():
    l.SetNDC()
    textStyle(l)
    l.Draw()

pandolfiniWantsThis.Draw()

can.SaveAs("plots/mt2ResVS"+vs+".eps")
can.SaveAs("plots/mt2ResVS"+vs+".png")



