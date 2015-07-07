CC = g++
CFLAGS = -Wall -c -g -fpermissive 

ROOFIT_INCLUDE := $(shell cd $(CMSSW_BASE); scram tool info roofitcore | grep INCLUDE= | sed 's|INCLUDE=||')
ROOFIT_LIBDIR := $(shell cd $(CMSSW_BASE); scram tool info roofitcore | grep LIBDIR= | sed 's|LIBDIR=||')

INCLUDES = -I. -I.. -I$(ROOTSYS)/include -I$(ROOFIT_INCLUDE)/ -I$(CMSSW_BASE)/src -I$(CMSSW_BASE)/src/CommonTools -I$(CMSSW_BASE)/src/UserCode/mmasciov/ -I$(CMSSW_BASE)/src/CondFormats/JetMETObjects/interface 

ROOTSYS ?= ERROR_RootSysIsNotDefined

ROOTFLAG = `${ROOTSYS}/bin/root-config --cflags --libs`

EXTRALIBS := -L$(ROOTSYS)/lib -L$(ROOFIT_LIBDIR)/ -L$(CMSSW_BASE)/lib/slc6_amd64_gcc481/ -L/swshare/cms/slc6_amd64_gcc481/cms/cmssw/CMSSW_7_0_6_patch3/lib/slc6_amd64_gcc481/ -L/swshare/cms/slc6_amd64_gcc481/cms/cmssw-patch/CMSSW_7_0_6_patch3/lib/slc6_amd64_gcc481 -lHtml -lMathCore -lGenVector -lMinuit -lEG  -lPhysics -L/swshare/cms/slc6_amd64_gcc481/external/zlib/1.2.3/lib -lz -lTreePlayer -lRooFitCore -lRooFit




all: computeLostLepton_CR printYield printYield_CR fakeQCDestimate createDatacards createDatacards_fromPhotons regionEventYields zllPurityTrees  zllPurity  zllPurityTopSplit zllContolRegion computeZinvFromZll computeZllofBackground


computeLostLepton_CR: computeLostLepton_CR.cpp MT2Region.o MT2Sample.o MT2EstimateSyst.o MT2Estimate.o MT2DrawTools.o
	$(CC) -Wall $(INCLUDES) -o computeLostLepton_CR computeLostLepton_CR.cpp MT2Region.o MT2Sample.o MT2EstimateSyst.o MT2Estimate.o MT2DrawTools.o $(ROOTFLAG) $(EXTRALIBS) -fpermissive


computeZinvFromZll: computeZinvFromZll.cpp MT2Region.o MT2Sample.o MT2EstimateSyst.o MT2Estimate.o MT2EstimateTree.o MT2Config.o 
	$(CC) -Wall $(INCLUDES) -o computeZinvFromZll computeZinvFromZll.cpp MT2Region.o MT2Sample.o MT2EstimateSyst.o MT2Estimate.o MT2EstimateTree.o MT2Config.o $(ROOTFLAG) $(EXTRALIBS)


computeZllofBackground: computeZllofBackground.cpp MT2Region.o MT2Estimate.o MT2EstimateTree.o  MT2Sample.o MT2EstimateSyst.o MT2Config.o
	$(CC) -Wall $(INCLUDES) -o computeZllofBackground computeZllofBackground.cpp MT2Region.o MT2Estimate.o MT2DrawTools.o MT2EstimateTree.o  MT2Sample.o MT2EstimateSyst.o MT2Config.o $(ROOTFLAG) $(EXTRALIBS)

computeZllGammaRatio: computeZllGammaRatio.cpp MT2Region.o MT2Estimate.o MT2EstimateTree.o  MT2Sample.o MT2EstimateSyst.o MT2Config.o
	$(CC) -Wall $(INCLUDES) -o computeZllGammaRatio computeZllGammaRatio.cpp MT2Region.o MT2Estimate.o MT2DrawTools.o MT2EstimateTree.o  MT2Sample.o MT2EstimateSyst.o MT2Config.o $(ROOTFLAG) $(EXTRALIBS)


computeZllGammaRatio_HT_split: computeZllGammaRatio_HT_split.cpp MT2Region.o MT2Estimate.o MT2EstimateTree.o  MT2Sample.o MT2EstimateSyst.o MT2Config.o
	$(CC) -Wall $(INCLUDES) -o computeZllGammaRatio_HT_split computeZllGammaRatio_HT_split.cpp MT2Region.o MT2Estimate.o MT2DrawTools.o MT2EstimateTree.o  MT2Sample.o MT2EstimateSyst.o MT2Config.o $(ROOTFLAG) $(EXTRALIBS)

zllControlRegion: zllControlRegion.cpp MT2Region.o MT2Estimate.o MT2EstimateTree.o  MT2Sample.o MT2EstimateSyst.o  MT2DrawTools.o
	$(CC) -Wall $(INCLUDES) -o zllControlRegion zllControlRegion.cpp MT2Region.o MT2Estimate.o MT2EstimateTree.o  MT2Sample.o MT2DrawTools.o  MT2EstimateSyst.o  $(ROOTFLAG) $(EXTRALIBS)

zllPurity: zllPurity.cpp MT2Region.o MT2Estimate.o MT2EstimateTree.o  MT2Sample.o MT2EstimateSyst.o MT2Config.o
	$(CC) -Wall $(INCLUDES) -o zllPurity zllPurity.cpp MT2Region.o MT2Estimate.o MT2DrawTools.o MT2EstimateTree.o  MT2Sample.o MT2EstimateSyst.o MT2Config.o $(ROOTFLAG) $(EXTRALIBS)

zllPurityTopSplit: zllPurityTopSplit.cpp MT2Region.o MT2Estimate.o MT2EstimateTree.o  MT2Sample.o MT2EstimateSyst.o MT2Config.o
	$(CC) -Wall $(INCLUDES) -o zllPurityTopSplit zllPurityTopSplit.cpp MT2Region.o MT2Estimate.o MT2DrawTools.o MT2EstimateTree.o  MT2Sample.o MT2EstimateSyst.o MT2Config.o $(ROOTFLAG) $(EXTRALIBS)

zllPurityTrees: zllPurityTrees.cpp MT2Region.o MT2Estimate.o MT2EstimateTree.o  MT2Sample.o MT2EstimateSyst.o MT2Config.o
	$(CC) -Wall $(INCLUDES) -o zllPurityTrees zllPurityTrees.cpp MT2Region.o MT2Estimate.o MT2DrawTools.o MT2EstimateTree.o  MT2Sample.o MT2EstimateSyst.o MT2Config.o $(ROOTFLAG) $(EXTRALIBS)


drawZll: drawZll.cpp MT2Region.o MT2Estimate.o MT2EstimateTree.o  MT2Sample.o MT2EstimateSyst.o MT2Config.o
	$(CC) -Wall $(INCLUDES) -o drawZll drawZll.cpp MT2Region.o MT2Estimate.o MT2DrawTools.o MT2EstimateTree.o  MT2Sample.o MT2EstimateSyst.o MT2Config.o $(ROOTFLAG) $(EXTRALIBS)



computeQCD_CR: computeQCD_CR.cpp MT2Region.o MT2Sample.o MT2EstimateQCD.o MT2Estimate.o
	$(CC) -Wall $(INCLUDES) -o computeQCD_CR computeQCD_CR.cpp MT2Region.o MT2Sample.o MT2EstimateQCD.o MT2Estimate.o $(ROOTFLAG) $(EXTRALIBS)

fakeQCDdata: fakeQCDdata.cpp MT2Region.o MT2Sample.o MT2EstimateQCD.o MT2Estimate.o
	$(CC) -Wall $(INCLUDES) -o fakeQCDdata fakeQCDdata.cpp MT2Region.o MT2Sample.o MT2EstimateQCD.o MT2Estimate.o $(ROOTFLAG) $(EXTRALIBS)

computeZinvFromGamma: computeZinvFromGamma.cpp MT2Region.o MT2Sample.o MT2Estimate.o MT2EstimateTree.o MT2EstimateSyst.o  MT2EstimateZinvGamma.o MT2DrawTools.o
	$(CC) -Wall $(INCLUDES) -o computeZinvFromGamma computeZinvFromGamma.cpp MT2Region.o MT2Sample.o MT2Estimate.o MT2EstimateTree.o MT2EstimateSyst.o  MT2EstimateZinvGamma.o MT2DrawTools.o $(ROOTFLAG) $(EXTRALIBS)

regionEventYields: regionEventYields.cpp MT2Region.o MT2Sample.o MT2EstimateTree.o MT2Estimate.o MT2DrawTools.o
	$(CC) -Wall $(INCLUDES) -o regionEventYields regionEventYields.cpp MT2Region.o MT2Sample.o MT2EstimateTree.o MT2Estimate.o MT2DrawTools.o $(ROOTFLAG) $(EXTRALIBS) -fpermissive

createDatacards: createDatacards.cpp MT2Region.o MT2Sample.o MT2EstimateSyst.o MT2Estimate.o  MT2DrawTools.o
	$(CC) -Wall $(INCLUDES) -o createDatacards createDatacards.cpp MT2Region.o MT2Sample.o MT2EstimateSyst.o MT2Estimate.o  MT2DrawTools.o $(ROOTFLAG) $(EXTRALIBS) -fpermissive

createDatacards_ZllProjection: createDatacards_ZllProjection.cpp MT2Region.o MT2Sample.o MT2EstimateSyst.o MT2Estimate.o  MT2DrawTools.o
	$(CC) -Wall $(INCLUDES) -o createDatacards_ZllProjection createDatacards_ZllProjection.cpp MT2Region.o MT2Sample.o MT2EstimateSyst.o MT2Estimate.o  MT2DrawTools.o $(ROOTFLAG) $(EXTRALIBS) -fpermissive

createDatacards_fromPhotons: createDatacards_fromPhotons.cpp MT2Region.o MT2Sample.o MT2EstimateSyst.o MT2Estimate.o  MT2DrawTools.o
	$(CC) -Wall $(INCLUDES) -o createDatacards_fromPhotons createDatacards_fromPhotons.cpp MT2Region.o MT2Sample.o MT2EstimateSyst.o MT2Estimate.o  MT2DrawTools.o $(ROOTFLAG) $(EXTRALIBS) -fpermissive

createDatacards_pureZll: createDatacards_pureZll.cpp MT2Region.o MT2Sample.o MT2EstimateSyst.o MT2Estimate.o  MT2DrawTools.o
	$(CC) -Wall $(INCLUDES) -o createDatacards_pureZll createDatacards_pureZll.cpp MT2Region.o MT2Sample.o MT2EstimateSyst.o MT2Estimate.o  MT2DrawTools.o $(ROOTFLAG) $(EXTRALIBS) -fpermissive


mt2Optimization: mt2Optimization.cpp MT2Region.o MT2Sample.o MT2EstimateSyst.o MT2Estimate.o MT2EstimateTree.o MT2DrawTools.o
	$(CC) -Wall $(INCLUDES) -o mt2Optimization mt2Optimization.cpp MT2Region.o MT2Sample.o MT2EstimateSyst.o MT2Estimate.o MT2EstimateTree.o MT2DrawTools.o $(ROOTFLAG) $(EXTRALIBS) -fpermissive

inclusivePlots: inclusivePlots.cpp MT2Region.o MT2Sample.o MT2EstimateTree.o MT2Estimate.o MT2DrawTools.o
	$(CC) -Wall $(INCLUDES) -o inclusivePlots inclusivePlots.cpp MT2Region.o MT2Sample.o MT2EstimateTree.o MT2Estimate.o MT2DrawTools.o $(ROOTFLAG) $(EXTRALIBS) -fpermissive

prova: prova.cpp MT2Region.o MT2Sample.o MT2Estimate.o MT2EstimateSyst.o
	$(CC) -Wall $(INCLUDES) -o prova prova.cpp MT2Region.o MT2Sample.o MT2Estimate.o MT2EstimateSyst.o $(ROOTFLAG) $(EXTRALIBS)

provaRead: provaRead.cpp MT2Region.o MT2Sample.o MT2Estimate.o MT2EstimateSyst.o
	$(CC) -Wall $(INCLUDES) -o provaRead provaRead.cpp MT2Region.o MT2Sample.o MT2Estimate.o MT2EstimateSyst.o $(ROOTFLAG) $(EXTRALIBS)


printYield: printYield.cpp MT2Region.o MT2Sample.o MT2Estimate.o MT2EstimateSyst.o
	$(CC) -Wall $(INCLUDES) -o printYield printYield.cpp MT2Region.o MT2Sample.o MT2Estimate.o MT2EstimateSyst.o $(ROOTFLAG) $(EXTRALIBS)

printYield_CR: printYield_CR.cpp MT2Region.o MT2Sample.o MT2Estimate.o MT2EstimateSyst.o
	$(CC) -Wall $(INCLUDES) -o printYield_CR printYield_CR.cpp MT2Region.o MT2Sample.o MT2Estimate.o MT2EstimateSyst.o $(ROOTFLAG) $(EXTRALIBS)

compareYield: compareYield.cpp MT2Region.o MT2Sample.o MT2Estimate.o MT2EstimateSyst.o MT2DrawTools.o
	$(CC) -Wall $(INCLUDES) -o compareYield compareYield.cpp MT2Region.o MT2Sample.o MT2Estimate.o MT2EstimateSyst.o  MT2DrawTools.o $(ROOTFLAG) $(EXTRALIBS)

fakeQCDestimate: fakeQCDestimate.cpp MT2Region.o MT2Sample.o MT2Estimate.o MT2EstimateSyst.o
	$(CC) -Wall $(INCLUDES) -o fakeQCDestimate fakeQCDestimate.cpp MT2Region.o MT2Sample.o MT2Estimate.o MT2EstimateSyst.o $(ROOTFLAG) $(EXTRALIBS)

makeZinvGammaTemplates: makeZinvGammaTemplates.cpp MT2Region.o MT2Sample.o MT2Estimate.o MT2EstimateZinvGamma.o
	$(CC) -Wall $(INCLUDES) -o makeZinvGammaTemplates makeZinvGammaTemplates.cpp MT2Region.o MT2Sample.o MT2Estimate.o MT2EstimateZinvGamma.o $(ROOTFLAG) $(EXTRALIBS)

makePromptTemplateData: makePromptTemplateData.cpp MT2Region.o MT2Sample.o MT2Estimate.o MT2EstimateZinvGamma.o MT2DrawTools.o
	$(CC) -Wall $(INCLUDES) -o makePromptTemplateData makePromptTemplateData.cpp MT2Region.o MT2Sample.o MT2Estimate.o MT2EstimateZinvGamma.o MT2DrawTools.o $(ROOTFLAG) $(EXTRALIBS)

drawZinvGamma: drawZinvGamma.cpp MT2Region.o MT2Estimate.o MT2DrawTools.o MT2EstimateTree.o
	$(CC) -Wall $(INCLUDES) -o drawZinvGamma drawZinvGamma.cpp MT2Region.o MT2Estimate.o MT2DrawTools.o MT2EstimateTree.o  $(ROOTFLAG) $(EXTRALIBS)


drawGammaTemplates: drawGammaTemplates.cpp MT2Region.o MT2Estimate.o MT2EstimateZinvGamma.o MT2DrawTools.o 
	$(CC) -Wall $(INCLUDES) -o drawGammaTemplates drawGammaTemplates.cpp MT2Region.o MT2Estimate.o MT2EstimateZinvGamma.o MT2DrawTools.o  $(ROOTFLAG) $(EXTRALIBS)

drawGammaTemplatesFromData: drawGammaTemplatesFromData.cpp MT2Region.o MT2Estimate.o MT2EstimateZinvGamma.o MT2DrawTools.o 
	$(CC) -Wall $(INCLUDES) -o drawGammaTemplatesFromData drawGammaTemplatesFromData.cpp MT2Region.o MT2Estimate.o MT2EstimateZinvGamma.o MT2DrawTools.o  $(ROOTFLAG) $(EXTRALIBS)

drawPurityFits: drawPurityFits.cpp MT2Region.o MT2Estimate.o MT2EstimateZinvGamma.o MT2DrawTools.o MT2EstimateSyst.o 
	$(CC) -Wall $(INCLUDES) -o drawPurityFits drawPurityFits.cpp MT2Region.o MT2Estimate.o MT2EstimateZinvGamma.o MT2DrawTools.o MT2EstimateSyst.o $(ROOTFLAG) $(EXTRALIBS)

drawQG: drawQG.cpp MT2Region.o MT2Estimate.o MT2EstimateTree.o MT2DrawTools.o 
	$(CC) -Wall $(INCLUDES) -o drawQG drawQG.cpp MT2Region.o MT2Estimate.o MT2EstimateTree.o MT2DrawTools.o  $(ROOTFLAG) $(EXTRALIBS)

drawGammaControlRegion: drawGammaControlRegion.cpp MT2Region.o MT2Estimate.o MT2EstimateSyst.o MT2DrawTools.o 
	$(CC) -Wall $(INCLUDES) -o drawGammaControlRegion drawGammaControlRegion.cpp MT2Region.o MT2Estimate.o MT2EstimateSyst.o MT2DrawTools.o  $(ROOTFLAG) $(EXTRALIBS)

fitPurityGamma: fitPurityGamma.cpp  MT2Region.o MT2Estimate.o MT2EstimateZinvGamma.o  MT2EstimateSyst.o
	$(CC) -Wall $(INCLUDES) -o fitPurityGamma fitPurityGamma.cpp MT2Region.o MT2Estimate.o MT2EstimateZinvGamma.o MT2EstimateSyst.o $(ROOTFLAG) $(EXTRALIBS)

fitIsoVsSietaieta: fitIsoVsSietaieta.cpp  MT2Region.o MT2Estimate.o MT2EstimateTree.o MT2Sample.o MT2DrawTools.o
	$(CC) -Wall $(INCLUDES) -o fitIsoVsSietaieta fitIsoVsSietaieta.cpp MT2Region.o MT2Estimate.o MT2EstimateTree.o MT2Sample.o MT2DrawTools.o $(ROOTFLAG) $(EXTRALIBS)

gammaControlRegion: gammaControlRegion.cpp MT2Region.o MT2Sample.o MT2Estimate.o MT2EstimateZinvGamma.o MT2EstimateSyst.o MT2EstimateTree.o
	$(CC) -Wall $(INCLUDES) -o gammaControlRegion gammaControlRegion.cpp MT2Region.o MT2Sample.o MT2Estimate.o MT2EstimateZinvGamma.o MT2EstimateSyst.o MT2EstimateTree.o $(ROOTFLAG) $(EXTRALIBS)

checkGenIso: checkGenIso.cpp MT2Region.o MT2Sample.o MT2Estimate.o MT2EstimateTree.o
	$(CC) -Wall $(INCLUDES) -o checkGenIso checkGenIso.cpp MT2Region.o MT2Sample.o MT2Estimate.o MT2EstimateTree.o $(ROOTFLAG) $(EXTRALIBS)

drawIsoROC: drawIsoROC.cpp MT2DrawTools.o
	$(CC) -Wall $(INCLUDES) -o drawIsoROC drawIsoROC.cpp MT2DrawTools.o $(ROOTFLAG) $(EXTRALIBS)

draw100pb: draw100pb.cpp MT2DrawTools.o
	$(CC) -Wall $(INCLUDES) -o draw100pb draw100pb.cpp MT2DrawTools.o $(ROOTFLAG) $(EXTRALIBS)


provaTree: provaTree.cpp MT2Region.o MT2Sample.o MT2Estimate.o MT2EstimateTree.o
	$(CC) -Wall $(INCLUDES) -o provaTree provaTree.cpp MT2Region.o MT2Sample.o MT2Estimate.o MT2EstimateTree.o $(ROOTFLAG) $(EXTRALIBS)

provaZinv: provaZinv.cpp MT2Region.o MT2Sample.o MT2Estimate.o MT2EstimateZinvGamma.o
	$(CC) -Wall $(INCLUDES) -o provaZinv provaZinv.cpp MT2Region.o MT2Sample.o MT2Estimate.o MT2EstimateZinvGamma.o $(ROOTFLAG) $(EXTRALIBS)

provaRoofit: provaRoofit.cpp 
	$(CC) -Wall $(INCLUDES) -o provaRoofit provaRoofit.cpp $(ROOTFLAG) $(EXTRALIBS)

provaTemplate: provaTemplate.cpp 
	$(CC) -Wall $(INCLUDES) -o provaTemplate provaTemplate.cpp $(ROOTFLAG) $(EXTRALIBS)


printLowYields: printLowYields.cpp MT2Region.o MT2Sample.o MT2Estimate.o MT2EstimateZinvGamma.o MT2EstimateSyst.o MT2EstimateTree.o
	$(CC) -Wall $(INCLUDES) -o printLowYields printLowYields.cpp MT2Region.o MT2Sample.o MT2Estimate.o MT2EstimateZinvGamma.o MT2EstimateSyst.o MT2EstimateTree.o $(ROOTFLAG) $(EXTRALIBS)

printPurity: printPurity.cpp MT2Region.o MT2Sample.o MT2Estimate.o MT2EstimateZinvGamma.o MT2EstimateSyst.o MT2EstimateTree.o
	$(CC) -Wall $(INCLUDES) -o printPurity printPurity.cpp MT2Region.o MT2Sample.o MT2Estimate.o MT2EstimateZinvGamma.o MT2EstimateSyst.o MT2EstimateTree.o $(ROOTFLAG) $(EXTRALIBS)

printLatexBGTable: printLatexBGTable.cpp MT2Region.o MT2Estimate.o 
	$(CC) -Wall $(INCLUDES) -o printLatexBGTable printLatexBGTable.cpp MT2Region.o MT2Estimate.o $(ROOTFLAG) $(EXTRALIBS)


computeTemplateSyst: computeTemplateSyst.cpp MT2Region.o MT2Estimate.o MT2DrawTools.o
	$(CC) -Wall $(INCLUDES) -o computeTemplateSyst computeTemplateSyst.cpp MT2Region.o MT2Estimate.o MT2DrawTools.o $(ROOTFLAG) $(EXTRALIBS)


translate8TeV: translate8TeV.cpp
	$(CC) -Wall $(INCLUDES) -o translate8TeV translate8TeV.cpp $(ROOTFLAG) $(EXTRALIBS) 


checkGammaSignalContamination: checkGammaSignalContamination.cpp MT2Region.o MT2Estimate.o 
	$(CC) -Wall $(INCLUDES) -o checkGammaSignalContamination checkGammaSignalContamination.cpp MT2Region.o MT2Estimate.o $(ROOTFLAG) $(EXTRALIBS)




MT2Sample.o: ../src/MT2Sample.cc 
	$(CC) -Wall $(INCLUDES) -c ../src/MT2Sample.cc $(ROOTFLAG) $(EXTRALIBS)

MT2Region.o: ../src/MT2Region.cc
	$(CC) -Wall $(INCLUDES) -c ../src/MT2Region.cc $(ROOTFLAG) $(EXTRALIBS)

MT2Estimate.o: ../src/MT2Estimate.cc MT2Region.o
	$(CC) -Wall $(INCLUDES) -c ../src/MT2Estimate.cc MT2Region.o $(ROOTFLAG) $(EXTRALIBS)

MT2EstimateSyst.o: ../src/MT2EstimateSyst.cc MT2Region.o MT2Estimate.o
	$(CC) -Wall $(INCLUDES) -c ../src/MT2EstimateSyst.cc MT2Region.o MT2Estimate.o $(ROOTFLAG) $(EXTRALIBS)

MT2EstimateZinvGamma.o: ../src/MT2EstimateZinvGamma.cc MT2Region.o MT2Estimate.o
	$(CC) -Wall $(INCLUDES) -c ../src/MT2EstimateZinvGamma.cc MT2Region.o MT2Estimate.o $(ROOTFLAG) $(EXTRALIBS)

MT2EstimateQCD.o: ../src/MT2EstimateQCD.cc MT2Region.o MT2Estimate.o
	$(CC) -Wall $(INCLUDES) -c ../src/MT2EstimateQCD.cc MT2Region.o MT2Estimate.o $(ROOTFLAG) $(EXTRALIBS)

MT2EstimateTree.o: ../src/MT2EstimateTree.cc MT2Region.o MT2Estimate.o
	$(CC) -Wall $(INCLUDES) -c ../src/MT2EstimateTree.cc MT2Region.o MT2Estimate.o $(ROOTFLAG) $(EXTRALIBS)

MT2Config.o: ../src/MT2Config.cc MT2Region.o MT2Estimate.o
	$(CC) -Wall $(INCLUDES) -c ../src/MT2Config.cc MT2Region.o MT2Estimate.o $(ROOTFLAG) $(EXTRALIBS)

MT2DrawTools.o: ../src/MT2DrawTools.cc 
	$(CC) -Wall $(INCLUDES) -c ../src/MT2DrawTools.cc $(ROOTFLAG) $(EXTRALIBS)






localQGLikelihoodCalculator2.o: $(CMSSW_BASE)/src/QGDev/qgMiniTuple/localQGLikelihoodCalculator/localQGLikelihoodCalculator2.cc
	$(CC) -Wall $(INCLUDES) -c $(CMSSW_BASE)/src/QGDev/qgMiniTuple/localQGLikelihoodCalculator/localQGLikelihoodCalculator2.cc $(ROOTFLAG) $(EXTRALIBS)

localQGLikelihoodCalculator.o: $(CMSSW_BASE)/src/QGDev/qgMiniTuple/localQGLikelihoodCalculator/localQGLikelihoodCalculator.cc
	$(CC) -Wall $(INCLUDES) -c $(CMSSW_BASE)/src/QGDev/qgMiniTuple/localQGLikelihoodCalculator/localQGLikelihoodCalculator.cc $(ROOTFLAG) $(EXTRALIBS)


# this is a comment
clean:
	rm *.o computeLostLepton computeLostLepton_CR createDatacards fakeQCDestimate printYield printYield_CR regionEventYields 

