######echo "T1tttt 1900,1"
######python skimBabies.py /pnfs/psi.ch/cms/trivcat/store/user/mschoene/MT2production/80X/PostProcessed/signal2016_Nov02_pP_Nov16/skimAndPrune/ /scratch/mmasciov/signalPoints_2016/ "GenSusyMGluino==1900 && GenSusyMNeutralino==1"  --filter="T1tttt" --gfalProtocol="gsiftp"
######hadd /scratch/mmasciov/signalPoints_2016/T1tttt_1900_1.root /scratch/mmasciov/signalPoints_2016/SMS_T1tttt_*prune_skim.root 
######rm -f /scratch/mmasciov/signalPoints_2016/SMS_T1tttt_*prune_skim.root
###
###echo "T1tttt 1100,875"
###python skimBabies.py /pnfs/psi.ch/cms/trivcat/store/user/mschoene/MT2production/80X/PostProcessed/signal2016_Nov02_pP_Nov16/skimAndPrune/ /scratch/mmasciov/signalPoints_2016/ "GenSusyMGluino==1100 && GenSusyMNeutralino==875"  --filter="T1tttt" --gfalProtocol="gsiftp"
###hadd /scratch/mmasciov/signalPoints_2016/T1tttt_1100_875.root /scratch/mmasciov/signalPoints_2016/SMS_T1tttt_*prune_skim.root 
###rm -f /scratch/mmasciov/signalPoints_2016/SMS_T1tttt_*prune_skim.root
###
###echo "T2bb 1100,1"
###python skimBabies.py /pnfs/psi.ch/cms/trivcat/store/user/mschoene/MT2production/80X/PostProcessed/signal2016_Nov02_pP_Nov16/skimAndPrune/ /scratch/mmasciov/signalPoints_2016/ "GenSusyMSbottom==1100 && GenSusyMNeutralino==1"  --filter="T2bb" --gfalProtocol="gsiftp"
###hadd /scratch/mmasciov/signalPoints_2016/T2bb_1100_1.root /scratch/mmasciov/signalPoints_2016/SMS_T2bb_*prune_skim.root 
###rm -f /scratch/mmasciov/signalPoints_2016/SMS_T2bb_*prune_skim.root
###
###echo "T2bb 500,475"
###python skimBabies.py /pnfs/psi.ch/cms/trivcat/store/user/mschoene/MT2production/80X/PostProcessed/signal2016_Nov02_pP_Nov16/skimAndPrune/ /scratch/mmasciov/signalPoints_2016/ "GenSusyMSbottom==500 && GenSusyMNeutralino==475"  --filter="T2bb" --gfalProtocol="gsiftp"
###hadd /scratch/mmasciov/signalPoints_2016/T2bb_500_475.root /scratch/mmasciov/signalPoints_2016/SMS_T2bb_*prune_skim.root 
###rm -f /scratch/mmasciov/signalPoints_2016/SMS_T2bb_*prune_skim.root
###
###echo "T2qq 1500,1"
###python skimBabies.py /pnfs/psi.ch/cms/trivcat/store/user/mschoene/MT2production/80X/PostProcessed/signal2016_Nov02_pP_Nov16/skimAndPrune/ /scratch/mmasciov/signalPoints_2016/ "GenSusyMSquark==1500 && GenSusyMNeutralino==1"  --filter="T2qq" --gfalProtocol="gsiftp"
###hadd /scratch/mmasciov/signalPoints_2016/T2qq_1500_1.root /scratch/mmasciov/signalPoints_2016/SMS_T2qq_*prune_skim.root 
###rm -f /scratch/mmasciov/signalPoints_2016/SMS_T2qq_*prune_skim.root
###
###echo "T2qq 450,425"
###python skimBabies.py /pnfs/psi.ch/cms/trivcat/store/user/mschoene/MT2production/80X/PostProcessed/signal2016_Nov02_pP_Nov16/skimAndPrune/ /scratch/mmasciov/signalPoints_2016/ "GenSusyMSquark==450 && GenSusyMNeutralino==425"  --filter="T2qq" --gfalProtocol="gsiftp"
###hadd /scratch/mmasciov/signalPoints_2016/T2qq_450_425.root /scratch/mmasciov/signalPoints_2016/SMS_T2qq_*prune_skim.root 
###rm -f /scratch/mmasciov/signalPoints_2016/SMS_T2qq_*prune_skim.root



echo "T1tttt 1900,1"
python skimBabies.py /pnfs/psi.ch/cms/trivcat/store/user/mschoene/MT2production/80X/PostProcessed/signal2016_Nov02_pP_Nov16/ /scratch/mmasciov/signalPoints_2016_unskimmed/ "GenSusyMGluino==1900 && GenSusyMNeutralino==1"  --filter="T1tttt" --gfalProtocol="gsiftp"
hadd /scratch/mmasciov/signalPoints_2016_unskimmed/T1tttt_1900_1.root /scratch/mmasciov/signalPoints_2016_unskimmed/SMS_T1tttt_*skim.root 
rm -f /scratch/mmasciov/signalPoints_2016_unskimmed/SMS_T1tttt_*skim.root

echo "T1tttt 1100,875"
python skimBabies.py /pnfs/psi.ch/cms/trivcat/store/user/mschoene/MT2production/80X/PostProcessed/signal2016_Nov02_pP_Nov16/ /scratch/mmasciov/signalPoints_2016_unskimmed/ "GenSusyMGluino==1100 && GenSusyMNeutralino==875"  --filter="T1tttt" --gfalProtocol="gsiftp"
hadd /scratch/mmasciov/signalPoints_2016_unskimmed/T1tttt_1100_875.root /scratch/mmasciov/signalPoints_2016_unskimmed/SMS_T1tttt_*skim.root 
rm -f /scratch/mmasciov/signalPoints_2016_unskimmed/SMS_T1tttt_*skim.root

echo "T2bb 1100,1"
python skimBabies.py /pnfs/psi.ch/cms/trivcat/store/user/mschoene/MT2production/80X/PostProcessed/signal2016_Nov02_pP_Nov16/ /scratch/mmasciov/signalPoints_2016_unskimmed/ "GenSusyMSbottom==1100 && GenSusyMNeutralino==1"  --filter="T2bb" --gfalProtocol="gsiftp"
hadd /scratch/mmasciov/signalPoints_2016_unskimmed/T2bb_1100_1.root /scratch/mmasciov/signalPoints_2016_unskimmed/SMS_T2bb_*skim.root 
rm -f /scratch/mmasciov/signalPoints_2016_unskimmed/SMS_T2bb_*skim.root

echo "T2bb 500,475"
python skimBabies.py /pnfs/psi.ch/cms/trivcat/store/user/mschoene/MT2production/80X/PostProcessed/signal2016_Nov02_pP_Nov16/ /scratch/mmasciov/signalPoints_2016_unskimmed/ "GenSusyMSbottom==500 && GenSusyMNeutralino==475"  --filter="T2bb" --gfalProtocol="gsiftp"
hadd /scratch/mmasciov/signalPoints_2016_unskimmed/T2bb_500_475.root /scratch/mmasciov/signalPoints_2016_unskimmed/SMS_T2bb_*skim.root 
rm -f /scratch/mmasciov/signalPoints_2016_unskimmed/SMS_T2bb_*skim.root

echo "T2qq 1500,1"
python skimBabies.py /pnfs/psi.ch/cms/trivcat/store/user/mschoene/MT2production/80X/PostProcessed/signal2016_Nov02_pP_Nov16/ /scratch/mmasciov/signalPoints_2016_unskimmed/ "GenSusyMSquark==1500 && GenSusyMNeutralino==1"  --filter="T2qq" --gfalProtocol="gsiftp"
hadd /scratch/mmasciov/signalPoints_2016_unskimmed/T2qq_1500_1.root /scratch/mmasciov/signalPoints_2016_unskimmed/SMS_T2qq_*skim.root 
rm -f /scratch/mmasciov/signalPoints_2016_unskimmed/SMS_T2qq_*skim.root

echo "T2qq 450,425"
python skimBabies.py /pnfs/psi.ch/cms/trivcat/store/user/mschoene/MT2production/80X/PostProcessed/signal2016_Nov02_pP_Nov16/ /scratch/mmasciov/signalPoints_2016_unskimmed/ "GenSusyMSquark==450 && GenSusyMNeutralino==425"  --filter="T2qq" --gfalProtocol="gsiftp"
hadd /scratch/mmasciov/signalPoints_2016_unskimmed/T2qq_450_425.root /scratch/mmasciov/signalPoints_2016_unskimmed/SMS_T2qq_*skim.root 
rm -f /scratch/mmasciov/signalPoints_2016_unskimmed/SMS_T2qq_*skim.root

