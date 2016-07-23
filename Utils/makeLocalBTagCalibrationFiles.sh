#!/bin/bash

# This is a workaround to avoid conflicts between the BTagCalibration classes in these files and those
# which comes with the release

cat $CMSSW_RELEASE_BASE/src/CondTools/BTau/test/BTagCalibrationStandalone.cpp | \
    sed 's#BTagCalibration#BTagCalibrationStandalone#g' | sed 's#BTagEntry#BTagEntryStandalone#g'  > BTagCalibrationStandalone.cc

sed -i 's#StandaloneStandalone#Standalone#g' BTagCalibrationStandalone.cc

cat $CMSSW_RELEASE_BASE/src/CondTools/BTau/test/BTagCalibrationStandalone.h | \
    sed 's#BTagCalibration#BTagCalibrationStandalone#g' | sed 's#BTagEntry#BTagEntryStandalone#g'  > BTagCalibrationStandalone.h