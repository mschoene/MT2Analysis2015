#!/bin/env python

import os

def copyToSE(inputPath, sePath):
    if os.path.isdir(inputPath):
        fileList = [ f for f in os.listdir(inputPath) if ( os.path.isfile(os.path.join(inputPath, f)) and '.root' in f ) ]
        print '# files to be copied to SE: ', len(fileList)
        for f in fileList:
            if os.path.isfile(os.path.join("/pnfs/psi.ch/cms/trivcat/store/user/"+sePath, f)):
                print 'File already exist on SE! Going to next one...'
                continue
            else:
                print '\nCopying file ', inputPath+f, 'to /pnfs/psi.ch/cms/trivcat/store/user/'+sePath
                os.system("env -i X509_USER_PROXY=~/.x509up_u`id -u` gfal-copy file://"+inputPath+f+" srm://t3se01.psi.ch/pnfs/psi.ch/cms/trivcat/store/user/"+sePath+f)

    else: 
        print 'Input directory does NOT exist! Exiting...'
        exit
    
    return 1
