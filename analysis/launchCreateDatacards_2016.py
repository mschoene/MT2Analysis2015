import os
import sys
from sys import argv

if len(argv)>1:
    model = argv[1]
else:
    model = "T1tttt"

if len(argv)>2:
    label = argv[2]
else:
    label = ""


cfg = "dataETH_SnTMC_35p9ifb" #"data_2016_SnT_36p8_FixedWJets"  #data_2016_SnTMC_362ifb_18ifbUB_Signal" #"data_Run2016_7p7ifb"
stepSize = 5 if "T2cc" in model else 25

M = range(600 if "T1" in model else 100 if "T2tt"==model else 100 if "T2cc"==model else 300,2301 if "T1" in model else 1801 if "T2qq"==model else 500,stepSize)

os.system("mkdir jobs_"+model+"_2016")

for i,m in enumerate(M):
    m1 = m
    m2 = m1+stepSize
    #    m2 = m1+25
    #    Y = range(400,m1+1,25)
    Y = range(0,m1+1,stepSize)
    for j,y in enumerate(Y):
        m11=y
        m22=Y[j+1] if j+1 < len(Y) else y+stepSize    
    #m11=0
    #m22=m1+1
        if m11 >= m2: 
            continue
        
        command="qsub -l h_vmem=6g -q short.q -o `pwd`/jobs_"+model+"_2016/log_"+str(m1)+"_"+str(m2)+"_"+str(m11)+"_"+str(m22)+".out -e `pwd`/jobs_"+model+"_2016/log_"+str(m1)+"_"+str(m2)+"_"+str(m11)+"_"+str(m22)+".err -N creatingDatacards_"+model+"_"+str(m1)+"_"+str(m2)+"_"+str(m11)+"_"+str(m22)+" createDatacards_batch_2016.sh "+cfg+" "+model+" "+str(m1)+" "+str(m2)+" "+str(m11)+" "+str(m22)+" "+label+" "
        
        print command
        os.system(command)

       # command_mkdir="env --unset=LD_LIBRARY_PATH gfal-mkdir -p  gsiftp://t3se01.psi.ch/pnfs/psi.ch/cms/trivcat/store/user/`whoami`/EventYields_"+cfg+"/datacards_"+model+"_"+label+"/datacards_"+str(m1)+"_"+str(m11)
       # print command_mkdir
       # os.system(command_mkdir)


    # for j,y in enumerate(Y):
    #     m11=y
    #     m22=Y[j+1] if j+1 < len(Y) else y+25

    #     if m11 >= m2: 
    #         continue
    
    #     command="qsub -l h_vmem=6g -q short.q -o `pwd`/jobs_"+model+"_2016/log_"+str(m1)+"_"+str(m2)+"_"+str(m11)+"_"+str(m22)+".out -e `pwd`/jobs_"+model+"_2016/log_"+str(m1)+"_"+str(m2)+"_"+str(m11)+"_"+str(m22)+".err -N creatingDatacards_"+model+"_"+str(m1)+"_"+str(m2)+"_"+str(m11)+"_"+str(m22)+" createDatacards_batch_2016.sh "+cfg+" "+model+" "+str(m1)+" "+str(m2)+" "+str(m11)+" "+str(m22)+" "+label+" "
    #     print command
    #     os.system(command)

