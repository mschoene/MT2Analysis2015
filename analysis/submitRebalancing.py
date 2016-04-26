import os

cfg="mc_QCD4RS"
mc="mc"
#sampleIDs = range(151,158,1)
#sampleIDs = range(154,158,1)
sampleIDs = [154]
Njobs = 30
toSE = "true"

queue = "all.q"
jobsDir = "rebalanceJobsApr19/"
os.system("mkdir "+jobsDir)

wdir=os.environ["PWD"]

# make sure we have valid credentials
os.system("voms-proxy-init -voms cms")

for sID in sampleIDs:
    for job in range(Njobs):
        command = "qsub -q {q} -o {wdir}/{jDir}/job_id{sid}_{j}of{N}.out -e {wdir}/{jDir}/job_id{sid}_{j}of{N}.err rebalanceJob.sh {cfg} {mc} {sid} {j} {N} {toSE}".format(q=queue, wdir=wdir, jDir=jobsDir, sid=sID, j=job, N=Njobs, cfg=cfg, mc=mc, toSE=toSE)
        #print command
        os.system(command)
