import os

cfg="mc_QCD4RS"
mc="mc"
sampleIDs = range(151,158,1)
Njobs = 25
isBatch = "true"

queue = "all.q"
jobsDir = "smearJobsApr05/"
os.system("mkdir "+jobsDir)

wdir=os.environ["PWD"]

for sID in sampleIDs:
    for job in range(Njobs):
        command = "qsub -q {q} -o {wdir}/{jDir}/job_id{sid}_{j}of{N}.out -e {wdir}/{jDir}/job_id{sid}_{j}of{N}.err smearJob.sh {cfg} {mc} {sid} {j} {N} {isbatch}".format(q=queue, wdir=wdir, jDir=jobsDir, sid=sID, j=job, N=Njobs, cfg=cfg, mc=mc, isbatch=isBatch)
        #print command
        os.system(command)
