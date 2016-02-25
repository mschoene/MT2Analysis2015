import os

cfg="mc_QCD4RS"
mc="mc"
sampleIDs = range(154,158,1)
Njobs = 20

queue = "long.q"
jobsDir = "rebalanceJobs/"
os.system("mkdir "+jobsDir)

wdir=os.environ["PWD"]

for sID in sampleIDs:
    for job in range(Njobs):
        command = "qsub -q {q} -o {wdir}/{jDir}/job_id{sid}_{j}of{N}.out -e {wdir}/{jDir}/job_id{sid}_{j}of{N}.err rebalanceJob.sh {cfg} {mc} {sid} {j} {N}".format(q=queue, wdir=wdir, jDir=jobsDir, sid=sID, j=job, N=Njobs, cfg=cfg, mc=mc)
        #print command
        os.system(command)
