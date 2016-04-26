import os

cfg="mc_QCD4RS"
mc="mc"
#sampleIDs = range(151,158,1)
sampleIDs = [157]
Njobs = 30

#queue = "long.q"
queue = "all.q"
jobsDir = "closureJobsApr20_gen/"
os.system("mkdir "+jobsDir)

wdir=os.environ["PWD"]

# make sure we have valid credentials
os.system("voms-proxy-init -voms cms")

for sID in sampleIDs:
    for job in range(Njobs):
        command = "qsub -q {q} -o {wdir}/{jDir}/job_id{sid}_{j}of{N}.out -e {wdir}/{jDir}/job_id{sid}_{j}of{N}.err closureJob.sh {cfg} {mc} {sid} {j} {N}".format(q=queue, wdir=wdir, jDir=jobsDir, sid=sID, j=job, N=Njobs, cfg=cfg, mc=mc)
        #print command
        os.system(command)
