import os, subprocess

cfg="mc_QCD4RS"
mc="mc"
sampleIDs = range(154,158,1)
#sampleIDs = [154]
Njobs = 30
label = "genjets"

#queue = "long.q"
queue = "all.q"
jobsDir = "closureJobsMay12_pandolfiGen_%s/" % (label) 
os.system("mkdir "+jobsDir)

wdir=os.environ["PWD"]

# make sure we have valid credentials
out = subprocess.check_output("voms-proxy-info")
hoursleft = int([x.split(":")[1] for x in out.strip().replace(" ","").split("\n") if x.split(":")[0]=='timeleft'][0])
if hoursleft < 5:
    os.system("voms-proxy-init -voms cms")

for sID in sampleIDs:
    for job in range(Njobs):
        command = "qsub -q {q} -o {wdir}/{jDir}/job_id{sid}_{j}of{N}.out -e {wdir}/{jDir}/job_id{sid}_{j}of{N}.err closureJob.sh {cfg} {mc} {sid} {j} {N} {lbl}".format(q=queue, wdir=wdir, jDir=jobsDir, sid=sID, j=job, N=Njobs, cfg=cfg, mc=mc, lbl=label)
        #print command
        os.system(command)
