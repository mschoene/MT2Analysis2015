import os, subprocess

cfg="mc_QCD4RS"
mc="mc"
sampleIDs = range(154,158,1)
#sampleIDs = [153]
Njobs = 30
toSE = "true"
fromSE = "true"
label = "genjets"

queue = "all.q"
jobsDir = "smearingJobsMay11_pandolfis_nogen_%s/" % label
os.system("mkdir "+jobsDir)

wdir=os.environ["PWD"]

# make sure we have valid credentials
out = subprocess.check_output("voms-proxy-info")
hoursleft = int([x.split(":")[1] for x in out.strip().replace(" ","").split("\n") if x.split(":")[0]=='timeleft'][0])
if hoursleft < 5:
    os.system("voms-proxy-init -voms cms")

for sID in sampleIDs:
    for job in range(Njobs):
        #if os.path.isfile("/pnfs/psi.ch/cms/trivcat/store/user/casal/analysisRS_genjets/EventYields_mc_QCD4RS_dummy/smearedTrees_rebGaus_smearHist_smearGen/mc_id%d_job%dof30.root"%(sID,job)):
        #    continue
        command = "qsub -q {q} -o {wdir}/{jDir}/job_id{sid}_{j}of{N}.out -e {wdir}/{jDir}/job_id{sid}_{j}of{N}.err smearJob.sh {cfg} {mc} {sid} {j} {N} {toSE} {fromSE} {lbl}".format(q=queue, wdir=wdir, jDir=jobsDir, sid=sID, j=job, N=Njobs, cfg=cfg, mc=mc, toSE=toSE, fromSE=fromSE, lbl=label)
        #print command
        os.system(command)
