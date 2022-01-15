import os,sys,re,fileinput,string,shutil
from datetime import date


energy = int(sys.argv[1])
type_ = sys.argv[2]

version = ""

if(type_ == "data"):
    sample = ["data"]
    version = "v16_v8"

elif(type_ == "sim"):
    sample = ["FTFP_BERT_EMN","QGSP_FTFP_BERT_EMN"]
    #version = "v44_VtxBeam_v3_correctFH10"
    version = "v46_patchMIP"
else:
    print "Incorrect sample"
    sys.exit()
    
particle = "pion"
if(int(energy) == 420):
    particle = "muon"
    energy = "200"

for i,samp in enumerate(sample):
    condorSubmit = "condor_submit/%s_submitCondor_%s_%03dgev_%s_%s"%(type_,particle,energy,version,samp)
    if(type_ == "data"):
        file_ = "./file_list/data/v16_v8/%s%d.txt"%(particle,energy)
        output_ = "/home/work/spandey/public/HGCAL/CERN_TB/octTB/CMSSW_9_3_0/src/pion_analysis/newTree_showerShapes_inclusive_SS_faster/condor_jobs/condor_output/%s_%s%d_%s.root"%(type_,particle,energy,version)
    else:
        file_ = "./file_list/sim/%s/%s/%s%d.txt"%(version,samp,particle,energy)    
        output_ = "/home/work/spandey/public/HGCAL/CERN_TB/octTB/CMSSW_9_3_0/src/pion_analysis/newTree_showerShapes_inclusive_SS_faster/condor_jobs/condor_output/%s_%s%d_%s_%s.root"%(type_,particle,energy,samp,version)

    sample_ = type_


    print (str(date.today()))
    today = str(date.today()).split('-')[2]+str(date.today()).split('-')[1]+str(date.today()).split('-')[0]
    print (today)
    
    jobName      = today+"_"+particle+str(energy)+"_job_"+"_"+samp+"_"+version
    

    
    shutil.copyfile("proto_condor_submit",condorSubmit)
    for line in fileinput.FileInput(condorSubmit, inplace=1):
        line=line.replace("FILELIST", file_)
        line=line.replace("OUTPUT", output_)
        line=line.replace("SAMPLE", sample_)
        line=line.replace("ENERGY", str(energy))
        line=line.replace("JOBNAME", jobName)
        
        print line.rstrip()
            
        submitCommand = "condor_submit "+condorSubmit
        i = i+1

# print "Final NFilesDone ", NFilesDone
