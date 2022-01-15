import os,sys,re,fileinput,string,shutil

type_ = sys.argv[1]


version = ""

if(type_ == "data"):
    sample = ["data"]
    version = "v16_v8"

elif(type_ == "sim"):
    sample = ["FTFP_BERT_EMN","QGSP_FTFP_BERT_EMN"]
    # version = "v44_VtxBeam_v3_correctFH10"
    version = "v46_patchMIP"
else:
    print "Incorrect sample"
    sys.exit()

energy = [20,50,80,100,120,200,250,300]
#energy = [50,80,120,250]

count = 0

for en in energy:
    for samp in sample:
        cmd = "condor_submit condor_submit/%s_submitCondor_pion_%03dgev_%s_%s"%(type_,en,version,samp)
        print cmd
        os.system(cmd)
        count += 1


print "\n\n Total jobs submitted :  ",count
