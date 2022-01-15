import os,sys


failedTxt = open("here.txt")
file_list = list()
for line in failedTxt:
    file_ = (line.split(" ")[-1])[:-1]
    print file_
    file_list.append(file_)
    cmd = "rm -f condor_output/%s"%(file_)
    print cmd
    os.system(cmd)


print "Total file deleted = %d"%(len(file_list))
# sys.exit(0)

for file_ in file_list:
    #print file_
    file_attr = file_.split("_")
    cmd_log = ""
    cmd_sub = ""
    if(len(file_attr) == 4):
        cmd_log = "condor_13042021_%s_job__data_v16_v8.*"%(file_attr[1])
        cmd_sub = "data_submitCondor_pion_%03dgev_v16_v8_data"%(int(file_attr[1].strip("pion")))
    elif(len(file_attr) == 7):
        cmd_log = "condor_13042021_%s_job__FTFP_BERT_EMN_v46_patchMIP.*"%(file_attr[1])
        cmd_sub = "sim_submitCondor_pion_%03dgev_v46_patchMIP_FTFP_BERT_EMN"%(int(file_attr[1].strip("pion")))
    elif(len(file_attr) == 8):
        cmd_log = "condor_13042021_%s_job__QGSP_FTFP_BERT_EMN_v46_patchMIP.*"%(file_attr[1])
        cmd_sub = "sim_submitCondor_pion_%03dgev_v46_patchMIP_QGSP_FTFP_BERT_EMN"%(int(file_attr[1].strip("pion")))
    else:
        print "Unidentified file!!!"
        sys.exit(0)
    cmd = "rm -f condor_output/condor_logs/%s"%(cmd_log)
    print cmd
    os.system(cmd)
    cmd = "condor_submit condor_submit/%s"%(cmd_sub)
    print cmd
    os.system(cmd)
