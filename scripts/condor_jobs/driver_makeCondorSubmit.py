import os,sys,re,fileinput,string,shutil

energy = [20,50,80,100,120,200,250,300]
# energy = [50,100,200,300]
type_ = ["data","sim"]
# type_ = ["sim"]

for en in energy:
    for ty in type_:
        cmd = "python makeCondorSubmit.py %d %s"%(en,ty)
        print cmd
        os.system(cmd)
