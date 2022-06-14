import os
import shutil
import json


def readConfigFile(fn):
    """
    
    """
    
    cfgr = json.load(open(fn))
    return cfgr

def createProject(fn):
    """
    
    """
        
    cfgr = readConfigFile(fn)
    print("creating Directory "+ cfgr["experiment_id"])
    os.makedirs(cfgr["experiment_id"],exist_ok=True)
    
    print("Creating Directory "+ cfgr["experiment_id"]+cfgr["folder_output"])
    os.makedirs(cfgr["experiment_id"]+cfgr["folder_output"],exist_ok=True)
    
    print("Creating Directory "+ cfgr["experiment_id"]+cfgr["folder_semivariances"])
    os.makedirs(cfgr["experiment_id"]+cfgr["folder_semivariances"],exist_ok=True)
    
    print("Copy "+ fn + " -> " + cfgr["experiment_id"]+fn)
    shutil.copyfile(fn, cfgr["experiment_id"]+fn)
    