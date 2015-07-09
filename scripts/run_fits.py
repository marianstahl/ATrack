'''
Submit with python submit.py <jobname>
(<jobname> is something like v1-simple_nsim_dg_wf_cheb3_PT_1D , i.e. no '' or "")
'''
import os
import sys
import logging
import time
import socket
import subprocess

jobname = str(sys.argv[1])
version, systematic = jobname.split("-")

for i in range(1,13):

    logfile = '/afs/cern.ch/work/m/mstahl/public/dumpdir/tracking_asymmetry_raw/{0}/Logs/{1}_{2}.log'.format(version,systematic,str(i))

    logging.basicConfig(level=logging.INFO,
                    format='%(message)s',
                    filename=logfile,
                    filemode='w')

    logger = logging.getLogger(str(i))
    ch = logging.StreamHandler()#sys.stdout
    formatter = logging.Formatter('%(message)s')
    ch.setFormatter(formatter)
    logger.addHandler(ch)
    
    fh = logging.FileHandler(logfile)
    fh.setFormatter(formatter)
    logger.addHandler(fh)
    
    logger.info( "================================================================================================================================")
    logger.info( "Starting on\t\t: " + time.strftime("%c"))
    logger.info( "Running on node\t\t: " + socket.gethostname())
    logger.info( "Current directory\t: "  + os.getcwd())
    logger.info( "Current job ID \t\t: manual submission")
    logger.info( "Current job name\t: " + jobname)
    logger.info( "Task index number\t: " + str(i))
    logger.info( "================================================================================================================================")
    logger.info( subprocess.check_output(["./bin/ATrack.exe",str(i),jobname]) )
    logger.info( "================================================================================================================================")
    logger.info( "Finished on\t\t: " + time.strftime("%c"))
    logger.info( "================================================================================================================================")
