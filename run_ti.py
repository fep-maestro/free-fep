# Created 2024
import logging
import os
import sys
#import re
#import time
#import subprocess
#from src import sge
from utils import mdsetup as Set
import utils as utils


def MD_TI(runRecord, args):
    """
    MAIN ROUTINE to run MD and TI trajectories
    :param runRecord: record/line from the runFile with jobnm, lignm, phase and etc
    :param args: arguments of gmxfe command
    :return:
    """

    tokens = runRecord.split()
    logging.debug(f'Start working on the input line: "{runRecord}"')
    if len(tokens) < 3:
        logging.error(f'Record {runRecord} is empty or has wrong format. '
                          f'It must have 4 tokens, see gmxFE manual.')
        return
    # Parse input record:
    jobnm  = tokens[0]
    sLigNm = tokens[1]
    sPhase = tokens[2]
    utils.env.TaskNm = jobnm
    # read XML config template names for jobs in the same order as lines order
    #utils.env.Template = compare_lists(utils.env.Templates, list_lines, idx)

    ##
    ##  Setup Simulation Box
    ##
    cfgType = utils.env.cfgType
    #cfgType = args.gmxconfig.get('cfgType')
    if cfgType == 'byPhase':
        JobArgs        = Set.job_init(jobnm, sLigNm, sPhase, args)
        
        job_dir = JobArgs['job_dir']
        os.chdir(job_dir)  # Change dir to the job dir
    
        JobArgs        = Set.top_coor(sLigNm, sPhase, JobArgs, args)
        #sys.exit(0)

    elif cfgType == 'bySolvent':
        logging.error(f'cfgType == \'bySolvent\' is not implemented, yet. Exit')
        sys.exit(0)
    else:
        logging.error(f'cfgType == \'{cfgType}\' is not implemented, yet. Exit')
        sys.exit(0)


    JobArgs        = Set.mdp_Tasks(sLigNm, sPhase, job_dir, JobArgs, args)
    workflow       = Set.Tasks_to_workflow_script(JobArgs, args)

    # Generate folder structure for MD or TI runs
    if args.type == 'run_md':
        JobArgs        = Set.MD(JobArgs, args)
    elif args.type == 'run_ti':
        JobArgs        = Set.TI_grid(JobArgs, args)

    jobscripts = Set.submit_scripts(workflow, JobArgs, args)
    Set.schedule_job(jobscripts, JobArgs, args)

    os.chdir(utils.env.InitPath)
    #sys.exit(0)
    return()

def RunTI(jobnm, sLigNm, sPhase, args):
    """
    MAIN RunTI ROUTINE to run TI trajectory
    :param jobnm: name of the computation job
    :param sLigNm: ligand name
    :param sPhase: phase defined in config (e.g. Gas, Water, Protein)
    :param args: arguments of gmxfe command
    :return:
    """
    JobArgs        = Set.job_init(jobnm, sLigNm, sPhase, args)
    
    job_dir = JobArgs['job_dir']
    os.chdir(job_dir)  # Change dir to the job dir

    JobArgs        = Set.top_coor(sLigNm, sPhase, JobArgs, args)
    #sys.exit(0)

    JobArgs        = Set.mdp_Tasks(sLigNm, sPhase, job_dir, JobArgs, args)
    workflow       = Set.Tasks_to_workflow_script(JobArgs, args)


    # If TI then
    JobArgs        = Set.TI_grid(JobArgs, args)



    jobscript = Set.gen_submit_scripts(workflow, JobArgs, args)
    Set.submit_job(jobscript, JobArgs, args)

    os.chdir(utils.env.InitPath)
    #sys.exit(0)
    return()


def RunMD(args):
    """
    ROUTINE to run any simulation with provided gmxFE submit script.
    Requires fully prepared compute folder. Path to files should be
    absolute or relative to the config's location. 
    
    :return:
    """

    tokens = runRecord.split()
    logging.debug(f'Submitting the job script ...')
