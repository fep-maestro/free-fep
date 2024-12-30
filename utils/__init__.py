# main file
# Created 2024

__all__ = ['main', 'mdsetup',  'env']


import os
import configparser
import tempfile
import sys
import logging

ROOT_DIR = os.path.dirname( os.path.dirname(__file__) )
print('ROOT_DIR: ', ROOT_DIR)
#print('package =',__package__)
#print('__name__ =',__name__)
sys.path.append(ROOT_DIR)   # Add Parent dir of the file to PYTHONPATH (Note, the package dir must be "OmiX")

from .mdsetup import check_gmx
# read config values
config = configparser.ConfigParser()

tmpdir_conf = config.get('gmx', 'TMPDIR')
tmpdir_env = os.getenv('TMPDIR')
gmxdir = os.environ.get('CONDA_PREFIX')

################## gmxFE ####################
gmxBin      = config.get('gmx', 'gmx')
gmxBin_d    = config.get('gmx', 'gmx_d')
gmxBinDir   = config.get('gmx', 'GMXBINDIR')

import toml
# Read/Write TOML config  https://towardsdatascience.com/managing-deep-learning-models-easily-with-toml-configurations-fb680b9deabe#:~:text=The%20concept%20of%20a%20TOML,there%20are%20multiple%20nested%20levels.
# YAML format verificator: https://yaml-online-parser.appspot.com/
def load_global_config( filepath : str = 'config_gmxpy.toml' ):
    import munch
    return munch.munchify( toml.load( filepath ) )
def save_global_config( new_config , filepath : str = 'project_config.toml' ):
    with open( filepath , 'w' ) as file:
        toml.dump( new_config , file )
#############################################

##
##  GMX FE
##
class configure:
    """
    Class for processing, auditing and sharing configuration parameters of GMX MD runs
    """

    def __init__(self):
        #self.sPhase = 'Example'

        if tmpdir_env:
            self.sTmpDir = tmpdir_env
        elif tmpdir_conf and os.path.exists(tmpdir_conf):
            self.sTmpDir = tmpdir_conf
        else:
            self.sTmpDir = tempfile.gettempdir()

        self.TaskNm = 'Test'
        self.TopDir  = os.path.abspath('INPUT/TOP')
        self.MdpDir  = os.path.abspath('INPUT/MDP')
        self.runFile = os.path.abspath('INPUT/runList.txt')
        self.OutDir  = os.path.abspath('OUTPUT')
        self.ResFile = os.path.abspath('OUTPUT/ResFile.txt')
        #self.runFile = ''
        self.RunTIstates = ['all']
        self.InitPath = './'
        self.bVerbose = False
        self.sge_job_id = 0
        self.Templates = []
        self.Template = ''
        self.JobNames = []
        self.bReplex = 0                      # replace 'sReplex' -> 'bReplex'
        self.Type = 'run_ti'
        self.queue = 'none'
        self.SolvCoor = ''
        self.BoxSize = ''
        self.Lines = []

            
    ##
    ##  Read and do diagnostics of the CLI arguments and config file parameters
    ##  Parameters defined in CLI override those defined in the config
    ##
    def Init(self, args):
        """
        Process command line arguments

        :param args: CLI gmxfe arguments
        :return: none
        """
        if args.lines and (args.beginline or args.endline):
            logging.error(f'--lines and --beginline|--endline are mutually exclusive options')
            sys.exit()
        if args.queue == 'none' and hasattr(args, 'partition'):
            logging.error(f'The option --partition is used only with --queue (slurm|pbs|sge) but ignored with --queue=none')
            logging.error(f'Ignore the option --partition')
        self.InitPath = os.getcwd()
 
        # CLI arguments have higher priority, i.e. overwriting config params
        if not (hasattr(args, 'config') and args.config):
            logging.error('The mandatotory parameter --config is not defined. Exit!')
            sys.exit(0)
        cfgPath = os.path.realpath(args.config)
        if not os.path.isfile(cfgPath):
            logging.error(f'The config {cfgPath} is not found. Exit!')
            sys.exit(0)
        logging.info(f'Loading config: {cfgPath}')
        cfg = load_global_config(cfgPath)
        print('type(config): ', type(cfg))
        # print( toml.load( "config_gmxpy.toml" ) )
        #print( cfg.molsys1.LigTitle )
        #args.gmxconfig = cfg
        self.ProcessConfig(cfg, cfgPath)
        #logging.info(f'Config: {cfg}')
    ##
    ## 1. SETUP config parameters
    ## Define GMX binary paths
        if hasattr(args, 'gmx') and args.gmx:
            if os.path.isfile(args.gmx):
                self.sGmx = args.gmx
            elif gmxdir and os.path.isfile(os.path.join(gmxdir, 'bin', args.gmx)):
                self.sGmx = os.path.join(gmxdir, 'bin', args.gmx)
            else:
                self.sGmx = args.gmx
        if not self.sGmx:
           logging.error('The mandatotory parameter - gmx executable command is not definded. Exit!')
           sys.exit(0)
        logging.info(f'gmx executable command: \"{self.sGmx}\"')
        check_gmx(self.sGmx  , logging)

        if hasattr(args, 'gmx_d') and args.gmx_d:
            if os.path.isfile(args.gmx_d):
                self.sGmx_d = args.gmx_d
            elif gmxdir and os.path.isfile(os.path.join(gmxdir, 'bin', args.gmx_d)):
                self.sGmx_d = os.path.join(gmxdir, 'bin', args.gmx)
            else:
                self.sGmx_d = args.gmx_d
        if not self.sGmx_d:
            logging.warning('The double precision gmx executable command is not defined. Will use the single precision executable for minimization (not recommended)')
            self.sGmx_d = self.sGmx
        logging.info(f'gmx double precision executable command: \"{self.sGmx_d}\"')
        check_gmx(self.sGmx_d, logging)

        if hasattr(args, 'queue') and args.queue:
            self.queue = args.queue
        if hasattr(args, 'tmpdir') and \
                args.tmpdir and \
                os.path.isdir(args.tmpdir): self.sTmpDir = args.tmpdir

        if hasattr(args, 'runfile') and args.runfile: self.runFile = os.path.abspath(args.runfile)
        if not os.path.isfile(self.runFile):
            logging.error(f'The SysList file {self.runFile} is not defined! Exit')
            sys.exit()
        if hasattr(args, 'mdpdir') and args.mdpdir:
            self.MdpDir = os.path.abspath(args.mdpdir)
        if not os.path.isdir(self.MdpDir):
            logging.error(f'Specified MdpDir {self.MdpDir} does not exist. Exit!')
            sys.exit()
        if hasattr(args, 'topdir') and args.topdir:
            self.TopDir = os.path.abspath(args.topdir)
        if not os.path.isdir(self.TopDir):
            logging.error(f'Specified TopDir {self.TopDir} does not exist. Exit!')
            sys.exit()
        if hasattr(args, 'outdir') and args.outdir: self.OutDir = os.path.abspath(args.outdir)
        if hasattr(args, 'resfile') and args.resfile: self.ResFile = os.path.abspath(args.resfile)

        if hasattr(args, 'lambdas') and args.lambdas: self.RunTIstates = args.lambdas
        if hasattr(args, 'type') and args.type: self.Type = args.type
        # Check if labda-grid is specified
        if self.Type == 'run_ti':
            vdw_lambdas  = self.Task.get('DEFAULT','')['mdp'].get('vdw-lambdas')
            fep_lambdas  = self.Task.get('DEFAULT','')['mdp'].get('fep-lambdas')
            # #slCoul_Lambdas = self.mdp.get('general','').get('coul-lambdas')
            # vdw_lambdas  = self.mdp.get('general','').get('vdw-lambdas')
            # fep_lambdas  = self.mdp.get('general','').get('fep-lambdas')
            if not vdw_lambdas and not fep_lambdas:
                logging.error(f'The mandotory \"fep-lambdas\" or \"vdw-lambdas\" lambda grid parameters are not defined in the \"Task.DEFAULT\" section of the config {cfgPath}. Exit')
                #logging.error(f"The mandotory \"fep-lambdas\" or \"vdw-lambdas\" lambda grid parameters are not defined in the \"Task.general\" section of the config {cfgPath}. Exit")
                sys.exit(0)
        if args.replex and args.jobname == 'gmxTIp': args.jobname = 'gmxTIall'
        if hasattr(args, 'solvcoor') and args.solvcoor:
            self.SolvCoor = os.path.abspath(args.solvcoor)
            logging.info(f'SolvCoor = {self.SolvCoor}')
        if hasattr(args, 'boxsize') and args.boxsize: self.BoxSize = args.boxsize
        if hasattr(args, 'lines') and args.lines: self.Lines = args.lines
        if hasattr(args, 'verbose') and args.verbose: self.bVerbose = True


    ##
    ##  Read and do diagnostics of the config file
    ##
    def ProcessConfig(self, cfg, cfgPath):
        """
        Process dictionary obtained from configure file processing

        :param cfg: TOML config dict
        :param cfgPath: path to config.toml
        :return: none
        """
        cfgSupported = {'byPhase','bySolvent'}
        if 'cfgType' in cfg and cfg['cfgType']:
            cfgType      = cfg['cfgType']
            if cfgType in cfgSupported:
                self.cfgType = cfgType
            else:
                logging.error(f'The cfgType = \"{cfgType}\" is not supported. Currently only values \"{cfgSupported}\" are implemented')
                sys.exit(0)
        else:
            logging.error(f'The mondatory parameter \"cfgType\" is not defined in the config {cfgPath}. Exit')
            sys.exit(0)

        if 'cfgTitle' in cfg and cfg['cfgTitle']:
            self.cfgTitle = cfg['cfgTitle']
        ##
        ## 1) Read mandotory section Phase
        ##
        Phase = cfg['Phase']
        if not any(Phase.values()):
            logging.error(f'The mandotory \"Phase\" section is missing/empty in the config {cfgPath}. Exit')
            sys.exit(0)
        self.Phase = Phase
        #print(f'Phase = {Phase}')
        ##
        ## Find all Tasks in workflows in all phase definitions
        ##
        # A Munch is a subclass of dict; it supports all the methods a dict does: https://github.com/Infinidat/munch
        usedTasks = set()
        for key, val in Phase.items():
            if key in {'workflow','mdp'}: continue # Exclude service sections which are not a phase
            workflow = val.get('workflow')
            if workflow and len(workflow)>0:
                usedTasks.update(workflow)
            print(f'Phase.{key}.workflow = {workflow}')
        print(f'all Tasks: {usedTasks}') # Print list of Tasks found in workflows in all phase definitions

        ##
        ## 2) Read mandotory section Task
        ##
        Task = cfg['Task']
        if not any(Task.values()):
            logging.error(f'The mandotory \"Task\" section is missing/empty in the config {cfgPath}. Exit') # Task=DEFAULT is mondatory for Lambdas
            sys.exit(0)
        if not hasattr(Task, 'DEFAULT') or not hasattr(Task['DEFAULT'],'mdp') or not any(Task['DEFAULT']['mdp'].values()):
            logging.error(f'The mandotory \"Task.DEFAULT.mdp\" section is missing/empty in the config {cfgPath}. Exit') # Task=DEFAULT is mondatory for Lambdas
            sys.exit(0)
        vdw_lambdas  = Task['DEFAULT']['mdp'].get('vdw-lambdas')
        fep_lambdas  = Task['DEFAULT']['mdp'].get('fep-lambdas')
        if not vdw_lambdas and not fep_lambdas:
            logging.error(f'The mandotory lambda grid \"fep-lambdas\" or \"vdw-lambdas\" is not defined in the \"Task.DEFAULT.mdp\" section of the config {cfgPath}. Exit')
            sys.exit(0)
        if len(usedTasks) > 0:
            if not any(Task.values()):
                logging.error(f'The definition of Tasks: {usedTasks} is not found in the \"mdp\" section of the config {cfgPath}. Exit')
                sys.exit(0)
            # Check that all used Tasks are defined in config
            # TODO: Comment out bellow to generilize the input (for the case when the specific Phase with missing Task definition is not called)
            for task in usedTasks:
                if not task in Task.keys():
                    logging.error(f'The Task \"{task}\" used in Phase workflows is not defined in \"Task\" section of the config {cfgPath}. Exit')
                    sys.exit(0)
        self.Task = Task


        # ##
        # ## 2) Read and Diagnostics of mdp.Task definitions -> re-defined in the Task-section
        # ##
        # mdp = cfg['mdp']
        # if len(usedTasks) > 0:
        #     if not any(mdp.values()):
        #         logging.error(f"The definition of Tasks: {usedTasks} is not found in the \"mdp\" section of the config {cfgPath}. Exit")
        #         sys.exit(0)
        # for task in usedTasks:
        #     mdp.get(task)
        #     if not mdp.get(task):
        #         logging.error(f"The definition of the Task \"{task}\" is not found in the \"mdp\" section of the config {cfgPath}. Exit")
        #         sys.exit(0)
        #     print(f'mdp.{task}')
        # self.mdp = mdp
        # #sys.exit()


        ##
        ##  Process "paths" section  in the config
        ##
        paths = cfg['paths']
        if 'gmx'   in paths and paths.get('gmx'):
            sGmx      = paths['gmx']
            self.sGmx = sGmx
            if os.path.exists(sGmx):
                self.sGmx = os.path.abspath(sGmx)
            elif gmxdir and os.path.exists(os.path.join(gmxdir, 'bin', sGmx)):
                self.sGmx = os.path.join(gmxdir, 'bin', sGmx)                   # in conda package we have gmx exe in PATH
        if 'gmx_d' in paths and paths.get('gmx_d'):
            sGmx_d      = paths['gmx_d']
            self.sGmx_d = sGmx_d
            if os.path.exists(sGmx_d):
                self.sGmx_d = os.path.abspath(sGmx_d)
            elif gmxdir and os.path.exists(os.path.join(gmxdir, 'bin', sGmx_d)):
                self.sGmx = os.path.join(gmxdir, 'bin', sGmx_d)                  # in conda package we have gmx exe in PATH

        if 'TmpDir'  in paths and paths['TmpDir'] : self.sTmpDir  = paths['TmpDir']
        if 'OutDir'  in paths and paths['OutDir'] : self.OutDir  = os.path.abspath(paths['OutDir'])
        if 'MdpDir'  in paths and paths['MdpDir'] : self.MdpDir  = os.path.abspath(paths['MdpDir'])
        if 'TopDir'  in paths and paths['TopDir'] : self.TopDir  = os.path.abspath(paths['TopDir'])
        if 'ResFile' in paths and paths['ResFile']: self.ResFile = os.path.abspath(paths['ResFile'])
        if 'SysList' in paths and paths['SysList']: self.sSysList = os.path.abspath(paths['SysList'])
        if 'TIdir'   in paths and paths['TIdir']  : self.sTIdir   = paths['TIdir']
        if 'outnm'   in paths and paths['outnm']  : self.sOutnm   = paths['outnm']

        if 'lambdas' in cfg and cfg['lambdas']: self.RunTIstates = cfg['lambdas']
        if 'solvcoor' in cfg and cfg['solvcoor']: self.SolvCoor = cfg['solvcoor']
        if 'boxsize' in cfg and cfg['boxsize']: self.BoxSize = cfg['boxsize']

        #if "templates" in cfg and cfg["templates"]: self.Templates = cfg["templates"]
        #if "jobnames" in cfg and cfg["jobnames"]: self.JobNames = cfg["jobnames"]
        #if "lines" in cfg and cfg["lines"]: self.Lines = cfg["lines"]
        if 'verbose' in cfg: self.bVerbose = cfg['verbose']


env = configure()
