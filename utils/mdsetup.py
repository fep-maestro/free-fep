# Created 2024
import logging
import os
import re
import sys
import time
import xml.etree.ElementTree as ET
from utils import scheduler
import utils as utils
import subprocess

def read_file(path):                        # Read file to list of lines
    try:
        f = open(path, 'r')
    except OSError:
        logging.error(f"Error: {path} doesn't exist! Exit")
        sys.exit()
    lines = f.readlines()                    # Read lines with "\n" at the end
    # with open(path, 'r') as f:
    #     lines = f.read().splitlines()      # Read lines without "\n" at the end
    f.close()
    return lines
def write_file(lines, path):                        # Save lines to file
    try:
        f = open(path, 'w')
    except OSError:
        logging.error(f"Error: cannot open file for writing {path}\nExit")
        sys.exit()
    f.writelines(lines)                    # Read lines with "\n" at the end
    # with open(MdpTemplatePath, 'r') as f:
    #     lines = f.read().splitlines()      # Read lines without "\n" at the end
    f.close()

def copy_file_loc_2dist(file, location, destination):
    """
    Copy file from location to destination and returns the new path.

    :param        file: file path (full or relative) to be copied
    :param    location: full directory path of the file to be copied
    :param destination: full destination path of the file to be copied
    :return: full destination path of the copied file
    """
    fnm = os.path.basename(file)
    if not os.path.exists(file) and location: # it is a relative path
        file = os.path.join(location, file)
    try:
        import shutil
        shutil.copy2(file, destination)
    except FileNotFoundError:
        logging.error(f'The file {file} is not found.')
        sys.exit(0)
    except:
        logging.error(f'The file {file} cannot be copied.')
        sys.exit(0)
    return os.path.join(destination, fnm)
##
## USE --qargs to set all extra options and flags
##

def key2re(key):
    return re.sub(r'([^a-zAZ0-9 ])', r'\1', key)
def setMdpParams(parDict, lines):
    if len(parDict) == 0: return lines
    for key, value in parDict.items():
        sVal = str(value)
        new_lines =[]
        for l in lines:
            new_l = re.sub(r'^([ \t]*' + re.escape(key) + '[ \t]*\=[ \t]*)[^;\n]+', '\g<1>' + sVal + ' ', l)
            new_lines.append(new_l)
            #if (re.match(r'^([ \t]*' + re.escape(key) + '[ \t]*\=[ \t]*)[^;\n]+', l)) :
            #    print(f'\'{key}\': {value} : {l}')
            #    print(f'\'{key}\': {value} : {new_l}')
        lines = new_lines
    return new_lines
def setMdpParam(key, value, lines):
    sVal = str(value)
    new_lines =[]
    # print('^([ \t]*' + re.escape(key) + '[ \t]*\=[ \t]*)[^ \t;]+', '\\1' + re.escape(sVal))
    for l in lines:
        new_l = re.sub(r'^([ \t]*' + re.escape(key) + '[ \t]*\=[ \t]*)[^;\n]+', '\g<1>' + sVal, l)
        new_lines.append(new_l)
    #     matchline = re.match(r'^([ \t]*' + re.escape(key) + '[ \t]*\=[ \t]*)[^; \t]+', l)
    #     if matchline:
    #         print ('MATCH: ', matchline )
    #         print('SUBSTITUTED: ', re.sub(r'^([ \t]*' + re.escape(key) + '[ \t]*\=[ \t]*)[^; \t]+', '\g<1>' + sVal, l) )
    # return [re.sub(r'^([ \t]*' + re.escape(key) + '[ \t]*\=[ \t]*)[^; \t]+', '\g<1>' + sVal, l) for l in lines]
    return new_lines

def copy_parse_topology(sLigNm,TopPath,job_dir,sPhase,args,log):
    ##
    ## Set System Topology
    ##
    phase_cfg = utils.env.Phase.get(sPhase)
    TopTemplatePath =  os.path.abspath(utils.env.TopDir + '/' + phase_cfg.get('TopTemplate','')) 
    #TopTemplatePath =  os.path.abspath(utils.env.TopDir + '/' + args.gmxconfig.Phase[sPhase].get('TopTemplate','')) 
    log.info(f"Setting up the system topology from template {TopTemplatePath}")
    top_template = read_file(TopTemplatePath) 
    top          = resub({'{ligand}': sLigNm}, top_template)
    write_file(top, TopPath)
    # Copy all include files
    #print("Molecules:", args.gmxconfig.Phase.get(sPhase).get('Molecule'))
    for mol in phase_cfg.get('Molecule',''):
        #print(mol)
        #if len(mol.get('include_files')) < 1:
        lIncludeFiles = mol.get('include_files','')
        print(f"Copying include files for molecule ", mol.get('name'),': ', lIncludeFiles)
        if len(lIncludeFiles) > 0 and type(lIncludeFiles) is not list: lIncludeFiles = [lIncludeFiles]
        #print ("lIncludeFiles:",lIncludeFiles)
        for fnm in lIncludeFiles:
            if fnm == 'ligand.itp' and mol.get('type') == 'ligand': fnm = sLigNm +'.itp'  # Substitute the current ligand name
            #print ("fnm:",fnm)
            copy_file_loc_2dist(fnm, utils.env.TopDir, job_dir)

def resub(subDict, lines):
    for key, value in subDict.items():
        sVal = str(value)
        #print("Pattern:",re.escape(key))
        new_lines =[]
        for l in lines:
            new_l = re.sub(r'' + re.escape(key), sVal, l)
            new_lines.append(new_l)
        lines = new_lines
    return new_lines

def parse_lig_coor(ligand_lines, fullsys_lines, re_pattrn:str):
    #print('re_pattrn: ', re_pattrn)
    new_lines =[]
    after_lig =[]
    bBeforeLig = True
    # Collect lines from fullsys_lines before ligand resnm
    for l in fullsys_lines:
        if re.search(r'' + re_pattrn, l):
            bBeforeLig = False
            if len(after_lig) > 0:
                logging.error("Found the second ligand entry in the template structure:\n%s\nThere should be only single ligand entry.",l)
                sys.exit(0)
        else:
            if bBeforeLig: new_lines.append(l)
            else:          after_lig.append(l)
    if bBeforeLig:
        logging.error("Found no ligand entry in the structure. There should be one ligand entry.")
        sys.exit(0)

    # Insert ligand_lines
    for l in ligand_lines:
        if re.search(r'' + re_pattrn, l): new_lines.append(l)
    # append after_lig lines
    for l in after_lig:
        new_lines.append(l)

    return new_lines
def resnm_after_lig(coor_all, lignm:str, re_pattrn:str):
    #print('re_pattrn: ', re_pattrn)
    resnames_after_lig =[]
    bAfterLig = False
    # Collect lines from fullsys_lines before ligand resnm
    for l in coor_all:
        (resnm, nsubs) = re.subn(r''+re_pattrn, '\g<1>', l)  # (?=ATOM|HETATM) - POSITIVE LOOKAHEAD Zero-Length Assertions https://www.regular-expressions.info/lookaround.html
        if nsubs == 0: continue # skip the record in no substitutions
        resnm = resnm.strip()
        if resnm == lignm:
            bAfterLig = True
            #print('Resnm:', resnm, l)
            if len(resnames_after_lig) > 0:
                logging.error("Found the second ligand entry in the template structure:\n%s\nThere should be only single ligand entry.",l)
                sys.exit(0)
        elif bAfterLig:
            resnames_after_lig.append(resnm)
    #print (resnames_after_lig)
    uniqres = list(set(resnames_after_lig)) # sort| unique - sorted(set(list)) https://stackoverflow.com/questions/2931672/what-is-the-cleanest-way-to-do-a-sort-plus-uniq-on-a-python-list
    res_after_lig = {resnm:resnames_after_lig.count(resnm) for resnm in uniqres} # Form dictionary of resnm:nat
    #print ('Unique resnames after lig:', res_after_lig)
    return res_after_lig

def subst_lig_coor(coor_all,re_pattrn:str,lignm:str, coor_lig,logging):
    coor_new = []
    coor_after_lig =[]
    resnames_after_lig =[]
    bAfterLig = False
    resnm_prev = ''
    # Collect lines from fullsys_lines before ligand resnm
    for l in coor_all:
        (resnm, nsubs) = re.subn(r''+re_pattrn, '\g<1>', l)  # (?=ATOM|HETATM) - POSITIVE LOOKAHEAD Zero-Length Assertions https://www.regular-expressions.info/lookaround.html
        if nsubs == 0:
            if resnm_prev == lignm: continue  # Skip ligand records
            if bAfterLig: coor_after_lig.append(l)
            else: coor_new.append(l)
            continue # skip resnm analysis if no resnm substitution
        resnm = resnm.strip()
        if resnm == lignm:
            bAfterLig = True
            #print('Resnm:', resnm, l)
            if len(resnames_after_lig) > 0:
                logging.error("Found the second ligand entry in the template structure:\n%s\nThere should be only single ligand entry.",l)
                sys.exit(0)
        else:
            if bAfterLig:
                coor_after_lig.append(l)
                resnames_after_lig.append(resnm)
            else:
                coor_new.append(l)
        resnm_prev = resnm

    if not bAfterLig:
        logging.error("Found no entry with the ligand resnm \"{lignm}\" in the template structure. There must be one ligand entry.")
        sys.exit(0)
    uniqres = list(set(resnames_after_lig)) # sort| unique - sorted(set(list)) https://stackoverflow.com/questions/2931672/what-is-the-cleanest-way-to-do-a-sort-plus-uniq-on-a-python-list
    res_after_lig = {resnm:resnames_after_lig.count(resnm) for resnm in uniqres} # Form dictionary of resnm:nat
    #print ('Unique resnames after lig:', res_after_lig)
    if len(res_after_lig) == 0:
        logging.warning('Found no resnm after the ligand in the structure. This structure sequence ordering is uncommon for the \"parse-template\" SetupType.')
        #sys.exit(0)
    if len(res_after_lig) > 1:
        logging.warning('Found more than one resnm after the ligand in the structure. Make sure all these residues are included in tc-groups in MDP.')
    if len(res_after_lig) > 0:
        logging.info('Residues found in the structure after the ligand: %s',res_after_lig)

    # Insert ligand_lines
    for l in coor_lig:
        if re.search(r'' + re_pattrn, l): coor_new.append(l)
    # append after_lig lines
    for l in coor_after_lig:
        coor_new.append(l)
    return coor_new

def GetMolNm(MolType, Molecules, bChkUniq = False):
    molnm = []
    for mol in Molecules:
        if mol.get('type','') == MolType:
            molnm.append( mol.get('name','').strip() )
    if bChkUniq:
        if len(molnm) == 0:
            logging.error('Found no definition of molecule with the type \"{MolType}\" for the phase {sPhase} in the config')
            sys.exit(0)
        if len(molnm) > 1:
            logging.error('Found multiple molecules {molnm} for the phase {sPhase} that have type {MolType}')
            sys.exit(0)
    return molnm

def parse_template_coor(lignm, CoorLigandPath, CoorTemplatePath, outCoorPath,log):
    ##
    ## Set System Structure
    ##
    ligFileNm = os.path.basename(CoorLigandPath)
    log.info(f"Substituting ligand coordinates {ligFileNm} to the system structure template {CoorTemplatePath}")

    coor_ext = os.path.splitext(CoorTemplatePath)[1]   # Extract file extension: prefix, ext = os.path.splitext(path)
    coor_template   = read_file(CoorTemplatePath)  # Read structure of the whole system as a template
    coor_ligand     = read_file(CoorLigandPath)  # Read structure of the ligand to be parsed into template
    
    #print('Ligans resname: ', lignm)
    if coor_ext == '.pdb':
        coor = subst_lig_coor(coor_template,'^(?=ATOM|HETATM).{17}(.{3}).+\n$',lignm, coor_ligand,log)  # (?=ATOM|HETATM) - POSITIVE LOOKAHEAD Zero-Length Assertions https://www.regular-expressions.info/lookaround.html
    if coor_ext == '.gro':   # TODO: correct regex pattern. What if ext for lig and sys is different?
        coor = subst_lig_coor(coor_template,'^(?=ATOM|HETATM)To Be Done for GRO-format\n$',lignm, coor_ligand,log)  # (?=ATOM|HETATM) - POSITIVE LOOKAHEAD Zero-Length Assertions https://www.regular-expressions.info/lookaround.html
    #OutCoorPath = outCoorPrefix + coor_ext
    write_file(coor, outCoorPath)
    #return OutCoorPath

def  check_gmx(bingmx:str, log): #  call gmx without args
    """
    Call gmx without input and extract from output its version
    
    :param         bingmx: gmx binary
    :param            log: log channel
    :return:       versio: gmx version as printed in its headers
    """
    ## https://stackoverflow.com/questions/163542/how-do-i-pass-a-string-into-subprocess-popen-using-the-stdin-argument
    ## On Python 3.5+ (3.6+ for encoding), you could use subprocess.run, to pass input as a string to an external command
    ## and get its exit status, and its output as a string back in one call
    p = subprocess.run(bingmx, shell=True, capture_output=True)

    # Check if the structure file was generated
    sOut = ''.join( re.findall(r'([\n]GROMACS[^\n]+)[\n]', p.stderr.decode('utf-8') ) )
    version = ''.join( re.findall(r'[\n](GROMACS\:[ \t]+[^\n]+)[\n]', sOut ) )
    if len(sOut) > 0 :
        log.info(f'Found executable {version} called by command:\"{bingmx}\"')
    else:
        log.error(f"PROBLEM CALLING GMX EXECUTABLE \"{bingmx}\":\n"    \
                   + "#"*80                   + '\n'               \
                   + p.stdout.decode('utf-8') + '\n'               \
                   + p.stderr.decode('utf-8') + '\n'               \
                   + "#"*80 )
        sys.exit(0)

def make_ndx_fromTemplate(bingmx: str, template_ndx: str, structure: str, cmd_ndx:str):
    """
    Generate system NDX file from template
    
    :param   template_ndx: NDX template file full path
    :param      structure: The full system structure file
    :param        cmd_ndx: stdin command for "gmx make_ndx" routine; by defult cmd_ndx="r LIG\\nr SOL\\nq\\n"
    :return:
    """
    # cmd_ndx = 'r LIG\nr SOL\nq\n'  # Default value is good for most cases
    #print ('make_ndx_cmd:', cmd_ndx)
    cmd = bingmx + ' make_ndx -f ' + structure + ' -o index.ndx -n ' + template_ndx
    logging.info(f"Generating the system NDX-file from Template")
    logging.info(f"Executing command: {cmd}")
    logging.info(f"make_ndx stdin commands:\n{cmd_ndx}")
    ## https://stackoverflow.com/questions/163542/how-do-i-pass-a-string-into-subprocess-popen-using-the-stdin-argument
    ## On Python 3.5+ (3.6+ for encoding), you could use subprocess.run, to pass input as a string to an external command
    ## and get its exit status, and its output as a string back in one call
    p = subprocess.run(cmd, input=cmd_ndx.encode('ascii'), shell=True, capture_output=True)

    # Analize output of make_ndx - extract generated groups
    ndx_groups = ''
    for line in re.findall(r'[^\n][0-9]+[ \t]+[^ \t\:]+[ \t]*\:[ \t]*[0-9]+[ \t]+atoms[^\n]*\n', p.stdout.decode('utf-8')):
        ndx_groups = ndx_groups + line
    if os.path.isfile('./index.ndx') and len(ndx_groups) > 0:
        ndx_groups = 'Generated NDX Groups:\n' + ndx_groups
        logging.info(ndx_groups)
    else:
        logging.error(f"PROBLEM AT GENERATING INDEX FILE")
        logging.error(p.stdout.decode('utf-8'))
        logging.error(p.stderr.decode('utf-8'))
        sys.exit(0)

def  center_box(bingmx:str, box, inStruc:str,  outStruc:str, log): #  generate shifted to center ligand stucture by "gmx editconf"
    """
    Generate box-centered structure
    
    :param            box: box size (str with 3 float dimensions)
    :param        inStruc: Input structure file
    :param       outStruc: Output structure file
    :return:
    """
    #print ('make_ndx_cmd:', cmd_ndx)
    cmd = bingmx + ' editconf -princ -box ' + box + ' -f ' + inStruc + ' -o ' + outStruc
    log.info(f"Create the structure with the ligand at the ceneter of the box \"{box}\"")
    log.info(f"Executing command: {cmd}")
    ## https://stackoverflow.com/questions/163542/how-do-i-pass-a-string-into-subprocess-popen-using-the-stdin-argument
    ## On Python 3.5+ (3.6+ for encoding), you could use subprocess.run, to pass input as a string to an external command
    ## and get its exit status, and its output as a string back in one call
    p = subprocess.run(cmd, input='0'.encode('ascii'), shell=True, capture_output=True)

    # Check if the structure file was generated
    lsOut = ''.join( re.findall(r'([^:\n]+[:][^:\n]+[\n])', p.stdout.decode('utf-8') ) )
    if len(lsOut) > 10 and os.path.isfile(outStruc) and os.path.getsize(outStruc) > 0:
        log.info(f'The centered structure \"{outStruc}\" was generated.')
        log.debug('EDITCONF OUTPUT:\n%s', "#"*80 +'\n' + ''.join(lsOut) + "#"*80 +'\n')
    else:
        log.error(f"PROBLEM AT GENERATING CENTERED STRUCTURE FILE \"{outStruc}\":")
        log.error(p.stdout.decode('utf-8'))
        log.error(p.stderr.decode('utf-8'))
        sys.exit(0)

def  solvate_box(bingmx:str, box, inStruc:str, outStruc:str, Top:str, log): #  generate shifted to center ligand stucture by "gmx editconf"
    # For Water solvent use solvate without "-cs" option (unless set gmxpth_solvate=gmx5.1) because "-cs" results in seg-fault crashes with gmx2016+ ver.
    # -cs option is probably required for non-water solvents only
    #echo -e "${gmxpth_solvate}/gmx_d solvate -scale 0.45 -cp $MolNmOut'_centered.pdb' -box $box $box $box -o $outnm'_0.pdb' -p $outnm'.top'\nrm '#'$outnm'.top'*" >> $sge_script0

    """
    Generate solvated box structure
    
    :param            box: box size (str with 3 float dimensions)
    :param            Top: Input (also output) topology file
    :param        inStruc: Input structure file
    :param       outStruc: Output structure file
    :return:
    """
    #print ('make_ndx_cmd:', cmd_ndx)
    cmd = bingmx + ' solvate -scale 0.45 -box ' + box + ' -cp ' + inStruc + ' -o ' + outStruc + ' -p ' + Top
    log.info(f"Solvate the structure with the ligand at the ceneter of the box \"{box}\"")
    log.info(f"Executing command: {cmd}")
    ## https://stackoverflow.com/questions/163542/how-do-i-pass-a-string-into-subprocess-popen-using-the-stdin-argument
    ## On Python 3.5+ (3.6+ for encoding), you could use subprocess.run, to pass input as a string to an external command
    ## and get its exit status, and its output as a string back in one call
    p = subprocess.run(cmd, shell=True, capture_output=True)

    # Check if the structure file was generated
    lsOut = ''.join( re.findall(r'([^:\n]+[:][^:\n]+[\n])', p.stderr.decode('utf-8') ) )
    if len(lsOut) > 10 and os.path.isfile(outStruc) and os.path.getsize(outStruc) > 0:
        log.info(f'The centered structure \"{outStruc}\" was generated.')
        log.debug('EDITCONF OUTPUT:\n%s', "#"*80 +'\n' + ''.join(lsOut) + "#"*80 +'\n')
    else:
        log.error(f"PROBLEM AT GENERATING CENTERED STRUCTURE FILE \"{outStruc}\":")
        log.error(p.stdout.decode('utf-8'))
        log.error(p.stderr.decode('utf-8'))
        sys.exit(0)

def addTask_to_workflow(workflow, Task:str, coorFileNm:str, MdpFile:str, JobArgs, args):
##
## Each MD Task requires execution of "grompp" and "mdrun"
##
    #sRunNm = sLigNm + "_" + ShortPrefixDic.get(sPhase)
    #outnm  = args.gmxconfig.paths.get('outnm')
    bingmx = JobArgs.get('gmx')
    if utils.env.Task[Task].get('type') == 'minimization':
        bingmx = JobArgs.get('gmx_d')
    # bingmx = args.gmxconfig.paths.get('gmx')
    # if args.gmxconfig.Task.get(Task).get('type') == 'equilibration':
    #     bingmx = args.gmxconfig.paths.get('gmx_d')
    # grompp_flags = args.gmxconfig.Phase.get(sPhase).get('grompp_flags') # ' -maxwarn 2000' - in Protein to suppress an error due to multiple warrnings
    # if args.gmxconfig.Phase.get(sPhase).get('UseNdx'):
    #      ndx_option   = ' -n ${JobDir}/index.ndx'
    # if args.gmxconfig.Phase.get(sPhase).get('PosreStructure'):
    #     posre_option = ' -r ${JobDir}/{sRunNm}_Xray.pdb'
    outnm        = JobArgs.get('outnm')
    # sRunNm       = JobArgs.get('sRunNm')
    # grompp_flags = JobArgs.get('grompp-flags') if hasattr(JobArgs, 'grompp-flags')   else ''
    # ndx_option   = JobArgs.get('ndx-option')   if hasattr(JobArgs, 'ndx-option')   else ''
    # posre_option = JobArgs.get('posre-option') if hasattr(JobArgs, 'posre-option') else ''

    # $gmxpth/gmx_d grompp -maxwarn 2000 -o $runnm'.tpr' ${ndx_option} -f ${mdpfile} -p ../../$outnm'.top' -c ../../$iniconf  ${restr_option} >> ../all_tpr.out 2>&1
    # if [ ! -r $runnm'.tpr' ]; then
    #    echo ERROR: $runnm.tpr for Lamb=$iLamb was not created, exit the job.
    #    exit
    # fi
    #  echo -e "cd ${SUBDIR}" >>  $script
    #  tab='    '
    #  echo -e "for ((iLamb=$FirstTIPointID; iLamb <= $LastTIPointID; iLamb++)); do" >>  $script
    #     #
    #     #  Add creation of tpr-files for all lambdas into the sumsission script
    #     #
    #     echo -e "${tab}lamDirNm=\${LamDirPrefix}\$(printf '%02i' \$iLamb)"  >>  $script
    #     echo -e "${tab}cd \${lamDirNm}" >>  $script
    #     [ ${type} == "nvt" ] &&  echo -e "${tab}cp ../../MIN/\${lamDirNm}/${outnm}_min.gro ./" >>  $script
    #     [ ${type} == "npt" ] &&  echo -e "${tab}cp ../../NVT/\${lamDirNm}/${outnm}_nvt.gro ./" >>  $script
    #     [[ (${type} == "prod") && (${phase} != "InG") ]] &&  echo -e "${tab}cp ../../NPT/\${lamDirNm}/${outnm}_npt.gro ./" >>  $script     # g=1; c=123; [[ ( "$g" -eq 1 && "$c" = "123" ) || ( "$g" -eq 2 && "$c" = "456" ) ]] && echo "g = $g"
    #     [[ (${type} == "prod") && (${phase} == "InG") ]] &&  echo -e "${tab}cp ../../MIN/\${lamDirNm}/${outnm}_min.gro ./" >>  $script     # g=1; c=123; [[ ( "$g" -eq 1 && "$c" = "123" ) || ( "$g" -eq 2 && "$c" = "456" ) ]] && echo "g = $g"
    #
    #     echo -e "${tab}$gmxpth/gmx_d grompp -maxwarn 2000 -o $runnm'.tpr' ${ndx_option} -f ${mdpfile} -p ../../$outnm'.top' -c $iniconf  ${restr_option} >> ../all_tpr.out 2>&1"  >> $script
    #     #exit
    #     echo -e "${tab}if [ ! -r $runnm'.tpr' ]; then\n${tab}${tab}echo ERROR: $runnm.tpr  for Lamb=\$iLamb was not created, exit the job.\n${tab}${tab}exit\n${tab}fi"  >> $script
    #     #echo -e "${tab}if [ ! -r $runnm'.tpr' ]; then\n${tab}${tab}echo ERROR: $runnm.tpr  for Lamb=\$iLamb was not created, exit the job.\n${tab}${tab}mv ../$JobDir ../$JobDir'_FAILED'\n${tab}${tab}rm -r \$TmpDir\n${tab}${tab}exit\nfi\n${tab}rm mdout.mdp"  >> $script
    #     echo -e "${tab}cd ../" >>  $script
    #  echo -e "done # End of loop over Lamb"  >>  $script
    #fi
    #echo -e "" >>  $script
    workflow.append('\n##########################################')
    workflow.append(  '############## Task: ' + Task + ' #################')
    workflow.append(  '##########################################')
    workflow.append('cd ' + Task + '/${subdir}/')
    # grompp commands
    #print(JobArgs)
    line = bingmx + ' grompp ' + JobArgs.get('grompp-flags')       \
         + ' -f ' + MdpFile                                        \
         + ' -o ' + outnm + '.tpr'                                 \
         + ' -p ${JobDir}/' + JobArgs.get('sRunNm') + '.top'       \
         + ' -c ' + coorFileNm                                     \
         + JobArgs.get('ndx-option')                               \
         + JobArgs.get('posre-option')                             \
         + ' >> ${JobDir}/' + Task + '/tpr_out 2>&1'
    workflow.append(line)
    #workflow.append('    check_file '+ outnm +'.tpr \'ERROR: ' + outnm +'.tpr was not created for Task: '+Task+', exit the job.\'')
    workflow.append('if [ ! -r \'' + outnm + '.tpr\' ]; then')
    workflow.append('    echo ERROR: \'' + outnm +'.tpr\' for Task: '+Task+' was not created, exit the job.')
    workflow.append('    exit')
    workflow.append('fi')
    # mdrun command
    line = bingmx + ' mdrun ' +  JobArgs.get('mdrun-flags')         \
           + ' -nt ' + str(args.ncpu)                               \
           + ' -s ' + outnm + '.tpr'                                \
           + ' -dhdl dhdl.xvg'                                      \
           + ' -c ' + outnm +'_' + Task.lower() + '.gro'            \
           + ' -deffnm ' + outnm                                    \
           + ' -g ' + JobArgs.get('job_dir')+'/'+Task+'/${subdir}/'+outnm+'.log'
    workflow.append(line)
    workflow.append('cd ${JobDir}')
    #workflow.append('    move_files_to_jobdir \''+ JobArgs.get('sJobDirFull')+'/'+Task+'/\'')
    return workflow


def job_init(sJobNm, sLigNm, sPhase, args):
    # Check that Phase and ligand files are defined
    phase_cfg = utils.env.Phase.get(sPhase)
    #phase_cfg = args.gmxconfig.Phase.get(sPhase)
    if len(phase_cfg) == 0:
        logging.error(f"The phase \"{sPhase}\" is not defined in the config {args.config}")
        sys.exit(0)
    top = utils.env.TopDir + '/' + sLigNm + '.itp'
    pdb = utils.env.TopDir + '/' + sLigNm + '.pdb' 
    if not (os.path.isfile(top) and os.path.getsize(top) > 0):  # Check that ligand TOP exists
        logging.error(f"The ligand topology does not exist: {top}")
        sys.exit(0)
    if not (os.path.isfile(pdb) and os.path.getsize(pdb) > 0):  # Check that ligand TOP exists
        logging.error(f"The ligand structure does not exist: {pdb}")
        sys.exit(0)

    # Initialize Job arguments
    bingmx   = utils.env.sGmx
    bingmx_d = utils.env.sGmx_d
    #bingmx   = args.gmxconfig.paths.get('gmx')
    #bingmx_d = args.gmxconfig.paths.get('gmx_d')
    #if not bingmx:
    #    logging.error("gmx binary path is not defined in config.")
    #    sys.exit(0)
    #if not bingmx_d:
    #    logging.warning("\"gmx_d\" double precision binary path is not defined in config. Will use the single precision \"gmx\" for minimization (not recommended)")
    #    bingmx_d = bingmx
    #check_gmx(bingmx  , logging)
    #check_gmx(bingmx_d, logging)
    # if not (os.path.isfile(bingmx) and os.path.getsize(bingmx) > 0):  # Check that ligand TOP exists
    #     logging.error(f"GMX binary does not exist: \"{bingmx}\" ")
    #     sys.exit(0)
    logging.info(f"gmx   binary: {bingmx}")
    logging.info(f"gmx_d binary: {bingmx_d}")
    
    FileLabel = phase_cfg.get('FileLabel','')
    DirPrefix = phase_cfg.get('DirPrefix','')
    if FileLabel == '':
        logging.error(f"The \"FileLabel\" parameter for the Phase {sPhase} is not defined in config. Exit!")
        sys.exit(0)
    if DirPrefix == '':
        logging.error(f"The \"DirPrefix\" parameter for the Phase {sPhase} is not defined in config. Exit!")
        sys.exit(0)

    sRunNm = sLigNm + "_" + FileLabel
    # JobArgs['sJobDirFull'] = job_dir
    # JobArgs['sRunNm'] = sRunNm
    job_id = str(1 + scheduler.get_last_job_id(utils.env.queue))   # TODO rename vars at sge 
    job_dir_name = job_id + '.' + DirPrefix + '_' + sJobNm
    logging.debug(f"job dir name: {job_dir_name}")
    # create job directory
    job_dir = os.path.abspath(utils.env.OutDir + "/" + job_dir_name)
    logging.info("path to job dir: %s", job_dir)
    #  Create Working Folder
    os.makedirs(job_dir)
    if utils.env.Type == 'run_ti':
        logging.info(f"Starting TI simulation of {sJobNm} in {sPhase}")
    elif utils.env.Type == 'run_md':
        logging.info(f"Starting MD simulation of {sJobNm} in {sPhase}")
    else:
        logging.info(f"Starting Mutation={sJobNm} Phase={sPhase}")
    logging.info(f"JobDir: {job_dir}")

    # Initialize Dictionary of Job Script parameters
    JobArgs = {'gmx'             : bingmx,
               'gmx_d'           : bingmx_d,
               'outnm'           : utils.env.sOutnm,
               'sLigNm'          : sLigNm,
               'sPhase'          : sPhase,
               'sJobNm'          : sJobNm,
               'sRunNm'          : sRunNm,
               'job_dir'         : job_dir,
               'jobid'           : job_id,
               'grompp-flags'    : phase_cfg.get('grompp_flags',''),    # ' -maxwarn 2000' - in Protein to suppress an error due to multiple warrnings
               'mdrun-flags'     : phase_cfg.get('mdrun_flags',''),     # ' -ntmpi' - To limit number of MPI-threads per mdrun
               'ndx-option'      : '',
               'posre-option'    : '',
               'workflow'        : []
               }
    return JobArgs

def top_coor(sLigNm, sPhase, JobArgs, args):
    """
    MAIN ROUTINE to setup MD input Topology and Coordinate files
    :param sJobNm: name of the job
    :param sLigNm: ligand name
    :param sPhase: phase in template (Gas, Water, Protein)
    :param args: arguments of gmxFE command
    :return:
    """
    ##
    ## SETUP SIMULATION BOX
    ##
    ##
    ## 1. SETUP config parameters
    ## Define GMX binary paths
    sRunNm = JobArgs['sRunNm']
    job_dir = JobArgs['job_dir']
    bingmx = JobArgs['gmx']
    TopPath = sRunNm + '.top'
    CoorPath = sRunNm + '.pdb'
    #phase_cfg = args.gmxconfig.Phase.get(sPhase)
    phase_cfg = utils.env.Phase.get(sPhase)
    ##
    ## Set System Topology
    ##
    copy_parse_topology(sLigNm,TopPath,job_dir,sPhase,args,logging)
    ##
    ## Set System Structure
    ##
    SetupType = phase_cfg.get('SetupType')
    if SetupType: logging.info(f"Do \"{SetupType}\" for setting up the \"{sPhase}\" phase system BOX.")
    else:
        logging.error(f'SetupType for phase {sPhase} was not found in the config.')
        sys.exit(0)
    if SetupType == 'parse-template':
        lignm = GetMolNm('ligand', phase_cfg.get('Molecule'), True)[0]
        CoorTemplatePath = os.path.abspath(utils.env.TopDir + '/' + phase_cfg.get('IniStructure'))
        CoorLigandPath   = os.path.abspath(utils.env.TopDir + '/' + sLigNm + '.pdb')
        parse_template_coor(lignm, CoorLigandPath, CoorTemplatePath, CoorPath,logging)
        ##
        ## Generate Ndx File from Template
        ##
              # cp $indir/index.ndx index.ndx
              # outndx=$(echo -e "r LIG\nr SOL\nq\n"| $gmxpth/gmx_d make_ndx -n index -f $outnm'_ini.pdb' 2>/dev/null| sed -n '/LIG\|SOL/p')
              # echo -e "Make index file for a given mutation of ligands:\n$outndx"
        if phase_cfg.get('UseNdx'):
            JobArgs['ndx-option'] = ' -n ${JobDir}/index.ndx'
            NdxTemplatePath =  os.path.abspath(utils.env.TopDir + '/' + phase_cfg.get('NdxTemplate'))
            logging.info(f"Generating the system NDX-file from Template: {NdxTemplatePath}")
            make_ndx_cmd = phase_cfg.get('Make_ndx_CMD', 'r LIG\\nr SOL\\nq\\n') # Read User defined or set default make_ndx_cmd = 'r LIG\nr SOL\nq\n', NOTE '\\n' because https://stackoverflow.com/questions/38401450/n-in-strings-not-working
            #make_ndx_cmd = make_ndx_cmd if len(make_ndx_cmd) else 'r LIG\nr SOL\nq\n' # set default make_ndx_cmd = 'r LIG\nr SOL\nq\n'
            # All resnm after ligand entry to be generated by make_ndx
            coor   = read_file(CoorPath)
            lignm = GetMolNm('ligand', phase_cfg.get('Molecule'), True)[0]
            res_after_lig = resnm_after_lig(coor, lignm, '^(?=ATOM|HETATM).{17}(.{3}).+\n$')  # PDB: resnm position 18-20  # (?=ATOM|HETATM) - POSITIVE LOOKAHEAD Zero-Length Assertions https://www.regular-expressions.info/lookaround.html
            if len(res_after_lig) > 0:
                logging.info('Merging the make_ndx command \'%s\' with residues found in the structure after the ligand: %s',make_ndx_cmd,list(res_after_lig.keys()))
                for resnm in res_after_lig.keys():
                    #print('re.search(r\'r '+resnm+'\\\\n\'):',re.search(r'r '+resnm+'\\\\n', make_ndx_cmd))
                    if not re.search(r'r[ \t]+'+resnm+'\\\\n', make_ndx_cmd): make_ndx_cmd = 'r '+resnm+'\\n'+make_ndx_cmd  # '\\n' because https://stackoverflow.com/questions/38401450/n-in-strings-not-working
            logging.info('The command passed to the make_ndx routine: \'%s\'',make_ndx_cmd)
            # NOTE, before this line all '\n' chars in make_ndx_cmd should be actually 2-char strings '\\n'. Bellow we convert to '\n' new line symb.
            make_ndx_cmd = re.sub(r'\\n','\n',make_ndx_cmd) # Subs '\\n' by the new line symbol '\n':
            make_ndx_fromTemplate(bingmx, NdxTemplatePath, CoorPath, make_ndx_cmd) #  generate index goups for "LIG" and "SOL" by "gmx make_ndx"
        #sys.exit(0)
    elif SetupType == 'box-center+solvate':
        logging.info(f"Place the ligand \"{sLigNm}\" at the center of the simulaion box and solvate")
        CoorLigandPath  =  os.path.abspath(utils.env.TopDir + '/' + sLigNm + '.pdb')
        box = phase_cfg.get('boxsize')
        center_box(bingmx, box, CoorLigandPath,  sRunNm + '_centered.pdb', logging) #  generate shifted to center ligand stucture by "gmx editconf"
        solvate_box(bingmx, box, sRunNm + '_centered.pdb',  sRunNm + '.pdb', sRunNm + '.top', logging) #  solvate the ligand stucture by "gmx solvate" and add solvent to topology
        #sys.exit(0)
    elif SetupType == 'box-center':
        logging.info(f"Place the ligand \"{sLigNm}\" at the center of the simulaion box.")
        CoorLigandPath  =  os.path.abspath(utils.env.TopDir + '/' + sLigNm + '.pdb')
        box = phase_cfg.get('boxsize')
        center_box(bingmx, box, CoorLigandPath,  sRunNm + '.pdb', logging) #  generate shifted to center ligand stucture by "gmx editconf"
        #sys.exit(0)
    elif SetupType == 'none':
        #CoorLigandPath  =  os.path.abspath(utils.env.TopDir + '/' + sLigNm + '.pdb')
        logging.info(f"Use the ligand initial structure \"{sLigNm}\"")
        copy_file_loc_2dist(sLigNm + '.pdb', utils.env.TopDir, job_dir + '/' + sRunNm + '.pdb')
    else:
        logging.error(f"The \"{SetupType}\" SetupType is not supported yet.")
        sys.exit(0)

    # Set position restrain reference file for ALL phases. But the use of POSRE is controlled by keywords defined in MDP
    PosreTemplate = phase_cfg.get('PosreStructure')
    if PosreTemplate:
        logging.info(f"Use the user defined structure for the POSRE reference coordinates: {PosreTemplate}")
        PosreTemplate = os.path.abspath(utils.env.TopDir + '/' + PosreTemplate)
        PosrePath = sRunNm + '_Xray.pdb'
        parse_template_coor(lignm, CoorLigandPath, PosreTemplate, PosrePath,logging)
        JobArgs['posre-option'] = ' -r ${JobDir}/' + PosrePath
    else:
       logging.info(f"Use the initial system structure for the POSRE reference coordinates")
       JobArgs['posre-option'] = ' -r ${JobDir}/' + CoorPath
    logging.info(f"Topology and Structure files for the ligand \"{sLigNm}\" in the \"{sPhase}\" phase are successfully generated.")
    #sys.exit(0)
    return JobArgs



def mdp_Tasks(sLigNm, sPhase, job_dir, JobArgs, args):
    ##
    ## Set MDP config files
    ##
    # Pick correct config template
    logging.info(f"\nSetting up the MD Parameter (MDP) files for all Tasks")
    phase_cfg = utils.env.Phase.get(sPhase)
    if utils.env.Type == 'run_ti':
        MdpTemplateName = phase_cfg.get('MdpTemplate')
        #MdpTemplateName = args.gmxconfig.Phase[sPhase].get('MdpTemplate')
        #MdpTemplateName = MdpTemplateDic.get(sPhase)
        if args.replex:
            MdpTemplateName = phase_cfg.get('MdpTemplate')
            #MdpTemplateName = args.gmxconfig.Phase[sPhase].get('MdpTemplate')
            #MdpTemplateName = MdpTemplateDic.get(sPhase)
            #MdpTemplatePath = os.path.abspath(utils.env.MdpDir + '/' + MdpTemplateDic.get(sPhase + '_HREX'))
    else:
        MdpTemplateName = phase_cfg.get('MdpTemplate')
        #MdpTemplateName = args.gmxconfig.Phase[sPhase].get('MdpTemplate')
        #MdpTemplateName = MdpTemplateDic.get('Solvent')
        if args.replex:
            MdpTemplateName = phase_cfg.get('MdpTemplate')
            #MdpTemplateName = args.gmxconfig.Phase[sPhase].get('MdpTemplate')
            #MdpTemplateName = MdpTemplateDic.get('Solvent')
            #MdpTemplatePath = os.path.abspath(utils.env.MdpDir + '/' + MdpTemplateDic.get('Solvent_HREX'))
    MdpTemplatePath = os.path.abspath(utils.env.MdpDir + '/' + MdpTemplateName)
    if utils.env.Template:
       custom_template = os.path.abspath(utils.env.MdpDir + '/' + utils.env.sMdpTemplate)
       if os.path.isfile(custom_template):
          MdpTemplatePath = custom_template
       else:
          logging.warning('User defined mdp-template %s not found. Using template from config', custom_template)
    logging.info("Making mdp files from Template: %s", MdpTemplatePath)
    mdp_template = read_file(MdpTemplatePath)  # Read mdp template

    ##
    ##   SETUP WORKFLOW
    ##
    ## 1. Create separate folder and MDP file for each of the workflow Task: "MIN", "NVT", "NPT", "PROD"
    ## 2. Add execution commands to the JOB script for each of the workflow Task
    #print(args.gmxconfig)
    Task_cfg = utils.env.Task if hasattr(utils.env, 'Task') else ''
    mdp_template = setMdpParams( Task_cfg.get('DEFAULT','').get('mdp','') , mdp_template)   # Set config mdp parameters common for all
    mdp_template = setMdpParams(phase_cfg.get('mdp','')     , mdp_template)   # Set config mdp parameters specific for a given Phase
    prefix = phase_cfg.get('FileLabel')
    CoorPath = JobArgs['sRunNm'] + '.pdb'
    TaskIniCoor =  '${JobDir}/' +  CoorPath
    for Task in phase_cfg.get('workflow',''):
    #for Task in args.gmxconfig.Phase[sPhase].get('workflow',''):
        #print('Task:', Task)
        #print('Phase.mdp:', Task, args.gmxconfig.mdp.get(Task))
        mdp_lines = setMdpParams(Task_cfg.get(Task,'').get('mdp',''), mdp_template)   # Set config mdp parameters specific for a given Stage
        #mdp_lines = setMdpParams(args.gmxconfig.mdp.get(Task,''), mdp_template)   # Set config mdp parameters specific for a given Stage
        #print(f'args.gmxconfig.[\'Phase\'][{sPhase}][\'mdp\'].get(\'{Task}\',\'\')]: ', args.gmxconfig['Phase'][sPhase]['mdp'].get(Task,''))
        mdp_lines = setMdpParams(phase_cfg.get('mdp','').get(Task,''), mdp_lines)   # Set config mdp parameters specific for a given Phase and Task
        #mdp_lines = setMdpParams(args.gmxconfig['Phase'][sPhase]['mdp'].get(Task,''), mdp_lines)   # Set config mdp parameters specific for a given Phase and Task
        os.makedirs(Task)
        #utils.makedir(job_dir + '/' + Task)
        # Save MDP for a given Task
        MdpFile = prefix + "_" + Task.lower() + ".mdp"
        write_file(mdp_lines, Task + "/" + MdpFile)

        JobArgs['workflow'].append(Task)
        JobArgs[Task] = [MdpFile, TaskIniCoor]
        #workflow_lines = addTask_To_Workflow(workflow_lines, Task, TaskIniCoor, MdpFile, JobArgs, args)

        TaskIniCoor = '${JobDir}/' + Task + '/${subdir}/' + utils.env.sOutnm + '_' + Task.lower() + '.gro'
        #TaskIniCoor = '${JobDir}/' + Task + '/${subdir}/' + args.gmxconfig.paths.get('outnm') + '_' + Task.lower() + '.gro'
    #for l in mdp_lines:
    #    print (l)
    #for l in workflow_lines:
    #    print (l)
    #sys.exit(0)
    return JobArgs

def Tasks_to_workflow_script(JobArgs, args):
#def Add_Tasks_Workflow_Script(JobArgs, args):
    ##
    ## JOB Workflow script file
    ##
    logging.info(f"Generating the JOB Workflow script")
    CoorPath = JobArgs['sRunNm'] + '.pdb'
    TaskIniCoor =  '${JobDir}/' +  CoorPath
    workflow = []
    for Task in JobArgs['workflow']:
        MdpFile     = JobArgs[Task][0]   # JobArgs[Task] = [MdpFile, TaskIniCoor] should be set in "Set_MDP_Tasks"
        TaskIniCoor = JobArgs[Task][1]    
        workflow = addTask_to_workflow(workflow, Task, TaskIniCoor, MdpFile, JobArgs, args)

        #TaskIniCoor = '${JobDir}/' + Task + '/${subdir}/' + args.gmxconfig.paths.get('outnm') + '_' + Task.lower() + '.gro'
    #for l in workflow_lines:
    #    print (l)
    #sys.exit(0)
    return workflow


def get_lambda_list(args):
    # Extract 'LambdaValues' from config (NOT TEMPLATE)
    # Mandotory parameter for 'run_ti'
    Task_cfg = utils.env.Task # The existence of 'Task' defition had been already checked in the Init
    vdw_lambdas  = Task_cfg['DEFAULT']['mdp'].get('vdw-lambdas')
    fep_lambdas  = Task_cfg['DEFAULT']['mdp'].get('fep-lambdas')
    # mdp_cfg = utils.env.mdp if hasattr(utils.env, 'mdp') else ''
    # vdw_lambdas  = mdp_cfg.get('general','').get('vdw-lambdas')
    # fep_lambdas  = mdp_cfg.get('general','').get('fep-lambdas')
    if not vdw_lambdas and not fep_lambdas:
        logging.error(f"The mandotory lambda grid \"fep-lambdas\" or \"vdw-lambdas\" is not defined in the \"mdp.general\" section of the config {cfgPath}. Exit")
        sys.exit(0)
    slFep_Lambdas = fep_lambdas.split() if fep_lambdas else vdw_lambdas.split()
    #slCoul_Lambdas = mdp_cfg['general'].get('coul-lambdas').split()
    #slVdw_Lambdas  = mdp_cfg['general'].get('vdw-lambdas').split()
    #slCoul_Lambdas = args.gmxconfig.mdp['general'].get('coul-lambdas').split()
    #slVdw_Lambdas  = args.gmxconfig.mdp['general'].get('vdw-lambdas').split()
    #logging.debug(f'sCoul_Lambdas: \"{slCoul_Lambdas}\"\nsVdw_Lambdas: \"{slVdw_Lambdas}\"')
    logging.debug(f'Fep_Lambdas: \"{slFep_Lambdas}\"')

    #logging.debug("Current Job Dir: %s", job_dir)
    idFirstLam = 0
    idLastLam  = len(slFep_Lambdas) - 1
    #idLastLam  = len(slVdw_Lambdas) - 1
    #logging.debug(f'idLastTIpoint: {idLastTIpoint}')
    logging.info("Running lambda set: %s", utils.env.RunTIstates)
    if "all" in utils.env.RunTIstates:
        Lambdas = list(range(idFirstLam,idLastLam+1))
    elif "edge" in utils.env.RunTIstates:
        Lambdas = [idFirstLam, idLastLam]
    else:
        if type(utils.env.RunTIstates) == list:
            # array comes from command line option
            Lambdas = list(map(int, utils.env.RunTIstates))
        else:
            # array comes from config as a string
            Lambdas = [int(x) for x in utils.env.RunTIstates.split(' ')] 
    logging.debug("Running lambdas: %s", Lambdas)
    return Lambdas

def TI_grid(JobArgs, args):
    # Create subdir and replicate MDP-file for each TI point and each Task.
    Lambdas = get_lambda_list(args)
    JobArgs['RunLambdas'] = Lambdas
    lam_prefix = utils.env.sTIdir
    #lam_prefix = args.gmxconfig.paths.get('TIdir','')
    for Task in JobArgs['workflow']:
        MdpFile = JobArgs[Task][0]
        mdp_template = read_file(Task + '/' + JobArgs[Task][0])  # Read mdp template for a given Task
        for iL in Lambdas:
            subdir = Task + '/' + lam_prefix + str(iL)
            os.makedirs(subdir)
            mdp_lines = setMdpParams({'init-lambda-state':iL}, mdp_template)   # Set config mdp parameters common for all
            write_file(mdp_lines, Task + '/' + lam_prefix + str(iL) + "/" + MdpFile)
    return JobArgs

def gen_md_submit_script(job_id, args):
    """
    Generate start script for MD job

    :param job_id:
    :param args:
    :return: script name
    """
    rerun = True if hasattr(args, 'rerun') and args.rerun else False
    mpi, run_local = False, False
    if (hasattr(args, 'mpi') and args.mpi) or (hasattr(args, 'replex') and args.replex):
        mpi = True
    if hasattr(args, 'queue') and args.queue == 'none':
        run_local = True
    run_script = 'jobscript.qsh'
    if rerun: run_script = 'rerun.qsh'
    # TODO PBS DIRECTIVES
    line_opt = f'#!/bin/sh\n' \
               f'# SGE DIRECTIVES\n' \
               f'#$ -N {args.jobname}\n' \
               f'#$ -S /bin/sh\n' \
               f'#$ -cwd\n' \
               f'#$ -j y\n' \
               f'#$ -V\n' \
               f'# SLURM DIRECTIVES\n' \
               f'#SBATCH -J {args.jobname}\n' \
               f'#SBATCH -o {args.jobname}.o%j\n\n'
        # cmd = f'python3 -u bin/mdrun.py ' \ ...""
    if mpi: line_cmd += ' --mpi'
    if rerun: line_cmd += ' --rerun'
    if run_local: line_cmd += ' --local'
    with open(run_script, 'w') as f:
        f.write(line_opt)
        f.write(line_cmd)
    os.chmod(run_script, stat.S_IRWXU)
    return run_script

def gen_ti_submit_script(wf_lines, JobArgs, lam_ids_list, args):
    """
    Generate submit scripts for TI points

    :param JobArgs:
    :param lam_ids_list: list of lambda numbers
    :param args:
    :return: array of script names
    """
    import stat
    submit_scripts = []
    dir_prefix = utils.env.sTIdir
    mpi       = True if (hasattr(args, 'mpi') and args.mpi)  or (hasattr(args, 'replex') and args.replex) else False
    run_local = True if hasattr(args, 'queue') and args.queue == 'none' else False
    rerun     = True if hasattr(args, 'rerun') and args.rerun else False
    for idTIp in lam_ids_list:
        job_script = 'runjob' + str(idTIp) + '.qsh'
        if rerun: job_script = 'rerunjob' + str(idTIp) + '.qsh'
        # TODO PBS DIRECTIVES
        line_opt = f'#!/bin/sh\n' \
                f'# SGE DIRECTIVES\n' \
                f'#$ -N {args.jobname}{str(idTIp)}\n' \
                f'#$ -S /bin/sh\n' \
                f'#$ -cwd\n' \
                f'#$ -j y\n' \
                f'#$ -V\n' \
                f'# SLURM DIRECTIVES\n' \
                f'#SBATCH -J {args.jobname}{str(idTIp)}\n' \
                f'#SBATCH -o {args.jobname}{str(idTIp)}.o%j\n\n'

        jobnm        = JobArgs['sJobNm']
        LastTaskDir  = JobArgs['workflow'][-1]
        ResFile = utils.env.ResFile
        ## NOTE: OSX 'sed -i' requires extra param to rename old file, see https://stackoverflow.com/questions/5694228/sed-in-place-flag-that-works-both-on-mac-bsd-and-linux
        ## NOTE: sh shall is OLD and limited compare to bash and has different quotation rules
        line_record = '############# WRITE JOB RECORD #############\n'\
                      '####### removing duplicating records #######\n'\
                      '############################################\n'\
                      f'ResFile=\'{ResFile}\'\n' \
                      f'outdir=$(pwd)/{LastTaskDir}\n' \
                      f'outdir_esc=$(echo \"${{outdir}}\" | sed \'s/\([^a-zA-Z0-9]\)/\\\\\\1/g\')\n' \
                      f'sed -i.bak \'/\'\"${{outdir_esc}}\"\'/d\' ${{ResFile}}\n' \
                      f'rm ${{ResFile}}.bak\n' \
                      f'printf \"{jobnm}\\t$(date +%d.%m.%y-%T )\\t${{outdir}}\\n\" >> ${{ResFile}}\n'
        with open(job_script, 'w') as f:
            f.write(line_opt)
            f.write('subdir=' + dir_prefix + str(idTIp))
            f.write('\n'.join(wf_lines) + '\n')
            f.write(line_record)
            #f.write(line_cmd)
        os.chmod(job_script, stat.S_IRWXU)
        submit_scripts.append(job_script)
    return submit_scripts

def submit_scripts(workflow, JobArgs, args):
        ##
        ##  Generate the JOB script files for each lambda (in case plain TI)
        ##
    job_id = JobArgs['jobid']
    #logging.info(f'Handling the simulation type: \"{utils.env.Type}"')
    if utils.env.Type == 'run_ti':
        # Specify the structure of lambda foders for TI simulations: JobDir = '../..'
        workflow = [re.sub(r'' + re.escape('${JobDir}'), '../..', line)  for line in workflow]
        for l in workflow:
            print (l)
        Lambdas = JobArgs['RunLambdas']
        ##  Generate the JOB script files for each lambda (in case plain TI)
        if args.replex:
            # Generate run script
            script = gen_md_submit_script(job_id, args)
            return script
        else:
            # Generate run script for each TI point
            scripts = gen_ti_submit_script(workflow, JobArgs, Lambdas, args)
            return scripts
    elif utils.env.Type == 'run_md':
        # Specify the structure of lambda foders for TI simulations: JobDir = '../..'
        workflow = [re.sub(r'' + re.escape('${JobDir}'), '..', line)  for line in workflow]
        for l in workflow:
            print (l)
        # PLAIN MD
        script = gen_md_submit_script(workflow, job_id, Lambdas, args)
        return script
    else:
        logging.error(f'The simulation type \"{utils.env.Type}" is not implemented')
        sys.exit(0)

def schedule_job(jobscript, JobArgs, args):
    if utils.env.Type == 'run_ti':  # If TI
        Lambdas = JobArgs['RunLambdas']
        nThreads = len(Lambdas)
        if args.replex > 1:
            nThreads = args.replex
            if nThreads > len(Lambdas):
                logging.error("The number of user defined MPI threads nthreads = %3d "
                              "must be in the range between 1 and len(Lambdas) = %3d.",nThreads, len(Lambdas))
                sys.exit(0)
        mpiCores = nThreads * args.ncpu
        if args.replex:
            script = jobscript
            # Submit job
            nTIpPerThread = int(len(Lambdas) / nThreads) + ((len(Lambdas) % nThreads) > 0)
            if nTIpPerThread == 1:
                logging.info("Starting simulation with replica exchange for all lambda points in parallel. "
                             "Works only in Open MPI environment.")
            else:
                logging.info("Starting lambda replica exchange simulation with %d parallel MPI threads and %d lambda points "
                             "per thread. Works only in Open MPI environment.",nThreads,nTIpPerThread)
            scheduler.submit_job(script, args.jobname,
                                 args.partition, args.queue,
                                 args.ncpu, args.gpu, mpiCores)
        else:
            # Submit separate job for each TI point
            for idTIp, script in zip(Lambdas, jobscript):
                logging.info("Starting simulation for TI point: %i",idTIp)
                if args.jobname == 'gmxTIp':
                    jobname = 'gmxTIp' + str(idTIp)
                else:
                    jobname = args.jobname
                scheduler.submit_job(script, jobname,
                                 args.partition, args.queue,
                                 args.ncpu, args.gpu)
                #sys.exit(0)
                time.sleep(0)
    elif utils.env.Type == 'run_md':  # If MD
    # Plain MD single thread job
        script = jobscript
        logging.info("Starting plain MD simulation job")
        jobname = args.jobname
        scheduler.submit_job(script, jobname,
                         args.partition, args.queue,
                         args.ncpu, args.gpu)
    else:
        logging.error(f'The simulation type \"{utils.env.Type}" is not implemented')
        sys.exit(0)
