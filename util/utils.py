#
#
# Copyright (C) 2020 Niko Popitsch.  All rights reserved.
#
# This file is part of Nanopanel2, a somatic variant caller for Nanopore sequencing data.
# 
# Commercial License Usage
# 
# Licensees holding valid commercial nanopanel2 licenses may use this file in
# accordance with the commercial license agreement provided with the
# Software. For licensing terms and conditions see the file LICENSE.
# 
# Alternatively, you may use this file under the terms of the 
# GNU Affero General Public License v3.0 license as set out in the LICENSE file
#
#     This program is free software: you can redistribute it and/or modify
#     it under the terms of the GNU Affero General Public License as published
#     by the Free Software Foundation, either version 3 of the License, or
#     (at your option) any later version.
# 
#     This program is distributed in the hope that it will be useful,
#     but WITHOUT ANY WARRANTY; without even the implied warranty of
#     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#     GNU Affero General Public License for more details.
# 
#     You should have received a copy of the GNU Affero General Public License
#     along with this program.  If not, see <https://www.gnu.org/licenses/>
#    
# @author: niko.popitsch

'''
    Nanopanel2 utility functions
'''

if __name__ == '__main__':
    pass

from argparse import ArgumentTypeError
import logging, subprocess, os, sys
from subprocess import check_output, CalledProcessError


def run_task(cmd, shell=False):
    """ Run a commandline task """
    logging.info(cmd)
    if shell:
        out = check_output(" ".join(cmd), shell=True, stderr=subprocess.STDOUT)
    else: 
        out = check_output(cmd, stderr=subprocess.STDOUT)
    return out
        
def files_exist(files):
    """ Check whether a (list of) file(s) exists """
    if (type(files) is list) :
        for f in files:
            if f is None:
                return False
            if not os.path.exists(f):
                return False
    else:
        if files is None:
            return False
        if not os.path.exists(files):
            return False
    return True


def existing_file(files):
    """ For usage in argument parser """
    if files_exist(files):
        return files
    if (type(files) is list) :
        raise ArgumentTypeError("Not all files exist ["+",".join(files)+"]")
    else:
        raise ArgumentTypeError("Not all files exist ["+(files)+"]")
    
def remove_file(files):
    """ remove a (list of) file(s) (if it/they exists) """
    if (type(files) is list) :
        for f in files:
            if os.path.exists(f):
                os.remove(f)
    else:
        if os.path.exists(files):
            os.remove(files)
   
def checkFile(files):
    """ Check whether a file exists and exit if not """
    if (type(files) is list) :
        for f in files:
            checkFile(f)
    else :
        if not files is None and not os.path.exists(files):
            print("Error: file", files, "was not found! Exiting...")
            sys.exit(1)       
                
def pipeline_step(inputfile, outFile, cmd, shell=False, stdout=None, append=False, logfile=None):
    """ Execute a pipeline step """
    try:
        if inputfile is not None:
            if (type(inputfile) is list):
                for f in inputfile:
                    checkFile(f)
            else:
                checkFile(inputfile) 
        
        if stdout is not None:
            if shell is False:
                raise ArgumentTypeError("When using the parameter stdout, the shell parameter must be set to True!")
            else: 
                if append:
                    cmd.append(">>")
                else:
                    cmd.append(">")
                cmd.append(stdout)                                                                                                                                                                                                                                                                 

        # Overwrite output file?!
        if outFile is not None:
            if (type(outFile) is list) :
                for f in outFile:
                    if files_exist(f):
                        logging.warn("Overwriting file %s", f)
                else:
                    if files_exist(outFile):
                        logging.warn("Overwriting file %s", outFile)

        out = run_task(cmd, shell)         
        if logfile is None:
            logging.info(out)
        else:
            with open(logfile, "a") as log:
                log.write(out)
        
        if outFile is not None:
            checkFile(outFile)
        return True
    except CalledProcessError as e:
        logging.error(e.output)
        logging.error("ERROR %s - removing outputfile %s", e, outFile)
        if outFile is not None:
            if (type(outFile) is list) :
                for f in outFile:
                    remove_file([f]) 
            else:
                remove_file([outFile]) 
        return False        

def index_bam(inFileBam, override=False, exe="samtools"):
    """ Index BAM file with samtoiols """
    success = True
    idxFile = inFileBam + ".bai"
    if(not files_exist(idxFile) or override):
        success = success and pipeline_step(inFileBam, idxFile, [exe, "index", inFileBam], shell=True)
    return success

def sam2bam(inFile, outFile, index=True, sort=True, delinFile=False, override=False, onlyUnique=False, onlyProperPaired=False, onlyPrimary=False, filterMQ=0, maxmem=None, exe="samtools"):
    """ Convert SAM to BAM using samtools """
    if(onlyUnique and filterMQ == 0):
        filterMQ = 1;
        
    success = True    
    if(not files_exist(outFile) or override):        
        cmd = [exe, "view", "-Sb", inFile, "-o", outFile]
        if filterMQ > 0:
            cmd += ["-q", str(filterMQ)]
        if onlyProperPaired:
            cmd += ["-f", "2"]
        if onlyPrimary:
            cmd += ["-F", "256"]
        success = success and pipeline_step(inFile, outFile, cmd, shell=True)
        
        if(sort):         
            tmp = outFile + "_tmp"
            os.rename(outFile, tmp)   
            cmd = [exe, "sort"]
            if not maxmem is None:
                cmd += ["-m", maxmem ]
            cmd+=["-o", outFile]    
            cmd += [tmp]             
            success = success and pipeline_step(tmp, outFile, cmd, shell=True)
            if(success):
                remove_file(tmp)
        if(index):
            success = success and index_bam(outFile, exe=exe)
            
        if(success and delinFile):
            remove_file(inFile)
    else:
        logging.info("Skipping sam2bam. " + outFile + " already exists.");
    return success

def sort_bam(inFile, outFile, index=False, override=False, delinFile=False, maxmem=None):
    """ sort BAM file using samtools """
    success = True    
    if(not files_exist(outFile) or override):     
        cmd = ["samtools", "sort", "-o", outFile]
        if not maxmem is None:
            cmd += ["-m", maxmem ]
        cmd += [ inFile ]             
        success = success and pipeline_step(inFile, outFile, cmd, shell=True)
        if(index):
            success = success and index_bam(outFile)          
        if(success and delinFile):
            remove_file(inFile)
    else:
        logging.info("Skipping sortbam. " + outFile + " already exists.");
    return success


def bgzip(inFile, outFile=None, override=False, delinFile=False, exe="bgzip"):
    """ bgzip a file """
    if outFile == None:
        outFile = inFile+".gz"
    success = True    
    if(not files_exist(outFile) or override):                        
        success = success and pipeline_step(inFile, outFile, [exe, "-c", inFile], shell=True, stdout=outFile)
        if(success and delinFile):
            remove_file(inFile)
    else:
        logging.info("Skipping bgzip. " + outFile + " already exists.");
    return success

