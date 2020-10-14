#
#
# Copyright (C) 2020 Niko Popitsch.  All rights reserved.
#
# This file is part of Nanopanel2
# 
# See the file LICENSE for redistribution information.
#
# @author: niko.popitsch

'''
    Main nanopanel2 file
'''


from argparse import ArgumentParser, RawDescriptionHelpFormatter
import binascii
from collections import OrderedDict
import datetime, time, sys, os, commentjson, pathlib
import itertools
import math
from multiprocessing import Pool, Process, Manager
import random
import shutil
import logging

import h5py
from matplotlib.backends.backend_pdf import PdfPages
import pysam
from scipy.stats import chi2_contingency, chisquare
from tqdm import tqdm
import vcfpy

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pyranges as pr
import seaborn as sns
from util.utils import files_exist, existing_file, pipeline_step, bgzip, sam2bam, sort_bam, remove_file

# bug in older pyranges version!
import pkg_resources
pkg_resources.require("pyranges>=0.0.77")

# Necessary for including python modules from a parent directory
sys.path.append(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))

LOGO = """
.=============================================.
                                          _ 
 ._ _  ___ ._ _  ___  ___  ___ ._ _  ___ | |
 | ' |<_> || ' |/ . \| . \<_> || ' |/ ._>| |
 |_|_|<___||_|_|\___/|  _/<___||_|_|\___.|_| 2
                     |_|                    
.============================ vVERS ==========.
"""

VERSION=open(os.path.dirname(os.path.abspath(__file__)) + '/VERSION', 'r').read()

# BAM flags, @see https://broadinstitute.github.io/picard/explain-flags.html
# @abstract the read is paired in sequencing, no matter whether it is mapped in a pair
BAM_FPAIRED=0x1
# ! @abstract the read is mapped in a proper pair 
BAM_FPROPER_PAIR=0x2
# ! @abstract the read itself is unmapped; conflictive with BAM_FPROPER_PAIR 
BAM_FUNMAP=0x4
# ! @abstract the mate is unmapped 
BAM_FMUNMAP=0x8
# ! @abstract the read is mapped to the reverse strand 
BAM_FREVERSE=0x10
# ! @abstract the mate is mapped to the reverse strand 
BAM_FMREVERSE=0x20
# ! @abstract this is read1 
BAM_FREAD1=0x40
# ! @abstract this is read2 
BAM_FREAD2=0x80
# ! @abstract not primary alignment 
BAM_FSECONDARY=0x100
# ! @abstract QC failure 
BAM_FQCFAIL=0x200
# ! @abstract optical or PCR duplicate 
BAM_FDUP=0x400
# ! @abstract optical or PCR duplicate 
BAM_SUPPLEMENTARY=0x800


# number of ref bases considered before current position for HP calculation
FP_CTX = 5

# maximum deletion length
MAX_DEL = 100

# ref bases
BASES = ['A', 'C', 'G', 'T']

# variant types
TYPES = ['DEL', 'INS', 'SNV']

# whether to use Yate's correction for SB calculation
YATES_CORRECTION=False

# default TSV headers
TSV_HEADERS = [
                "chr",
                "pos",
                "ref",
                "alt",
                "type", 
                "is_filtered",
                "AF",
                "DP",
                "BQ",
                "SB",
                "AQ1",
                "AQ2",
                "HP",
                "LQ",
                "SI",
                "CO",
                "quality",
                "vaf_raw",
                "aa_count_raw",
                "vaf_corr",
                "aa_count_corr",
                "aa_skew",
                "sb_pvalue",
                "hp_frac",
                "hp_len",
                "avg_base_qual",
                "ctx",
                "is_called"                 
                ]

# ###############################################################################
#                Default config values
# ###############################################################################
MANDATORY_CONF = [
        "dataset_name", "ref", "mappers"
    ]
DEF_CONF = {
        "threads":              1,
        "max_dp":               100000, # max read depth
        "min_dp":               10, # minimum read depth
        "min_base_quality":     0, # minimum base quality
        "suppress_snv":         [], # SNV calls with this filter will be suppressed in the VCF output
        "suppress_del":         [], # SNV calls with this filter will be suppressed in the VCF output
        "suppress_ins":         [], # SNV calls with this filter will be suppressed in the VCF output
        "basecall_grp":         "Basecall_1D_001", # the used basecalling grp
        "tag_pA":               "pA", # BAM TAG used to store pA values
        "tag_pC":               "pC", # BAM TAG used to store pC values
        "tag_pG":               "pG", # BAM TAG used to store pG values
        "tag_pT":               "pT", # BAM TAG used to store pT values
        "max_h5_cache":         10, # maximum number of cached H5 files
        
        "min_alt_dp_snv":       10, # min number of AA reads for SNVs
        "min_vaf_snv":          0.01, # min VAF for SNVs
        "max_aaskew_snv":       1.0, # max aa_skew for SNVs
        "min_plus_read_frac":   0.2, # for read strand imbalance filtering
        "max_plus_read_frac":   0.8, # for read strand imbalance filtering
        "min_corr_frac":        0.3, # for filtering calls with too many corrected reads

        "min_alt_dp_del":       100, # min number of AA reads for DEL
        "min_vaf_del":          0.05, # min VAF for DEL
        "max_aaskew_del":       0.2, # max aa skew for DEL
        
        "min_alt_dp_ins":       100, # min number of AA reads for INS
        "min_vaf_ins":          0.03, # min VAF for INS
        "max_aaskew_ins":       0.5, # max aa skew for INS
        
        "max_hp_run_snv":       3, # for filtering in homopolymer runs
        "max_hp_run_del":       3, # for filtering in homopolymer runs
        "max_hp_run_ins":       3, # for filtering in homopolymer runs

        "min_avg_bq_snv":       12, # for filtering low-quality alleles
        "min_avg_bq_del":       12, # for filtering low-quality alleles
        "min_avg_bq_ins":       12, # for filtering low-quality alleles
        
        "min_quality_del":      3, # for filtering low overall call quality
        "min_quality_ins":      3, # for filtering low overall call quality
        "min_quality_snv":      10, # for filtering low overall call quality
        
        "min_read_len":         100, # minimum read length
        "max_read_len":         'auto', # maximum read length (longer reads will be split). This will be estimated from the max amplicon length if set to 'auto'

        "report_all_snv":       False, # if true, all calls (for all AAs) per position will be reported
        
        "add_gt_field":         True, # if true, 0/1 genotype (GT) entries will be written to the output VCF for compatibility with downstream tools (NOTE: that this does not reflect the actual GT of somatic variants)
        "debug_mode":           True, # if set, additional files useful for debugging/QC are produced
        
        "random_seed":          12345,
        "max_reads":           'unlimited'  # if set to a number the alignment will be (randomly) subsampled to this maximum number of reads
    }

# ###############################################################################
#                util
# ###############################################################################
class MyArgumentParser(ArgumentParser):
    """ for parsing commandline arguments """
    def error(self, message):
        print(LOGO.replace("VERS", VERSION))
        sys.stderr.write('error: %s\n' % message)
        sys.stderr.write("usage: %s\n" % self.usage)
        sys.exit(1)
        
def validate_config(config, mandatory=[]):
    """ Validates configuration file and sets default config keys """
    mandatory += MANDATORY_CONF
        
    for k in mandatory:
        if k not in config:    
            logging.error("Missing mandatory config param %s" % (k))
            sys.exit(1)
    for k in DEF_CONF.keys():
        if k not in config:            
            config[k] = DEF_CONF[k]
    return config
def get_exe(name, config):
    """ get executable path from config if any """
    return config["exe"][name] if 'exe' in config and name in config["exe"] else name
def rec_to_str(rec):
    """ String representation of a alignment record """
    return rec.CHROM + ":" + str(rec.POS)+rec.REF + ">" + ",".join([alt.value if alt.value else "." for alt in rec.ALT])
def tsv_to_str(rec):
    """ String representation of a TSV record """
    return rec[0] + ":" + str(rec[1])+rec[2] + ">" + rec[3]
def is_pass_rec(rec):
    """ Check whether an alignment record is filtered """
    return ( rec.FILTER==[] or rec.FILTER==["."] or rec.FILTER==["PASS"] )
def is_pass_tsv(tsv):
    """ Check whether a TSV record is filtered """
    return ( not tsv[5] ) # 5th column contains 'is_filtered'
def in_roi( rec, roi):
    """ check whether an alignment record overlaps with the given region of interest encoded in a pandas dataframe """
    return len( roi[(roi['Chromosome']==rec.CHROM ) & (roi['Start']<rec.POS) & (roi['End']>=rec.POS) ] ) > 0
def get_var_type( rec ):
    """ estimate variant type [SNV, DEL, INS] for alignment record """
    ref = rec.REF
    alt = rec.ALT[0].value
    return "SNV" if (len(ref)==1 and len(alt)==1) else "DEL" if len(ref)>len(alt) else "INS"
def no_none(a):
    """ return 'NA' if None was passed or the passed value otherwise """
    if a is None:
        return "NA"
    return a
def save_log(a, def_val):
    """ return log10(val) if val != 0 or the passed default value otherwise """
    if a == 0:
        return def_val
    return math.log10(a)
def get_str(s):
    """ Returns a string representation of the passed variable. If None was passed, 'NA' will be returned. If a list was passed, a string representation of the first element will be returned """ 
    if s is None:
        return "NA"
    if (type(s) is list) :
        return str(s[0])
    return str(s)
    
def filter_max_ascore_alignments(in_bam, out_bam):
    """ Filter alignments with maximum alignment score: this will ensure that there is only one primary alignment per readname in the
    output BAM to work-around last bug. This method will chose the read with maximum alignment score (in AS tag) """
    rn2as = {}
    # get max AS score per read
    with pysam.AlignmentFile(in_bam, "rb") as bam:  # @UndefinedVariable
        for read in bam:
            ascore = read.get_tag('AS') if read.has_tag('AS') else 0
            rn2as[read.query_name] = max(rn2as[read.query_name], ascore) if read.query_name in rn2as else ascore        
    # write BAM
    with pysam.AlignmentFile(in_bam, "rb") as bam:  # @UndefinedVariable
        with pysam.AlignmentFile(out_bam, "wb", template=bam) as out: # @UndefinedVariable
            for read in bam:
                ascore = read.get_tag('AS') if read.has_tag('AS') else 0
                if read.query_name not in rn2as:
                    continue
                if ascore != rn2as[read.query_name]:
                    continue
                out.write(read)

def get_ref_seq(ref_fasta, chrom, pos1, left=FP_CTX, right=MAX_DEL):
    """ Extract padded ref seq """
    start = pos1-left-1
    end = pos1+right
    prefix=""
    if start < 0:
        prefix = "N"*(-start)
        start = 0
    return prefix + ref_fasta.fetch(chrom, start, end).upper()

def show_prog(q, max_value):
    """ Show progress """
    prog = tqdm(total=max_value, unit=' pos', desc='Calling variants')
    while 1:
        try:
            to_add = q.get(timeout=.1)
            q.task_done()
            if prog.n + to_add >= max_value:
                break
            prog.update(to_add)
        except:
            continue

def count_hp(seq):
    """ Counts HP length (any allele) from start """
    hp_char=None
    hp_count=0
    for base in seq:
        if hp_char == None:
            hp_char = base
        if base == hp_char:
            hp_count+=1
        else:
            break
    return hp_count, hp_char

def count_hp_alt(seq, allele):
    """ Counts the HP len of passed allele from start """
    hp_len=0
    while len(seq)>len(allele):
        if not seq.startswith(allele):
            break
        hp_len+=1
        seq = seq[len(allele):]
    return hp_len
    
def get_max_outliers(data, overlapping_data=[], min_n=5, min_count=0):
    """ find max outliers (i.e., potentially real insertions/deletions) or any allele with a given min count """
    test_data = data + overlapping_data
    if (len(test_data)<min_n):
        # not enough data to test for outliers: report all alleles with a min count
        out_max = [ (aa, c) for aa, c in data if c > min_count  ] 
    else:
        # report only outlier alleles with a min countt
        q75 = max( np.percentile(test_data, 75 ), min_count )
        out_max = [ (aa, c) for aa, c in data if c > q75  ] 
    return out_max


# ###############################################################################
#                Build readname / barcode index
# ###############################################################################  
def build_demux_idx(in_dir, out_file):
    """ Build readname / barcode index """
    stats={}
    with open(out_file, 'w') as out:
        for file in os.listdir(in_dir):
            if file.endswith(".fastq") and file.startswith("BC"):
                barcode = file[:-len(".fastq")]
                print("barcode\t%s" % (barcode), file=out)
                with open(in_dir + "/" + file, 'r') as fp:
                    i = 0
                    for line in fp:
                        if i % 4 == 0:
                            read_name = 'read_' + line.split(' ')[0][1:] # split by space, take 1st element, remove '@' prefix
                            print("%s" % (read_name), file=out)
                            stats[barcode] = stats[barcode]+1 if barcode in stats else 1
                        i+=1
    logging.info("Barcode stats:")
    logging.info(stats)
    logging.info("Finished building demux index")       
# 
def load_demux_idx(in_file): 
    """ Load the demux index file """
    demux_index={}
    bc = None
    with open(in_file, 'r') as fp:
        for line in fp:
            if line.startswith('barcode'):
                bc = line.split("\t")[1].rstrip()
                demux_index[bc]=set()
            else:
                demux_index[bc].add(line.rstrip()) 
    return demux_index

def demultiplex_reads(config, outdir):
    """ Demultiplexes reads with porechop """
    logging.info("-----------------------------------------------")
    logging.info("Demultiplexing reads")
    logging.info("-----------------------------------------------")
    config = validate_config(config, ['fast5_dir', 'fastq_dir'])

    demux_idx_file = outdir + config["dataset_name"]+".demux.idx"
    if files_exist([demux_idx_file]):
        logging.warning("Found existing demux index file %s. Will not re-create" % demux_idx_file)
        demux_idx = load_demux_idx(demux_idx_file)
    else:
        # mk tmp dir
        tmp_dir =  outdir + "/tmp/"       
        if not os.path.exists(tmp_dir):
            logging.info("Creating dir " + tmp_dir)
            os.makedirs(tmp_dir)
        cmd = [get_exe("porechop", config), "-i", config['fastq_dir'], "-t", str(config["threads"]), "-b", tmp_dir]
        pipeline_step(config['fastq_dir'], None, cmd, shell=True) # NOTE: here we could check whether there was any demux output!
        # create index
        build_demux_idx(tmp_dir, demux_idx_file )
        demux_idx = load_demux_idx(demux_idx_file)
        # delete temp dir
        try:
            shutil.rmtree(tmp_dir)
        except OSError as e:
            print ("Could not remove %s - %s." % (e.filename, e.strerror))
    return demux_idx


# ###############################################################################
#                Dataset samples
# ###############################################################################
class Sample():
    """ Dataset samples """
    def __init__(self, config, outdir, name, barcode, read_names ):
        self.config = config
        self.outdir = outdir
        self.name = name
        self.barcode = barcode
        self.mappers = config['mappers'].keys()
        self.read_names = read_names        
        self.out_streams = {}

        # create filenames
        self.fastq_file = outdir + self.config["dataset_name"]+( "."+name if barcode else "" ) +".fq"
        self.read_idx   = outdir + self.config["dataset_name"]+( "."+name if barcode else "" ) +".idx"
        self.sam_file   = OrderedDict()
        self.bam_file   = OrderedDict()
        self.decorated_bam_file = OrderedDict()
        self.downsampled_bam_file = OrderedDict()
        self.eval_file  = OrderedDict()
        self.stats_file = OrderedDict()
        self.vcf_file   = OrderedDict()
        self.tsv_file   = OrderedDict()
        self.hap_file_prefix = OrderedDict()
        self.hapmat_file= OrderedDict()
        for m in self.mappers:
            self.sam_file[m]  = outdir + self.config["dataset_name"]+( "."+name if barcode else "" ) + "." + m + ".sam"
            self.bam_file[m]  = outdir + self.config["dataset_name"]+( "."+name if barcode else "" ) + "." + m + ".bam"
            self.decorated_bam_file[m] = outdir + self.config["dataset_name"]+( "."+name if barcode else "" ) + "." + m + ".decorated.bam"
            self.downsampled_bam_file[m] = outdir + self.config["dataset_name"]+( "."+name if barcode else "" ) + "." + m + ".decorated.downsampled.bam"
            self.eval_file[m]  = outdir + self.config["dataset_name"]+( "."+name if barcode else "" ) + "." + m + ".performance.tsv"
            self.stats_file[m] = outdir + self.config["dataset_name"]+( "."+name if barcode else "" ) + "." + m + ".stats.tsv"
            self.vcf_file[m]   = outdir + self.config["dataset_name"]+( "."+name if barcode else "" ) + "." + m + ".vcf"
            self.tsv_file[m]   = outdir + self.config["dataset_name"]+( "."+name if barcode else "" ) + "." + m + ".tsv"
            self.hap_file_prefix[m]   = outdir + '/haplo/' + self.config["dataset_name"]+( "."+name if barcode else "" ) + "." + m + ".hap"
        
    def open(self, token, out_file):
        self.out_streams[token] = open(out_file, "w") 
    def close(self, token):
        if self.out_streams[token]:
            self.out_streams[token].close()
            
    
# ###############################################################################
#                extract FASTQ 
# ###############################################################################
def extract_fastq(config, demux_index, outdir):
    """ Extracts FASTQ from FAST5 files """
    logging.info("-----------------------------------------------")
    logging.info("Extract FASTQ data")
    logging.info("-----------------------------------------------")
    config = validate_config(config, ['fast5_dir', ])
    
    # calc max read len
    max_read_len = config['max_read_len']
    if max_read_len == 'auto':
        ref_fasta = pysam.FastaFile( config['ref'] ) # @UndefinedVariable
        max_read_len = max(ref_fasta.lengths)
        logging.info("Estimated max_read_len from chromosome length to %i " % (max_read_len) )
    
    bcg = config['basecall_grp']
    
    # create dataset samples
    samples=[]
    if demux_index is not None:
        for barcode in demux_index.keys():
            if barcode in config['demultiplex']:
                s = Sample(config, outdir, config['demultiplex'][barcode], barcode, demux_index[barcode])            
                samples+=[s]
            else:
                logging.warn("Barcode %s found in demux_index but not in configuration. Sample will be skipped." % (barcode) )
    else:
        s = Sample(config, outdir, 'all', None, None)   
        samples+=[s]
    
    logging.info("Samples: %s" % ([s.name for s in samples]))
    
    if files_exist([a.fastq_file+'.gz' for a in samples] + [a.read_idx for a in samples]):
        logging.warning("Found existing fastq + index files %s and %s. Will not re-create" % ([a.fastq_file+'.gz' for a in samples], [a.read_idx for a in samples]) )
    else:
        filenames = set()
        logging.info("Reading FAST5 files from " + config['fast5_dir'])
        for file in os.listdir(config['fast5_dir']):
            if file.endswith(".fast5"):
                filenames.add(file)
        logging.info("Loaded " + str(len(filenames)) + " FAST5 files")

        h={}
        for s in samples:
            logging.info("Extracting fastq for %s" % s.name )
            s.open('fastq_file', s.fastq_file)
            s.open('read_idx', s.read_idx)
            
            print("read_name\tfast5_file", file=s.out_streams['read_idx'])
            for fn in tqdm(filenames, desc='Extracting FASTQ data for %s' % s.name):
                # read HFD5 file into memory
                fast5 = os.path.abspath(config['fast5_dir'] + "/" + fn )
                if not files_exist(fast5):
                    logging.error("ERROR: cannot find " + fast5)
                    sys.exit(1)
                #print("reading "+fast5)
                with h5py.File(fast5, 'r', driver="core",backing_store=False) as f:   
                    reads=set(f.keys())
                    if s.read_names is not None:
                        reads = set.intersection(reads, s.read_names)
                    for readname in reads:
                        fastq=f[readname]["Analyses"][bcg]["BaseCalled_template"]["Fastq"][()]
                        fastq=fastq.decode('ascii')
                        
                        # split
                        rn, bases, rn2, qualities = fastq.splitlines()
                        frag=1
                        while len(bases)>max_read_len:
                            print("%s\n%s\n%s\n%s" % (rn, bases[0:max_read_len], rn2, qualities[0:max_read_len]), file=s.out_streams['fastq_file'] )
                            bases=bases[max_read_len:]
                            qualities=qualities[max_read_len:]
                            frag+=1
                        if ( len(bases) > config['min_read_len']):
                            print("%s\n%s\n%s\n%s" % (rn, bases, rn2, qualities), file=s.out_streams['fastq_file'] )
                        h[frag]=h[frag]+1 if frag in h else 1 
                        print("%s\t%s" % (readname, fast5 ), file=s.out_streams['read_idx'])
            s.close('fastq_file')
            s.close('read_idx')
            bgzip(s.fastq_file, delinFile=True, exe = get_exe("bgzip", config))
        logging.info("Histogram of fragments:")
        logging.info(h)
    return samples

# ###############################################################################
#                map reads
# ###############################################################################
def map_reads(config, samples):
    """ Maps reads using the configured long-read mappers """
    logging.info("-----------------------------------------------")
    logging.info("Map reads")
    logging.info("-----------------------------------------------")
    config = validate_config(config)
    
    ref = config["ref"]
    for m in config['mappers'].keys():
        for s in samples:
            if files_exist([s.bam_file[m]]):
                logging.warning("Found existing BAM file %s. Will not re-create" % s.bam_file[m])
            else:
                logging.info("Mapping %s with %s" % ( s.name, m) )
                additional_param = config['mappers'][m]['additional_param'] if 'additional_param' in config['mappers'][m] else []
                # MINIMAP2
                if config['mappers'][m]['type'] == 'minimap2':
                    cmd = [get_exe("minimap2", config), "-ax", "map-ont", "-t", str(config["threads"]), ref, s.fastq_file+'.gz', "-o", s.sam_file[m]] + additional_param
                    logging.info("%s" % " ".join(cmd) )
                    success = pipeline_step(s.fastq_file+'.gz', s.sam_file[m], cmd, shell=True)
                    sam2bam(s.sam_file[m], s.bam_file[m], index=True, sort=True, delinFile=True, override=False, exe=get_exe("samtools", config))
                # NGMLR
                elif config['mappers'][m]['type'] == 'ngmlr':
                    cmd = [get_exe("ngmlr", config), "-r", ref, "-q", s.fastq_file+'.gz', "-x", "ont", "-t", str(config["threads"]), "-o", s.sam_file[m]] + additional_param
                    logging.info("%s" % " ".join(cmd) )
                    success = pipeline_step(s.fastq_file+'.gz', s.sam_file[m], cmd, shell=True)
                    sam2bam(s.sam_file[m], s.bam_file[m], index=True, sort=True, delinFile=True, override=False, exe=get_exe("samtools", config))
                # LAST
                elif config['mappers'][m]['type'] == 'last':
                    # get last_db. Usually this is named as the ref file
                    last_db = config['mappers'][m]['last_db'] if 'last_db' in config['mappers'][m] else os.path.splitext(ref)[0]+'.last_db'
                    # step1 : align
                    cmd = [get_exe("lastal", config), "-Q1", last_db, s.fastq_file+'.gz']
                    tmp_file = s.sam_file[m]+'.tmp'
                    success = pipeline_step(s.fastq_file+'.gz', tmp_file, cmd, shell=True, stdout=tmp_file)
                    # step 2: split and create MAF file
                    cmd = [get_exe("last-split", config), "<", tmp_file]
                    maf_file = s.sam_file[m]+'.maf'
                    success = pipeline_step(tmp_file, tmp_file, cmd, shell=True, stdout=maf_file)
                    # step 3: create SAM
                    with open(s.sam_file[m], 'w') as out:
                        print("@HD\tVN:1.0\tSO:coordinate\t", file = out)
                        print("@PG\tID:last\tPN:last\tVN:i1042\tCL:lastal -Q1 %s %s | last-split" % (last_db, s.fastq_file+'.gz') , file = out)
                    cmd = [get_exe("maf-convert", config), "SAM",maf_file]
                    success = pipeline_step(maf_file, s.sam_file[m], cmd, shell=True, stdout=s.sam_file[m], append=True)
                    if success:
                        remove_file([tmp_file, maf_file])
                    # step 4 : convert to bam and add SQ dict from ref.fai file
                    cmd = [get_exe("samtools", config), "view", "-bt", ref+'.fai', s.sam_file[m] ]
                    bam_tmp = s.bam_file[m]+".tmp.bam"
                    success = pipeline_step(s.sam_file[m], bam_tmp, cmd, shell=True, stdout=bam_tmp)
                    if success:
                        remove_file([s.sam_file[m]])
                    # step5: select only the max-scoring alignments as maf-convert will include all of them (not marked as secondary alignments)
                    bam_tmp2 = s.bam_file[m]+".tmp2.bam"
                    filter_max_ascore_alignments( bam_tmp, bam_tmp2)
                    remove_file([bam_tmp])
                    # step 6: sort bam
                    sort_bam( bam_tmp2, s.bam_file[m], index=True, delinFile=True)
                else:
                    logging.error("ERROR: unknown/unsupported mapper configured: %s" % m)
                    sys.exit(1)
    return samples 

def get_read_count(bam_file):
    """ Quick estimate of read count in a BAM file """
    lines = pysam.idxstats(bam_file).splitlines()  # @UndefinedVariable
    total_records = sum([int(l.split("\t")[2]) for l in lines if not l.startswith("#")])
    return total_records

def prob2ascii(dat_max):
    """ Convert a probability value to an ASCII representation """
    pA = binascii.hexlify(bytes(dat_max[:,0].tolist()))
    pC = binascii.hexlify(bytes(dat_max[:,1].tolist()))
    pG = binascii.hexlify(bytes(dat_max[:,2].tolist()))
    pT = binascii.hexlify(bytes(dat_max[:,3].tolist()))
    return pA,pC,pG,pT

def ascii2prob(s):
    """ Converts an ASCII representation to a probability value """
    return binascii.unhexlify(s)

# Table of reverse complement bases
COMP_TABLE = {
    "A": 'T', "C": 'G', "T": 'A', "G": 'C'
    }
def reverse_complement(seq):
    """ Calculate reverse complement DNA sequence """
    rev=[]
    for c in seq[::-1]:
        c=c.upper()
        rev+=[COMP_TABLE.get(c, 'N')]
    return ''.join(rev)

# ###############################################################################
#                decorate reads
# ###############################################################################
def decorate_reads(config, samples):
    """ Adds tags containing guppy probabilities to a BAM file """ 
    logging.info("-----------------------------------------------")
    logging.info("Extract guppy probabilities")
    logging.info("-----------------------------------------------")
    config = validate_config(config)
  
    h5_cache = {}
    for m in config['mappers'].keys():
        for s in samples:
            # iterate bam
            if files_exist([s.decorated_bam_file[m]]):
                logging.warning("Found existing decorated BAM file %s. Will not re-create..." % s.decorated_bam_file[m])
            else:
                logging.info("Decorating %s / %s" % ( s.name, m ) )
                bcg = config['basecall_grp']
                # read_name 2 fast5_file index
                r2f = pd.read_csv(s.read_idx, sep='\t').set_index('read_name').to_dict()['fast5_file']
                # progressbar: get readcount from bam
                pbar = tqdm(total=get_read_count(s.bam_file[m]), desc='Decorating reads in %s/%s' % (s.name, m))
                max_h5_cache = config['max_h5_cache']
                with pysam.AlignmentFile(s.bam_file[m], "rb") as bam:  # @UndefinedVariable
                    with pysam.AlignmentFile(s.decorated_bam_file[m], "wb", template=bam) as out:  # @UndefinedVariable
                        for read in bam:
                            rid = 'read_' + read.query_name
                            if rid not in r2f:
                                logging.warning("WARN: read %s not found in index. Skipping...")
                            else:
                                # load + cache h5 
                                if r2f[rid] not in h5_cache:
                                    if len(h5_cache) == max_h5_cache:
                                        # remove oldest entry
                                        h5_cache[next(iter(h5_cache))].close()
                                        del h5_cache[next(iter(h5_cache))]
                                    #print("%s Loading %s" % (read.reference_name, h5_filename) )
                                    h5_cache[r2f[rid]]=h5py.File(r2f[rid], 'r', driver="core",backing_store=False)
                                    
                                h5_file = h5_cache[r2f[rid]]
                                if rid not in h5_file.keys():
                                    logging.warning("WARN: read %s not found in fast5 (maybe opened by another process?). Skipping...")
                                fastq=h5_file[rid]["Analyses"][bcg]["BaseCalled_template"]["Fastq"][()].splitlines()[1].decode('ascii')
                                if fastq is None or read.query_sequence is None: # weird special case w/o read query sequence or missing fastq data
                                    tqdm.write('Could not extract fastq for %s. Skipping read.' % (rid) )
                                    continue
                                trace = np.array(h5_file[rid]["Analyses"][bcg]["BaseCalled_template"]["Trace"])
                                move = np.array(h5_file[rid]["Analyses"][bcg]["BaseCalled_template"]["Move"])
                                dat = trace[move == 1] # mask data
                                # get max flip/flop (A, C, G, T, A', C', G' and T')
                                if read.is_reverse:
                                    # reverse complement!
                                    rev = np.array(list(reversed(dat)))
                                    dat_max = np.column_stack(( np.maximum(rev[:,3], rev[:,7]), 
                                                                np.maximum(rev[:,2], rev[:,6]),
                                                                np.maximum(rev[:,1], rev[:,5]),
                                                                np.maximum(rev[:,0], rev[:,4]) ))
                                    fastq = reverse_complement(fastq)
                                else:
                                    dat_max = np.column_stack(( np.maximum(dat[:,0], dat[:,4]), 
                                                                np.maximum(dat[:,1], dat[:,5]),
                                                                np.maximum(dat[:,2], dat[:,6]),
                                                                np.maximum(dat[:,3], dat[:,7]) ))
                                    
                                # now correct for hard-clipped bases that are not in BAM file.
                                try:
                                    bam_start=str(fastq).index(read.query_sequence)
                                except ValueError:
                                    logging.error("ERROR: cannot correct for hard-clipped bases in " + rid)
                                    sys.exit(1)
                                bam_end=bam_start+len(read.query_sequence)   
                                
                                pA,pC,pG,pT = prob2ascii( dat_max[bam_start:bam_end] )
        
                                read.set_tag(config['tag_pA'], pA)
                                read.set_tag(config['tag_pC'], pC)
                                read.set_tag(config['tag_pG'], pG)
                                read.set_tag(config['tag_pT'], pT)
                            # write read
                            out.write(read)
                            pbar.update(1)
                                         
                # index bam file
                pysam.index(s.decorated_bam_file[m])   # @UndefinedVariable
            
    # close h5 files
    for h5_file in h5_cache.values():
        h5_file.close() 
    return samples

# ###############################################################################
#                downsample reads
# ###############################################################################
def downsample_reads(config, samples):
    """ Downsamples a BAMN file """

    if ('max_reads' in config) and (config['max_reads'] != 'unlimited'):
        logging.info("-----------------------------------------------")
        logging.info("Downsampling reads")
        logging.info("-----------------------------------------------")
        config = validate_config(config)
      
        for m in config['mappers'].keys():
            for s in samples:
                if files_exist([s.downsampled_bam_file[m]]):
                    logging.warning("Found existing downsampled BAM file %s. Will not re-create..." % s.downsampled_bam_file[m])
                else:
                    # get readcount from bam
                    nreads = get_read_count(s.decorated_bam_file[m])
                    max_reads = config['max_reads']
                    frac = max_reads / nreads
                    if frac >=1:
                        # no need to downsample!
                        s.downsampled_bam_file[m] = s.decorated_bam_file[m]
                    else:
                        written = 0
                        pbar = tqdm(total=nreads, desc='Downsampling reads in %s/%s' % (s.name, m))
                        with pysam.AlignmentFile(s.decorated_bam_file[m], "rb") as bam:  # @UndefinedVariable
                            with pysam.AlignmentFile(s.downsampled_bam_file[m], "wb", template=bam) as out:  # @UndefinedVariable
                                for read in bam:
                                    if ( random.random() < frac ) and ( written < max_reads ):
                                        # write read
                                        out.write(read)
                                        written+=1
                                    pbar.update(1)
                    # index bam file
                    pysam.index(s.downsampled_bam_file[m])   # @UndefinedVariable
    else:
        logging.info("No downsampling configured")
        # no downsampling
        for m in config['mappers'].keys():
            for s in samples:
                s.downsampled_bam_file[m] = s.decorated_bam_file[m]
        
    return samples

# ###############################################################################
#                Call variants 
# ###############################################################################

# .........................................................................
# Variant call 
# .........................................................................
class VariantCall():
    """ A variant call """
    
    # SNVs only
    # calc contingency matrix and strand bias 
    def calc_sb(self, alt_base, plus_read_frac):
        ref_base = self.s.ref_bases[FP_CTX]
        
        c_alt_all=self.s.counts[alt_base] if alt_base in self.s.counts else 0
        c_alt_plus=self.s.counts_plus[alt_base] if alt_base in self.s.counts_plus else 0
        
        c_ref_all=self.s.counts_ref[ref_base] if ref_base in self.s.counts_ref else 0
        c_ref_plus=self.s.counts_ref_plus[ref_base] if ref_base in self.s.counts_ref_plus else 0
        
        # add pseudo-counts. array of [+.-] reads
        c_ref=[ c_ref_plus + 0.0000000001, c_ref_all - c_ref_plus + 0.0000000001 ] 
        c_alt=[ c_alt_plus + 0.0000000001, c_alt_all - c_alt_plus + 0.0000000001 ] 

        # correct by plus_read frac
        c_ref_corr = [ c_ref[0] * ( 1 - plus_read_frac ) * 2, c_ref[1] * ( plus_read_frac ) * 2 ]
        c_alt_corr = [ c_alt[0] * ( 1 - plus_read_frac ) * 2, c_alt[1] * ( plus_read_frac ) * 2 ]
                
        # calc aa skew from corrected counts for expected 50:50 ration of read strands
        aa_skew = math.log2( c_alt_corr[0] / c_alt_corr[1] ) if c_alt_corr[0] > 0 and c_alt_corr[1] > 0 else None
        
        # calc SB
        if (c_alt_corr[0] > 0 and c_alt_corr[1] > 0):
            _, sb_pvalue, _, _ = chi2_contingency([c_ref_corr, c_alt_corr], correction=YATES_CORRECTION) 
        else:
            _, sb_pvalue, _, _ =  None, 0, None, None

        return aa_skew, sb_pvalue, c_ref, c_alt, c_ref_corr, c_alt_corr

    def calc_corrected_aa_count(self, alt_base):
        """ SNVs only. Calc corrected read counts and pX for non-AA and chi2 test """
        if not alt_base:
            return None, None
        AB = list(BASES)
        AB.remove(alt_base)
        corrected_aa_count = 0
        non_alt_avg_prob=[]
        for i, called_base in enumerate(self.s.meta['called_base']):
            if called_base == alt_base:
                prob=[]
                for b in AB:
                    prob += [self.s.meta['p'+b][i]]
                prob_alt = self.s.meta['p'+alt_base][i]
                #logging.info("prob %s / %i = %i" % (prob, prob_alt, sum(prob)+prob_alt))
                non_alt_avg_prob+=[sum(prob)/3/255]
                if sum(prob)>0 and sum(prob)+prob_alt > 200: # ignore low-prob reads
                    _,p = chisquare(prob)
                    if p > 0.05: # the probs of all other possible bases are equally distributed: consider as high quality alt_base
                        corrected_aa_count+=1
                else:
                    corrected_aa_count+=1
        return corrected_aa_count, corrected_aa_count/self.s.dp, np.mean(non_alt_avg_prob) if len(non_alt_avg_prob)>0 else 0
                   

    def calc_aa_skew(self, alt_base, plus_read_frac):
        """ INDELs only. Calculate AA skew """
        c_alt_all=self.s.counts[alt_base] if alt_base in self.s.counts else 0
        c_alt_plus=self.s.counts_plus[alt_base] if alt_base in self.s.counts_plus else 0
        # add pseudo-counts. array of [+.-] reads
        c_alt=[ c_alt_plus + 0.0000000001, c_alt_all - c_alt_plus + 0.0000000001 ] 

        # correct counts for expected 50:50 ration of read strands
        c_alt_corr = [
                c_alt[0] * (1-plus_read_frac) * 2,
                c_alt[1] * (plus_read_frac) * 2
            ]
        aa_skew = math.log2(c_alt_corr[0]/c_alt_corr[1]) if c_alt_corr[0] > 0 and c_alt_corr[1] > 0 else None
        #logging.info("%s:%i +: %i -: %i aas: %f" % ( self.chr, self.ref_pos1, plus, minus, aa_skew ) )        
        return aa_skew, c_alt, c_alt_corr
    
    def calc_hp_filter(self, max_hp_run):
        """ calculate HP filter """
        seq_eq_ins = False
        ctx=self.s.ref_bases
        if self.type =='SNV':
            # count left HP run
            hp_count_left, hp_char_left = count_hp(ctx[:FP_CTX][::-1])
            hp_count_right, hp_char_right = count_hp(ctx[FP_CTX+1:])
            filtered = ( hp_count_left>=max_hp_run or hp_count_right >=max_hp_run-1 ) and ( (hp_char_left ==  self.alt) or (hp_char_right ==  self.alt) )
            hp_len = max(hp_count_left,  hp_count_right)
        elif self.type=='INS':
            allele =  self.alt[len( self.ref):]
            if len(set(allele)) == 1: # for hp-alleles: search only for first char
                allele=allele[0]
            hp_count_left = count_hp_alt(ctx[:FP_CTX][::-1], allele[::-1])
            hp_count_right = count_hp_alt(ctx[FP_CTX+1:], allele)
            hp_len=hp_count_left + hp_count_right
            filtered = ( hp_len >=max_hp_run )    
            # inserted allele is exactly same as the one following the current position
            if self.alt == ctx[FP_CTX:FP_CTX+len(self.alt)]:
                seq_eq_ins=True   
        elif self.type=='DEL':
            allele =   self.ref[len( self.alt):]
            if len(set(allele)) == 1: # for hp-alleles: search only for first char
                allele=allele[0]
            hp_count_left = count_hp_alt(ctx[:FP_CTX][::-1], allele[::-1])
            hp_count_right = count_hp_alt(ctx[FP_CTX+1:], allele)
            hp_len=hp_count_left + hp_count_right
            #print(allele, ctx[:FP_CTX][::-1], hp_count_left, ctx[FP_CTX+1:], hp_count_right)
            filtered = ( hp_len >=max_hp_run )      
        else:
            filtered=False
            hp_len=0
        return filtered, hp_len, seq_eq_ins

    
    def __init__(self, pileup_stats, aa, vtype):
        self.s = pileup_stats
        self.config = self.s.config
        self.type = vtype
        self.filters=[]
        self.aa_skew, self.aa_count_corr, self.vaf_corr, self.hp_frac, self.hp_len, self.sb_pvalue = [0]*6
        self.c_ref, self.c_alt, self.c_ref_corr, self.c_alt_corr, self.non_alt_avg_prob = [None] * 5
        self.seq_eq_ins = False # true if ref sequence at this pos is equal to the inserted bases
        
        # get AA counts
        self.alt_allele, self.aa_count_raw = aa
    
        # calc raw VAF
        self.vaf_raw = self.aa_count_raw / self.s.dp if self.s.dp > 0 else 0

        # calculate expected fraction of + strand reads (plus_read_frac) 
        self.plus_read_frac = self.s.count_plus_reads / self.s.dp 
        
        # calculate average per-base quality (SNV and INS only)
        self.avg_base_qual = None
        if self.type == 'SNV' and self.alt_allele is not None:
            aa_qual=[]
            for i, called_base in enumerate(self.s.meta['called_base']):
                if called_base == self.alt_allele:
                    aa_qual+=[self.s.meta['qual'][i]]
            self.avg_base_qual = np.mean(aa_qual)
        elif self.type == 'INS' and self.alt_allele is not None:
            self.avg_base_qual = np.mean(self.s.ins_quals[self.alt_allele])

        # filter variant
        if self.alt_allele:
            
            # calc ref/alt allele
            self.ref, self.alt = self.get_ref_alt()
 
            # thresholds
            min_alt_dp=self.config ['min_alt_dp_del'] if self.type=='DEL' else self.config ['min_alt_dp_ins'] if self.type=='INS' else self.config ['min_alt_dp_snv']
            min_vaf=self.config ['min_vaf_del'] if self.type=='DEL' else self.config ['min_vaf_ins'] if self.type=='INS' else self.config ['min_vaf_snv']
            max_aaskew=self.config ['max_aaskew_del'] if self.type=='DEL' else self.config['max_aaskew_ins'] if self.type=='INS' else self.config ['max_aaskew_snv']
            max_hp_run=self.config['max_hp_run_del'] if self.type=='DEL' else self.config['max_hp_run_ins'] if self.type=='INS' else self.config['max_hp_run_snv']
            min_avg_bq=self.config['min_avg_bq_del'] if self.type=='DEL' else self.config['min_avg_bq_ins'] if self.type=='INS' else self.config['min_avg_bq_snv']
            min_quality=self.config['min_quality_del'] if self.type=='DEL' else self.config['min_quality_ins'] if self.type=='INS' else self.config['min_quality_snv']
           
            # filter if not enough reads or not enough vaf
            if self.vaf_raw < min_vaf:
                self.filters+=["AF"] # low AF
            
            # minimum depth filter
            if self.aa_count_raw < min_alt_dp:
                self.filters+=["DP"] # low DP
                
            # average base quality filter
            if self.avg_base_qual is not None and self.avg_base_qual < min_avg_bq:
                self.filters+=["BQ"] # low BQ
                
            # strand imbalance filter
            if self.plus_read_frac < config['min_plus_read_frac'] or self.plus_read_frac > config['max_plus_read_frac']:
                self.filters+=["SI"] # high imbalance of read strands

            # sequence context (+/- FP_CTX bases)
            self.CTX = self.s.ref_bases[:FP_CTX*2+1]

            # calculate other filters only if not filtered yet
            if (not self.is_filtered()) or (config['debug_mode'] is True):
                
                if self.type == 'SNV':
                    # calc strand bias   
                    self.aa_skew, self.sb_pvalue, self.c_ref, self.c_alt, self.c_ref_corr, self.c_alt_corr = self.calc_sb( self.alt_allele, self.plus_read_frac ) 
                            
                    # calculate corrected counts dfrom reads with expected pX distribution for non-ref reads
                    self.aa_count_corr, self.vaf_corr, self.non_alt_avg_prob = self.calc_corrected_aa_count(self.alt_allele)
                else: # INDELs
                    # aa skew   
                    self.aa_skew, self.c_alt, self.c_alt_corr = self.calc_aa_skew( self.alt_allele, self.plus_read_frac) 

                
                # filter if aa Skew is too large (note sb_pvalue will be 0 for INDELs)
                if ( self.aa_skew is not None) and ( abs(self.aa_skew) > max_aaskew ) and ( self.sb_pvalue <= 0.05 ):
                    self.filters += ["SB"]
                #logging.info("%s:%i %f %f %s" % ( self.s.chr, self.s.ref_pos1, self.aa_skew, self.sb_pvalue, self.filters) )
    
                # HP length
                hp_filter, self.hp_len, self.seq_eq_ins = self.calc_hp_filter( max_hp_run )
                        
                # homopolymer filter
                if hp_filter:
                    self.filters += ["HP"] 
                    
                if self.type == 'SNV':
                    # calc corrected counts and check AF again   
                    if ( self.vaf_corr is not None and self.vaf_corr < config['min_vaf_snv'] ) or ( self.aa_count_corr is not None and self.aa_count_corr < config['min_alt_dp_snv']):
                            self.filters += ["AQ1"] # alt-allele quality filter 1
                
                    # filter SNV calls with too many corrected calls
                    if self.aa_count_corr is not None and self.aa_count_raw > 0 and self.aa_count_corr / self.aa_count_raw < config['min_corr_frac']:
                        self.filters += ["AQ2"] # alt-allele quality filter 2

                # calc variant quality
                self.quality = self.calc_quality()
                if self.quality<min_quality:
                    self.filters += ["LQ"] 
            else:
                self.quality = None 
 
    def calc_quality(self):
        """ Calculates Phred-score call quality values """ 
        quality = None
        if (self.type == 'SNV') and (self.non_alt_avg_prob is not None):
            quality = min(100, -10 * math.log10( self.non_alt_avg_prob ) if (self.non_alt_avg_prob>0) else 0 )
        elif (self.type == 'INS'):
            score=[]
            score += [1 - min(self.hp_len, 5)/5] # long HP = low score
            score += [0] if self.seq_eq_ins else [1]  # ref seq equals inserted bases = low score
            score += [1-min(5, abs(self.aa_skew))/5] # high aaSkew = low score
            if self.avg_base_qual is not None:
                score += [min(30, self.avg_base_qual)/30] # high avg base quality = high score
            quality = min(100, -10 * math.log10( 1 - np.mean(score) ) )
        elif (self.type == 'DEL'):
            score=[]
            score += [1 - min(self.hp_len, 5)/5] # long HP = low score
            score += [1-min(5, abs(self.aa_skew))/5] # high aaSkew = low score
            if self.avg_base_qual is not None:
                score += [min(30, self.avg_base_qual)/30] # high avg base quality = high score
            quality = min(100, -10 * math.log10( 1 - np.mean(score) ) )
        return quality
 
    def is_filtered(self):
        """ True if call is filtered """
        return len(self.filters)>0
    
    def is_nocall(self):
        """ True if this is a no-call """
        return self.alt_allele is None
    
    def get_ref_alt(self):
        """ get reference and alternate allele strings """
        if self.type == 'DEL':
            ref = self.s.ref_bases[FP_CTX:FP_CTX+(-self.alt_allele)+1] # note: alt-allele is - number of deleted bases (e.g., -3)
            alt = self.s.ref_bases[FP_CTX]
        elif self.type == 'INS':
            ref = self.s.ref_bases[FP_CTX]
            alt = self.s.ins_bases[self.alt_allele] # note: alt-allele is number of inserted bases (e.g., +5)
        elif self.type == 'SNV':
            ref = self.s.ref_bases[FP_CTX]
            alt = self.alt_allele
        else:
            sys.exit("unknown variant type %s" % self.type)
        return ref, alt
    

    def to_vcf(self):
        """ Convert to a vcfpy record """
        INFO=vcfpy.OrderedDict()
        INFO['SOMATIC']=1
        INFO['AF'] =   (str(self.vaf_raw))
        INFO['CNT_RAW'] =   (str(self.aa_count_raw))
        INFO['VAF_CORR'] =  (str(self.vaf_corr))
        INFO['CNT_CORR'] =  (str(self.aa_count_corr))
        INFO['AA_SKEW'] =   (str(self.aa_skew))
        INFO['SB_P'] =      (str(self.sb_pvalue))
        INFO['HP_FRAC'] =   (str(self.hp_frac))        
        INFO['HP_LEN'] =    (str(self.hp_len))
        INFO['TYPE'] =      (self.type)
        INFO['BASE_QUAL'] = (str(self.avg_base_qual))
        INFO['CTX'] =       (self.CTX)


        # for debungging only
        if type == 'SNV':
            if self.c_alt:
                INFO['DEBUG'] = '|'.join( [str(y)+'_'+str(x[0])+'/'+str(x[1]) for x,y in zip([self.c_ref, self.c_alt, self.c_ref_corr, self.c_alt_corr], ['c_ref','c_alt','c_ref_corr','c_alt_corr'])] )
        else:
            if self.c_alt:
                INFO['DEBUG'] = '|'.join( [str(y)+'_'+str(x[0])+'/'+str(x[1]) for x,y in zip([self.c_alt, self.c_alt_corr], ['c_alt','c_alt_corr'])] )            
        
        # calls 
        vcf_calls = [vcfpy.Call(config['dataset_name']+'.'+self.s.sample_name, vcfpy.OrderedDict([('GT', '1'), ('DP', str(self.s.dp))]))] if config['add_gt_field'] else None

        #build vcf record
        vcf_record = vcfpy.Record(CHROM=self.s.chr, 
                                  POS=self.s.ref_pos1, 
                                  ID=[], 
                                  REF=self.ref, 
                                  ALT=[vcfpy.Substitution(type_=self.type, value=self.alt)] , 
                                  QUAL=self.quality, 
                                  FILTER=self.filters, 
                                  INFO=INFO, 
                                  FORMAT=["GT", "DP"] if config['add_gt_field'] else [],
                                  calls=vcf_calls )
        
        return vcf_record
 
    def to_tsv(self):
        """ Convert to a TSV record """
        return ([
            self.s.chr,
            self.s.ref_pos1,
            self.ref,
            self.alt,
            self.type, 
            1 if self.filters else 0,
            1 if "AF" in self.filters else 0,
            1 if "DP" in self.filters else 0,
            1 if "BQ" in self.filters else 0,
            1 if "SB" in self.filters else 0,
            1 if "AQ1" in self.filters else 0,
            1 if "AQ2" in self.filters else 0,
            1 if "HP" in self.filters else 0,
            1 if "LQ" in self.filters else 0,
            1 if "SI" in self.filters else 0,
            1 if "CO" in self.filters else 0,
            no_none(self.quality),
            no_none(self.vaf_raw),
            no_none(self.aa_count_raw),
            no_none(self.vaf_corr),
            no_none(self.aa_count_corr),
            no_none(self.aa_skew),
            no_none(self.sb_pvalue),
            no_none(self.hp_frac),
            no_none(self.hp_len),
            no_none(self.avg_base_qual),
            self.CTX,
            1 ])  
        

def create_missing_record(truth_var_rec):
    """ Creates a FN TSV entry for missing calls """
    return ([
            truth_var_rec.CHROM,
            truth_var_rec.POS,
            truth_var_rec.REF,
            ",".join([alt.value if alt.value else "." for alt in truth_var_rec.ALT]),
            get_var_type(truth_var_rec), 
            'NA',
            'NA',
            'NA',
            'NA',
            'NA',
            'NA',
            'NA',
            'NA',
            'NA',
            'NA',
            'NA',
            0,
            'NA',
            'NA',
            'NA',
            'NA',
            'NA',
            'NA',
            'NA',
            'NA',
            'NA',
            'NA',
            0 ])
              
# ...................................................................
class PileupStats():
    """ Pileup statistics """
    def __init__(self, sample_name, chrom, ref_pos1, ref_bases, dp, config):
        self.sample_name = sample_name
        self.chr = chrom
        self.ref_pos1 = ref_pos1 
        self.ref_bases = ref_bases
        self.dp = dp
        self.counts = {} # counts of alt-allele (for SNVs) or # of reads with n deletions/insertions. E.g., 'A': 123, '-4', 12, '+3': 33
        self.counts_plus = {} # as above but only + reads
        self.counts_ref = {} # counts of ref-allele (for SNVs)
        self.counts_ref_plus = {} # as above but only + reads
        self.meta = {}
        for tag in ['called_base', 'qual', 'pA', 'pC', 'pG', 'pT', 'strand']:
            self.meta[tag]=[]
        self.ins_bases = {}
        self.ins_quals = {}
        self.config=config
        self.count_plus_reads = 0 # number of plus-strand reads

    # ...................................................................
    def call_variants(self, overlapping_del=[], overlapping_del_pos=None):
        """ Calculate additional measures, should be called after all data raw is collected """
        
        # sort alt-alleles by counts, e.g., [('A', 130), ('-1', 30), ('+2', 20), ('-12', 10)]
        mfaa = [(k,v) for k, v in sorted({i:self.counts[i] for i in self.counts}.items(), key=lambda item: item[1], reverse=True)]
        
        mfaa_snv = [(k,v) for k,v in mfaa if k in BASES]
        mfaa_del = [(k,v) for k,v in mfaa if k not in BASES and k < 0]
        mfaa_ins = [(k,v) for k,v in mfaa if k not in BASES and k > 0]

        calls = {}
        mfaa_calls = {}

        
        # insertions: iterate all alleles where the counts are max_outliers wrt. to all other noted alleles or 
        ins_alleles = get_max_outliers(mfaa_ins, overlapping_data=[], min_n=5, min_count=config['min_alt_dp_ins'])
        #logging.info("INS %s:%i %s | %s" % ( self.chr, self.ref_pos1, mfaa_ins, ins_alleles ) )
        for aa in ins_alleles:
            vtype = 'INS'
            call = VariantCall(self, aa, vtype)
            if vtype not in mfaa_calls:
                mfaa_calls[vtype] = call
         
        # deletions
        # update overlapping deletions
        if overlapping_del_pos:
            new_overlapping = []
            for i in range(0, self.ref_pos1 - overlapping_del_pos): # usually this is 1 
                for aa,c in overlapping_del:
                    if aa < -1:
                        new_overlapping += [(aa+1, c)]
            overlapping_del = new_overlapping
        #logging.info("DEL %s:%i %s | %s" % ( self.chr, self.ref_pos1, mfaa_del, overlapping_del ) )
        del_alleles = get_max_outliers(mfaa_del, overlapping_data=overlapping_del, min_n=5, min_count=config['min_alt_dp_del'])
        for aa in del_alleles:
            vtype = 'DEL'
            call = VariantCall(self, aa, vtype)
            if vtype not in mfaa_calls:
                mfaa_calls[vtype] = call          
        # update overlapping_deletions
        overlapping_del += mfaa_del
                
        # substitutions
        for aa in mfaa_snv:
            vtype = 'SNV'
            call = VariantCall(self, aa, vtype)
            if vtype not in mfaa_calls:
                mfaa_calls[vtype] = call
            if not call.is_filtered():
                calls[vtype]=calls[vtype] + [call] if vtype in calls else [call]

        # report all unfiltered calls and at least one (MFA) call per vtype
        for vtype in TYPES:
            if vtype not in calls and vtype in mfaa_calls:
                calls[vtype] = [mfaa_calls[vtype]]

        return calls, overlapping_del, self.ref_pos1

# ......................................................................... 
def call_variants_in_roi(args):
    """ Call variants in ROI """
    sample_name, chrom, start, end, bam_file, config, queue = args
    vcf_records = []
    tsv_records=[]
    uncalled_pos=[]
    debug_signal={}
    avg_dp=[] # for calculating average coverage
    read_names={}
    read_names['+']=set()
    read_names['-']=set()
    chr_stats={} # chr-wide stats
    overlapping_del=[] # for keeping track of previous deletions overlapping the current position
    overlapping_del_pos = None
    #logging.info("Extracting reads from '%s:%i-%i'" % (chr, start, end ) )

    bam = pysam.AlignmentFile(bam_file, "rb") # @UndefinedVariable
    ref_fasta = pysam.FastaFile( config['ref'] ) # @UndefinedVariable
        
    # @see https://pysam.readthedocs.io/en/latest/api.html
    flag_filter = BAM_FUNMAP | BAM_FSECONDARY | BAM_FQCFAIL | BAM_FDUP | BAM_SUPPLEMENTARY
    biter = bam.pileup(chrom, start, end, flag_filter=flag_filter, max_depth=config['max_dp'], truncate=True, mark_matches=False, mark_ends=False, add_indels=False, min_base_quality=config['min_base_quality'])
    for pileup_column in biter:
        chrom = pileup_column.reference_name
        ref_pos1=pileup_column.reference_pos+1
        ref_bases = get_ref_seq( ref_fasta, chrom, ref_pos1 ) # get (padded) context
        ref_base = ref_bases[FP_CTX] # this is the refc at the current position
        dp = pileup_column.get_num_aligned()
        
        # count fwd HP len
        if config['debug_mode']:
            hp_len1, _ = count_hp(ref_bases[:FP_CTX][::-1]) 
            hp_len2, _ =  count_hp(ref_bases[FP_CTX+1:])
            hp_len=hp_len1+hp_len2
            debug_signal['hp_len']=debug_signal['hp_len'] + [[chrom, ref_pos1-1, hp_len]] if 'hp_len' in debug_signal else [[chrom, ref_pos1-1, hp_len]]
                    
        if dp < config['min_dp']:
            # TODO: write uncalled to bed
            uncalled_pos+=[[chrom, ref_pos1-1, ref_pos1]]
        else:
            avg_dp+=[dp]
            stats = PileupStats( sample_name, chrom, ref_pos1, ref_bases, dp, config )

            avg_qual=[]      

            # iterate reads and collect stats
            for _, rec in enumerate(pileup_column.pileups):

                read_name = rec.alignment.query_name
                pos_in_read = rec.query_position
                read_is_del = rec.is_del or rec.is_refskip
                read_strand = "-" if rec.alignment.is_reverse else "+"   
                if read_strand=='+':
                    stats.count_plus_reads+=1
                # for counting +/- reads 
                read_names[read_strand].add(read_name)

                # Insertion
                if rec.indel>0:
                    stats.counts[rec.indel] = stats.counts[rec.indel]+1 if rec.indel in stats.counts else 1
                    if read_strand == '+':
                        stats.counts_plus[rec.indel] = stats.counts_plus[rec.indel]+1 if rec.indel in stats.counts_plus else 1
                    stats.ins_bases[rec.indel] = rec.alignment.query_sequence[pos_in_read:pos_in_read+rec.indel+1]
                    stats.ins_quals[rec.indel] = stats.ins_quals[rec.indel] + rec.alignment.query_qualities[pos_in_read:pos_in_read+rec.indel+1] if rec.indel in stats.ins_quals else rec.alignment.query_qualities[pos_in_read:pos_in_read+rec.indel+1]
                # Deletion
                elif rec.indel<0:
                    stats.counts[rec.indel] = stats.counts[rec.indel]+1 if rec.indel in stats.counts else 1
                    if read_strand == '+':
                        stats.counts_plus[rec.indel] = stats.counts_plus[rec.indel]+1 if rec.indel in stats.counts_plus else 1
                # SNVs
                else:
                    if not read_is_del:
                        read_base = rec.alignment.query_sequence[pos_in_read]
                        if read_base != ref_base: # count only alt-allele
                            stats.counts[read_base] = stats.counts[read_base]+1 if read_base in stats.counts else 1
                            if read_strand == '+':
                                stats.counts_plus[read_base] = stats.counts_plus[read_base]+1 if read_base in stats.counts_plus else 1
                        else:
                            stats.counts_ref[read_base] = stats.counts_ref[read_base]+1 if read_base in stats.counts_ref else 1
                            if read_strand == '+':
                                stats.counts_ref_plus[read_base] = stats.counts_ref_plus[read_base]+1 if read_base in stats.counts_ref_plus else 1
                        
                        base_quality = rec.alignment.query_qualities[pos_in_read]
                        pA = ascii2prob(rec.alignment.get_tag(config['tag_pA']))[pos_in_read]
                        pC = ascii2prob(rec.alignment.get_tag(config['tag_pC']))[pos_in_read]
                        pG = ascii2prob(rec.alignment.get_tag(config['tag_pG']))[pos_in_read]
                        pT = ascii2prob(rec.alignment.get_tag(config['tag_pT']))[pos_in_read]
                        
                        stats.meta['called_base'] += [read_base]
                        stats.meta["strand"] += [read_strand]
                        stats.meta["qual"] += [base_quality]
                        stats.meta["pA"] += [pA]
                        stats.meta["pC"] += [pC]
                        stats.meta["pG"] += [pG]
                        stats.meta["pT"] += [pT]
                        
                        #logging.info("%s/%s A %i + C %i + T %i + G %i = %i"%(ref_base, read_base, pA, pC, pT, pG, pA + pC + pT + pG))

                        avg_qual+=[base_quality]
                        
            avg_qual = np.mean(avg_qual)   

            if config['debug_mode']:
                debug_signal['avg_qual']=debug_signal['avg_qual'] + [[chrom, ref_pos1, avg_qual]] if 'avg_qual' in debug_signal else [[chrom, ref_pos1, avg_qual]]
            calls, overlapping_del, overlapping_del_pos = stats.call_variants(overlapping_del, overlapping_del_pos)

            # NOTE: write indel records before SNV record as it refers to previous position!
            for vtype in TYPES:
                if vtype in calls:
                    for call in calls[vtype]:
                        suppress = [ x.upper() for x in ( config['suppress_del'] if vtype == 'DEL' else config['suppress_ins'] if vtype == 'INS' else config['suppress_snv']) ]
                        # write to VCF?
                        if ( 'ALL' not in suppress ) and (len(set(call.filters).intersection( suppress ))==0):
                            vcf_records+=[call.to_vcf()]
                        # Write to TSV
                        tsv_records+=[call.to_tsv()]
                        # stats
                        if not call.is_filtered():
                            chr_stats[vtype] = chr_stats[vtype] + 1 if vtype in chr_stats else 1
        queue.put(1)
    chr_stats['plus_reads'] = len(read_names['+'])
    chr_stats['minus_reads'] = len(read_names['-'])
    chr_stats['avg_dp'] = np.mean(avg_dp) if len(avg_dp)>0 else 0
    return (chrom+"_"+str(start)+"_"+str(end), vcf_records, tsv_records, uncalled_pos, debug_signal, chr_stats)

# .........................................................................
def load_truth_set(config, s):
    """ Load truth set VCF """
    if 'truth_vcf' in config:
        truth_vcf_file = config['truth_vcf']
        # demultiplex samples.
        if ('demultiplex' in config) and (s.name in config['truth_vcf']):
            truth_vcf_file = config['truth_vcf'][s.name]
        if truth_vcf_file is not None:
            reader = vcfpy.Reader.from_path(truth_vcf_file)
            truth_records={}
            truth_filtered={}
            truth_stats={}
            for record in reader:
                if is_pass_rec(record):
                    var_type = get_var_type(record)
                    truth_stats[var_type] = truth_stats[var_type] + 1 if var_type in truth_stats else 1
                    truth_records[rec_to_str(record)]=record
                else:
                    truth_filtered[rec_to_str(record)]=record
            logging.info("Loaded %i records from truth vcf %s: %s" % (len(truth_records), truth_vcf_file, truth_stats) )        
            return truth_records, truth_filtered
        else:
            logging.warn("No truth vcf configured for sample %s !" % (s.name))
            return None,None
    else:
        return None,None

# .........................................................................
#                CALL variants
# ......................................................................... 
def call_variants(config, samples, outdir):
    """ call variants """
    logging.info("-----------------------------------------------")
    logging.info("Call variants")
    logging.info("-----------------------------------------------")
    config = validate_config(config)
    with open(outdir + config["dataset_name"] + ".effective_config.json", 'w') as out:
        print(commentjson.dumps(config, indent=4, sort_keys=True), file=out)

    # .........................................................................
    #                load ref seq and ROI
    # .........................................................................       
    ref_fasta = pysam.FastaFile( config["ref"] ) # @UndefinedVariable
    if 'roi_file' in config:
        roi=pr.read_bed(config["roi_file"]).merge().df # merge overlapping intervals
    elif 'roi_intervals' in config:
        rois=[]
        for ri in config['roi_intervals']: # TODO: check that this is a list!
            chrom = ri.split(":")[0]
            start = int(ri.split(":")[1].split("-")[0])
            end = int(ri.split(":")[1].split("-")[1])
            rois +=[pr.PyRanges(chromosomes=[chrom], starts=[start], ends=[end])]
        roi = pr.concat( rois ).df
    else:
        roi = pr.concat( [pr.PyRanges(chromosomes=[x], starts=[0], ends=[y]) for x,y in zip(ref_fasta.references, ref_fasta.lengths)] ).df
    total_roi_width = sum(roi.End-roi.Start) 
    logging.info("Calling variants in the following regions:")
    logging.info(roi)


    for m in config['mappers'].keys():
        for s in samples:
            if files_exist([s.vcf_file[m]]):
                logging.warning("Found existing VCF file %s. Will not re-create..." % s.vcf_file[m])
            else:
                logging.info("Calling %s / %s" % ( s.name, m ) )
                # .........................................................................
                #                Truth set
                # .........................................................................
                truth_records, truth_filtered = load_truth_set(config, s)

                # .........................................................................
                #                Parallel calling
                # .........................................................................
                try:
                    with Pool(processes=config['threads']) as pool:
                        queue = Manager().Queue()
                        progress = Process(target=show_prog, args=(queue, total_roi_width))
                        progress.start()
                        param = zip(
                            itertools.repeat(s.name),
                            list(roi.Chromosome),  
                            list(roi.Start),
                            list(roi.End), 
                            itertools.repeat(s.downsampled_bam_file[m]),
                            itertools.repeat(config),
                            itertools.repeat(queue)
                            )
                        results = pool.map(call_variants_in_roi, param)
                        logging.info("\nreceiving results...")
                        progress.terminate()
                        data={}
                        for c in results:
                            data[c[0]]= c
                except:
                    logging.info("Terminating progress thread.")
                    progress.terminate()
                    raise
            
                
                # .........................................................................
                #                Call variants
                # .........................................................................
                uncalled_pos=[] # uncalled positions due to low coverage
                debug_signal={} # debug signal
                vcf_records=OrderedDict() # VCF records
                tsv_records=OrderedDict() # TSV records
                stats=OrderedDict()
                for _, row in roi.iterrows():
                    chrom=row['Chromosome']
                    start=row['Start']
                    end=row['End']
                    rid = chrom+"_"+str(start)+"_"+str(end)
                    rid, vcf_records_chr, tsv_records_chr, uncalled_pos_chr, debug_signal_chr, chr_stats = data[rid]
                    logging.info("Got %i vcf/%i tsv/%i uncalled records for roi %s" % (len(vcf_records_chr), len(tsv_records_chr), len(uncalled_pos_chr), rid) )
                    for rec in vcf_records_chr:
                        vcf_records[rec_to_str(rec)] = rec
                    for rec in tsv_records_chr:
                        tsv_records[tsv_to_str(rec)] = rec
                    uncalled_pos += uncalled_pos_chr
                    for k in debug_signal_chr.keys():
                        debug_signal[k] = debug_signal[k] + debug_signal_chr[k] if k in debug_signal else debug_signal_chr[k]
                    stats[chrom]=chr_stats
                        
                # check truth status
                if truth_records:
                    with open(s.eval_file[m], 'w') as out_eval:
            
                        # annotate existing records
                        for k,v in vcf_records.items():
                            status = ("FN" if k in truth_records else "TN") if not is_pass_rec(v) else ("TP" if k in truth_records else "FP")
                            if k in truth_records and 'AF' in truth_records[k].INFO:
                                v.INFO['VAF_EXP'] = ( get_str( truth_records[k].INFO['AF'] ) )
                            v.INFO['STATUS'] = ( status )
                        
                        performance=OrderedDict()
                        for x in ["TP", "TN", "FP", "FN"]:
                            for t in ["SNV", "DEL", "INS"]:
                                performance[x+"_"+t]=0
            
                        # add missing calls as FN 
                        for x in truth_records.keys():
                            if x not in tsv_records and in_roi(truth_records[x], roi):
                                performance["FN"+"_"+get_var_type(truth_records[x]) ]+=1
                                print("FN %s" % (x), file = out_eval )
                                tsv_records[x] = create_missing_record( truth_records[x] )
                        
                        # TSV    
                        for k,v in tsv_records.items():
                            status = ("FN" if k in truth_records else "TN") if not is_pass_tsv(v) else ("TP" if k in truth_records else "FP")
                            vaf_exp = "NA"
                            if k in truth_records and 'AF' in truth_records[k].INFO:
                                vaf_exp = get_str(truth_records[k].INFO['AF'])
                            filtered_in_truth = "1" if ((status=='FP') and (k in truth_filtered)) else 0
                            v += [ status, vaf_exp, filtered_in_truth ]
                            performance[status+"_"+v[4]]+=1
                    
                        print("------------ truth stats -------------", file = out_eval)
                        for t in ["SNV", "DEL", "INS"]:
                            for x in ["TP", "TN", "FP", "FN"]:
                                print("%s\t%s\t%i" % (t, x, performance[x+"_"+t]), file = out_eval )
                                
                
                # write stats
                with open(s.stats_file[m], 'w') as out_eval:
                    print("chr\tvariable\tvalue", file = out_eval )
                    for chrom in stats.keys():
                        for k in stats[chrom].keys():
                            print("%s\t%s\t%s" % (chrom, k, str(stats[chrom][k])), file = out_eval )
                            
                # write results
                # .........................................................................
                #                Create VCF
                # .........................................................................
                header = vcfpy.Header(samples=vcfpy.SamplesInfos([config['dataset_name']+"." + s.name] if config['add_gt_field'] else []))
                header.add_line(vcfpy.HeaderLine("fileformat", "VCFv4.2"))
                header.add_line(vcfpy.HeaderLine("fileDate", datetime.datetime.today().strftime('%Y-%m-%d')))
                header.add_line(vcfpy.HeaderLine("source", "nanopanel2"))
                header.add_line(vcfpy.HeaderLine("reference", config['ref']))
                
                # INFO fields
                header.add_info_line(vcfpy.OrderedDict([('ID', 'AF'), ('Number', '1'), ('Type', 'Float'), ('Description', 'Raw vaf')]))
                header.add_info_line(vcfpy.OrderedDict([('ID', 'VAF_CORR'), ('Number', '1'), ('Type', 'Float'), ('Description', 'Corrected vaf')]))
                header.add_info_line(vcfpy.OrderedDict([('ID', 'CNT_RAW'), ('Number', '1'), ('Type', 'Float'), ('Description', 'Raw AA count')]))
                header.add_info_line(vcfpy.OrderedDict([('ID', 'CNT_CORR'), ('Number', '1'), ('Type', 'Float'), ('Description', 'Corrected AA count')]))
                header.add_info_line(vcfpy.OrderedDict([('ID', 'SB_P'), ('Number', '1'), ('Type', 'Float'), ('Description', 'Strand bias Chi2 P-value')]))
                header.add_info_line(vcfpy.OrderedDict([('ID', 'AA_SKEW'), ('Number', '1'), ('Type', 'Float'), ('Description', 'Alt-allele skew')]))
                header.add_info_line(vcfpy.OrderedDict([('ID', 'HP_FRAC'), ('Number', '1'), ('Type', 'Float'), ('Description', 'Homopolymer fraction')]))
                header.add_info_line(vcfpy.OrderedDict([('ID', 'HP_LEN'), ('Number', '1'), ('Type', 'Float'), ('Description', 'Homopolymer length')]))
                header.add_info_line(vcfpy.OrderedDict([('ID', 'TYPE'), ('Number', '1'), ('Type', 'String'), ('Description', 'Variant type')]))
                header.add_info_line(vcfpy.OrderedDict([('ID', 'STATUS'), ('Number', '1'), ('Type', 'String'), ('Description', 'Truth status')]))
                header.add_info_line(vcfpy.OrderedDict([('ID', 'BASE_QUAL'), ('Number', '1'), ('Type', 'String'), ('Description', 'Average quality for the alt-allele (SNV and INDEL only)')]))
                header.add_info_line(vcfpy.OrderedDict([('ID', 'VAF_EXP'), ('Number', '1'), ('Type', 'Float'), ('Description', 'Expected allele frequency')]))
                header.add_info_line(vcfpy.OrderedDict([('ID', 'CTX'), ('Number', '1'), ('Type', 'String'), ('Description', 'Sequence context')]))
                header.add_info_line(vcfpy.OrderedDict([('ID', 'DEBUG'), ('Number', '1'), ('Type', 'String'), ('Description', 'Debug info')]))
                header.add_info_line(vcfpy.OrderedDict([('ID', 'SOMATIC'), ('Number', '1'), ('Type', 'Flag'), ('Description', 'Set to indicate somatic variants')]))
                
                # for consensus calling
                header.add_info_line(vcfpy.OrderedDict([('ID', 'PASS_CALLS'), ('Number', '1'), ('Type', 'Float'), ('Description', 'Number of pass calls in consensus calling mode')]))
                header.add_info_line(vcfpy.OrderedDict([('ID', 'FILTER_CALLS'), ('Number', '1'), ('Type', 'Float'), ('Description', 'Number of filtered calls in consensus calling mode')]))
                header.add_info_line(vcfpy.OrderedDict([('ID', 'ALL_CALLS'), ('Number', '1'), ('Type', 'Float'), ('Description', 'Number of all calls in consensus calling mode')]))
                
                # filters
                header.add_filter_line(vcfpy.OrderedDict([("ID", "AF"), ("Description", "Allele frequency too low")])) # NOTE this overwrites the previous 'AF' entry leading to a warning from vcfpy. ignore.
                header.add_filter_line(vcfpy.OrderedDict([("ID", "DP"), ("Description", "Coverage too low")]))
                header.add_filter_line(vcfpy.OrderedDict([("ID", "BQ"), ("Description", "Low per-base quality filter")]))
                header.add_filter_line(vcfpy.OrderedDict([("ID", "SB"), ("Description", "Strand bias")]))
                header.add_filter_line(vcfpy.OrderedDict([("ID", "AQ1"), ("Description", "Alt-allele quality filter 1")]))
                header.add_filter_line(vcfpy.OrderedDict([("ID", "AQ2"), ("Description", "Alt-allele quality filter 2")]))
                header.add_filter_line(vcfpy.OrderedDict([("ID", "HP"), ("Description", "Homopolymer filter")]))
                header.add_filter_line(vcfpy.OrderedDict([("ID", "LQ"), ("Description", "Loq overall quality filter")]))
                header.add_filter_line(vcfpy.OrderedDict([("ID", "SI"), ("Description", "Strong strand imbalance")]))
                header.add_filter_line(vcfpy.OrderedDict([("ID", "CO"), ("Description", "Not enough PASS calls for consensus")]))
                
                # format + samples
                header.add_format_line(vcfpy.OrderedDict([('ID', 'DP'), ('Number', '1'), ('Type', 'Integer'), ('Description', 'Depth')]))
                header.add_format_line(vcfpy.OrderedDict([('ID', 'GT'), ('Number', '1'), ('Type', 'String'), ('Description', 'Genotype')]))
                
                writer = vcfpy.Writer.from_path(s.vcf_file[m], header)
                for k,v in vcf_records.items():
                    writer.write_record(v)
            
                # .........................................................................
                #                Create TSV
                # .........................................................................
                with open(s.tsv_file[m], 'w') as out_tsv:
                    headers = TSV_HEADERS.copy()
                    if truth_records:
                        headers += ["truth_status", "vaf_exp", "filtered_in_truth"]
                    print("\t".join(headers), file = out_tsv)
                    for k,v in tsv_records.items():  
                        print('\t'.join([str(x) for x in v]), file = out_tsv)
            
                writer.close()
            
                # write BED file with uncalled positions
                if uncalled_pos:
                    with open(outdir + config["dataset_name"]+'.'+s.name +"."+m+".low_coverage.bed", "w") as out:
                        last_p = None
                        for p in uncalled_pos:
                            if last_p:
                                if last_p[0]==p[0] and last_p[2]==p[1]:
                                    p[1]=last_p[1]
                                else:
                                    print("%s\t%i\t%i\tlow_coverage" % (p[0], p[1], p[2]), file=out )
                                    p = None
                            last_p = p
                        if p:
                            print("%s\t%i\t%i\tlow_coverage" % (p[0], p[1], p[2]), file=out )
            
                # write bedgraph file with debug signals
                if config['debug_mode'] and debug_signal:
                    for k in debug_signal_chr.keys():
                        with open(outdir + config["dataset_name"]+'.'+s.name+"."+m+".debug."+k+".bedGraph", "w") as out:
                            for p in debug_signal[k]:
                                chrom, pos, sig = p
                                print("%s\t%i\t%i\t%f" % (chrom, pos, pos+1, sig), file=out )

    # compress and index vcf
    #bgzip(vcf_file, delinFile=True, index=True, override=True, exe = get_exe("bgzip", config))
    #vcf_file += ".gz"
    return samples


# ###############################################################################
#                Calculate consensus VCF 
# ###############################################################################
class VCFPOS():
    """ Unique position in a VCF file """ 
    def __init__(self, rec ):
        self.pos = rec.POS
        self.ref = rec.REF
        self.alt =  ",".join([alt.value if alt.value else "." for alt in rec.ALT])
    def __repr__(self):
        return "%i %s>%s" % (self.pos, self.ref, self.alt)
    def __eq__(self, other):
        if not isinstance(other, VCFPOS):
            return False
        return ((self.pos == other.pos) and (self.ref == other.ref)and (self.alt == other.alt))
    def __hash__(self):
        return hash(self.__repr__())

def consensus_call(calls, truth_records):
    """ Calculate a consensus (majority vote) call from the passed calls """
    nonull = [x for x in calls if (x is not None)]
    
    if len(nonull)==0:
        return None
    
    INFO=vcfpy.OrderedDict()
    INFO['SOMATIC']=1
    # calculate mean attributes
    for att in ['AF', 'CNT_RAW', 'SB_P', 'AA_SKEW', 'HP_FRAC', 'HP_LEN', 'VAF_CORR', 'CNT_CORR', 'BASE_QUAL']:
        INFO[att] =   (str(np.mean([float(x.INFO[att]) for x in nonull if (att in x.INFO and x.INFO[att] != 'None')])))
    # copy attributes from 1st entry
    for att in ['TYPE', 'CTX']:
        INFO[att] = nonull[0].INFO[att]
    # mean qual
    quality = np.mean([x.QUAL for x in nonull])   
    
    # filters
    pass_calls = len([x for x in nonull if (is_pass_rec(x))])
    filter_calls = len([x for x in nonull if (not is_pass_rec(x))])
    all_calls = len(calls)
    INFO['PASS_CALLS'] = (str(pass_calls))
    INFO['FILTER_CALLS'] = (str(filter_calls))
    INFO['ALL_CALLS'] = (str(all_calls))
    filters=[] if pass_calls/all_calls >= 0.5 else ['CO']
    
    # format
    avg_dp = int( np.mean([float(x.calls[0].data["DP"]) for x in nonull if len(x.FORMAT)>0]) )
    
    # calls 
    vcf_calls = [vcfpy.Call(config['dataset_name']+'.consensus', vcfpy.OrderedDict([('GT', '1'),('DP', str(avg_dp))]))] if config['add_gt_field'] else None

    
    #build vcf record
    vcf_record = vcfpy.Record(CHROM=nonull[0].CHROM,
                              POS=nonull[0].POS,
                              ID=[], 
                              REF=nonull[0].REF,
                              ALT=nonull[0].ALT, 
                              QUAL=quality, 
                              FILTER=filters, 
                              INFO=INFO, 
                              FORMAT=["GT", "DP"] if config['add_gt_field'] else [],
                              calls=vcf_calls )
    if config['add_gt_field']:
        vcf_record.add_format("DP", str(avg_dp))
    
    # calc (optional) truth status
    if truth_records:
        k = rec_to_str(vcf_record)
        v = vcf_record
        status = ("FN" if k in truth_records else "TN") if not is_pass_rec(v) else ("TP" if k in truth_records else "FP")
        if k in truth_records and 'AF' in truth_records[k].INFO:
            v.INFO['VAF_EXP'] = ( get_str( truth_records[k].INFO['AF'] ) )
        v.INFO['STATUS'] = ( status )

    return vcf_record    

def create_consensus_tsv(v, truth_records, truth_filtered ):
    """ Creates a consensus TSV file """
    r = [ 
        v.CHROM,
        v.POS,
        v.REF,
        ",".join([alt.value if alt.value else "." for alt in v.ALT]),
        v.INFO['TYPE']
        ]
    r+=[0 if is_pass_rec(v) else 1]
    for fil in ["AF","DP","BQ","SB","AQ1","AQ2","HP","LQ","SI","CO"]:
        r+=[1 if fil in v.FILTER else 0]
    r+=[no_none(v.QUAL)]
    for att in ['AF', 'CNT_RAW', 'VAF_CORR', 'CNT_CORR', 'AA_SKEW', 'SB_P', 'HP_FRAC', 'HP_LEN',  'BASE_QUAL', 'CTX']:
        r+=[no_none(v.INFO[att])]
        
    if truth_records:
        k = rec_to_str(v) 
        status = ("FN" if k in truth_records else "TN") if not is_pass_rec(v) else ("TP" if k in truth_records else "FP")
        vaf_exp = "NA"
        if k in truth_records and 'AF' in truth_records[k].INFO:
            vaf_exp = get_str(truth_records[k].INFO['AF'])
        filtered_in_truth = "1" if ((status=='FP') and (k in truth_filtered)) else 0
        r += [ status, vaf_exp, filtered_in_truth ]
    return r
            

def calc_consensus(config, samples, outdir):
    """ calculate consensus VCF  """
    logging.info("-----------------------------------------------")
    logging.info("Calculate consensus vcf")
    logging.info("-----------------------------------------------")
    for s in samples:
        logging.info("Calculating consensus for %s" % ( s.name ) )
        chrs = pysam.FastaFile( config['ref'] ).references # @UndefinedVariable

        # .........................................................................
        #                Truth set
        # .........................................................................
        truth_records, truth_filtered = load_truth_set(config, s)
        
        # read all records
        all_rec={}
        all_pos={}
        for c in chrs:
            all_rec[c] = {}
            all_pos[c] = set()
        for m in config['mappers']:
            for c in chrs:
                all_rec[c][m]={}
            reader = vcfpy.Reader.from_path(s.vcf_file[m])
            for r in reader:
                p = VCFPOS(r)
                all_rec[r.CHROM][m][p]=r
                all_pos[r.CHROM].add(p)
            reader.close()
        # sort all positions    
        for c in chrs:
            all_pos[c] = sorted(all_pos[c], key = lambda x: (x.pos, x.ref, x.alt))
        
        # write VCF + TSV
        out_vcf_file = outdir + config["dataset_name"]+( "."+s.name if s.barcode else "" ) + ".consensus.vcf"
        out_tsv_file = outdir + config["dataset_name"]+( "."+s.name if s.barcode else "" ) + ".consensus.tsv"
        with open(out_tsv_file, 'w') as out_tsv:
            headers = TSV_HEADERS.copy()
            if truth_records:
                headers += ["truth_status", "vaf_exp", "filtered_in_truth"]
            print("\t".join(headers), file = out_tsv)
            vcf_header = reader.header.copy()
            if config['add_gt_field']:
                vcf_header.samples = vcfpy.SamplesInfos([config['dataset_name']+".consensus"])
            writer = vcfpy.Writer.from_path(out_vcf_file, vcf_header)
            for c in chrs:
                for p in all_pos[c]:
                    cc = consensus_call([all_rec[c][m][p] if p in all_rec[c][m] else None for m in config['mappers']], truth_records)
                    writer.write_record(cc)
                    print('\t'.join([str(x) for x in create_consensus_tsv(cc, truth_records, truth_filtered)]), file = out_tsv)       
            writer.close()
                 
    logging.info("Done.")
    
# ###############################################################################
#                Post-filter variants 
# ###############################################################################
def post_filter(vcf_file, filter_config, out_file):
    """ Post-filtering of variants """
    if not files_exist(vcf_file):
        sys.exit("Could not find nanopanel2 output vcf file at %s" % (vcf_file) )
    FILTERS = ['AF','CNT_RAW','SB_P','AA_SKEW','HP_FRAC','HP_LEN','TYPE','VAF_CORR','CNT_CORR','BASE_QUAL']
    if 'roi_file' in filter_config:
        roi=pr.read_bed(filter_config["roi_file"]).merge().df # merge overlapping intervals
    elif 'roi_intervals' in filter_config:
        rois=[]
        for ri in filter_config['roi_intervals']: # TODO: check that this is a list!
            chrom = ri.split(":")[0]
            start = int(ri.split(":")[1].split("-")[0])
            end = int(ri.split(":")[1].split("-")[1])
            rois +=[pr.PyRanges(chromosomes=[chrom], starts=[start], ends=[end])]
        roi = pr.concat( rois ).df
    else:
        roi = None

    reader = vcfpy.Reader.from_path(vcf_file)
    writer = vcfpy.Writer.from_path(out_file, reader.header)
    n_filtered = 0
    for r in reader:
        filtered = False
        # filter min/max values
        for f in FILTERS:
            if f+"_min" in filter_config and f in r.INFO and r.INFO[f]<filter_config[f+"_min"]:
                filtered = True
            if f+"_max" in filter_config and f in r.INFO and r.INFO[f]>filter_config[f+"_max"]:
                filtered = True
        # filter roi
        if roi is not None and not in_roi(r, roi):
            filtered = True
        if filtered:
            n_filtered += 1
        else:
            writer.write_record(r)
    reader.close()
    writer.close()
    logging.info("Done. Filtered %i variants" % (n_filtered))
    print("Done. Filtered %i variants" % (n_filtered))

# ###############################################################################
#                Calculate haplotypes 
# ###############################################################################
def calc_all_haplotypes(config, samples):
    """ Calculate all haplotypes for all passed samples """
    for s in samples:
        for m in config['mappers']:
            calc_haplotypes (s.vcf_file[m], s.bam_file[m], s.hap_file_prefix[m], config['max_dp'], config['min_base_quality'], m)
            
def calc_haplotypes(vcf_file, bam_file, hap_file_prefix, max_dp = DEF_CONF['max_dp'], min_base_quality = DEF_CONF['min_base_quality'], mapper=None):
    """ Calculate haplotypes for the passed sample """

    logging.info("-----------------------------------------------")
    logging.info("Calculating haplotypes")
    logging.info("-----------------------------------------------")

    pdf_out_file = hap_file_prefix+'.haplotypes.pdf'
    
    if files_exist([pdf_out_file]):
        logging.warning("Found existing PDF file %s. Will not re-create..." % pdf_out_file)
    else:
        logging.info("Calculating haplotype map %s " % ( pdf_out_file ) )
    
        # create output dir if non existing
        haplo_dir = pathlib.Path(hap_file_prefix).parent
        if not os.path.exists(haplo_dir):
            print("Creating dir " + str(haplo_dir))
            os.makedirs(haplo_dir)
    
        # read variants
        variants=OrderedDict()
        reader = vcfpy.Reader.from_path(vcf_file)
        for rec in reader:
            if is_pass_rec(rec):
                if rec.CHROM not in variants:
                    variants[rec.CHROM] = OrderedDict()
                pos = rec.CHROM + ":" + str(rec.POS)
                variants[rec.CHROM][pos] = variants[rec.CHROM][pos] + [rec] if ( pos in variants[rec.CHROM] ) else [rec]    
        bam = pysam.AlignmentFile(bam_file, "rb") # @UndefinedVariable
    
        # for storing readname / variant / support data        
        read2var=OrderedDict()
        read2strand=OrderedDict()
        # @see https://pysam.readthedocs.io/en/latest/api.html
        flag_filter = BAM_FUNMAP | BAM_FSECONDARY | BAM_FQCFAIL | BAM_FDUP | BAM_SUPPLEMENTARY
        bam_iter = bam.pileup(contig=None, start=None, stop=None, region=None, flag_filter=flag_filter, max_depth=max_dp, truncate=True, mark_matches=False, mark_ends=False, add_indels=False, min_base_quality=min_base_quality)
        for pileup_column in bam_iter:
            chrom = pileup_column.reference_name
            ref_pos1=pileup_column.reference_pos+1
            pos = chrom+":"+str(ref_pos1)
            if (chrom not in variants) or (pos not in variants[chrom].keys()):
                continue
            
            if chrom not in read2var:
                read2var[chrom]=OrderedDict()
                read2strand[chrom]=OrderedDict()
            # iterate reads and collect stats
            for _, rec in enumerate(pileup_column.pileups):
                read_name = rec.alignment.query_name
                pos_in_read = rec.query_position
                read_is_del = rec.is_del or rec.is_refskip
                read_strand = "-" if rec.alignment.is_reverse else "+"  
                for v in variants[chrom][pos]:
                    ref, alt = v.REF, [alt.value if alt.value else "." for alt in v.ALT][0] # TODO: only 1st alt allele considered
                    vtype = 'SNV' if (len(alt)==len(ref) and len(alt)==1) else 'DEL' if (len(alt)<len(ref)) else 'INS' # TODO: MNVs not considered
    
                    # read overlaps variant!
                    if read_name not in read2var[chrom]:
                        read2var[chrom][read_name]={}
                    read2strand[chrom][read_name]=read_strand
                    
                    # does this read support the variant?
                    if vtype == 'INS' and rec.indel>0:
                        ins_bases = rec.alignment.query_sequence[pos_in_read:pos_in_read+rec.indel+1]
                        if ins_bases == alt:
                            read2var[chrom][read_name][rec_to_str(v)] = 1
                        elif ins_bases == ref:
                            read2var[chrom][read_name][rec_to_str(v)] = 0
                    if vtype == 'DEL' and rec.indel<0:    
                        if abs(rec.indel) == len(ref)-1:
                            read2var[chrom][read_name][rec_to_str(v)] = 1
                        else:
                            read2var[chrom][read_name][rec_to_str(v)] = 0
                    if vtype == 'SNV' and not read_is_del:
                        read_base = rec.alignment.query_sequence[pos_in_read]
                        if read_base == alt:
                            read2var[chrom][read_name][rec_to_str(v)] = 1
                        elif read_base == ref:
                            read2var[chrom][read_name][rec_to_str(v)] = 0
    
        medians={}
        for t in ["TP","FP","TN","FN"]:
            medians[t]=[]
        
        # print per haplotypes
        plots=[]
        for chrom in variants.keys():
            with open(hap_file_prefix+'.' + chrom + '.read_info', 'w') as out:
                print("read_name\tstrand\t%s" % ("\t".join([rec_to_str(v)+'|'+ (v.INFO['STATUS'] if 'STATUS' in v.INFO else 'NA') for _, v in enumerate([x for v in variants[chrom].values() for x in v])])), file=out )
                for read_name in read2var[chrom].keys():
                    bit_map = [ str(read2var[chrom][read_name][rec_to_str(x)]) if rec_to_str(x) in read2var[chrom][read_name] else '.' for v in variants[chrom].values() for x in v  ]
                    print("%s\t%s\t%s" % (read_name, read2strand[chrom][read_name], "\t".join(bit_map)), file = out)
    
            # count haplotype support
            ht2count = {}
            ht2count_plus = {}
            for read_name in read2var[chrom].keys():
                bit_map = "\t".join([ str(read2var[chrom][read_name][rec_to_str(x)]) if rec_to_str(x) in read2var[chrom][read_name] else '.' for v in variants[chrom].values() for x in v  ])
                ht2count[bit_map] = ht2count[bit_map] + 1 if bit_map in ht2count else 1
                if read2strand[chrom][read_name] == '+':
                    ht2count_plus[bit_map] = ht2count_plus[bit_map] + 1 if bit_map in ht2count_plus else 1
            # sort by counts
            ht2count = {k: v for k, v in sorted(ht2count.items(), key=lambda item: item[1], reverse=True)}
            max_count = ht2count[next(iter(ht2count))]
            y_axis_labels = [rec_to_str(v)+'|'+ (str(round(v.INFO['AF'], 3)) if 'AF' in v.INFO else 'NA') for _, v in enumerate([x for v in variants[chrom].values() for x in v])]
            y_axis_colors = [('green' if v.INFO['STATUS'] == 'TP' else 'red') if 'STATUS' in v.INFO else 'grey' for _, v in enumerate([x for v in variants[chrom].values() for x in v])]
            x_axis_labels = []
            
            plot_data=[]
            max_haplotypes = 25
            with open(hap_file_prefix+'.' + chrom + '.counts', 'w') as out:
                print("ht\tcounts\tstrand_skew\t%s" % ("\t".join([rec_to_str(v)+'|'+ (str(round(v.INFO['AF'], 3)) if 'AF' in v.INFO else 'NA') for _, v in enumerate([x for v in variants[chrom].values() for x in v])])), file=out )
                i,shown_ht=0,0
                for ht, count in ht2count.items():
                    i+=1
                    plus_count = ht2count_plus[ht] if ht in ht2count_plus else 0
                    minus_count = count-plus_count
                    strand_skew = math.log2(plus_count/minus_count) if ( minus_count>0 and plus_count>0) else None
                    print("ht_%03d\t%i\t%s\t%s" % (i, count, no_none(strand_skew), ht), file=out)
                    if (strand_skew is not None) and (abs(strand_skew)<2) and (count/max_count>0.001) and (count>10) and i<=max_haplotypes:
                        shown_ht+=1
                        x_axis_labels += ["ht_%03d: %.3f" % (i, round(count/len(read2var[chrom]),3))]
                        plot_data += [[save_log(int(x) * count, 0) if (x != '.') else np.nan for x in ht.split('\t')]]
            more_than_max_hp = shown_ht>max_haplotypes   
            # print heatmap
            if ( len(plot_data) > 0 ):
                pal=[(0.95,0.95,0.95)]+sns.color_palette("Blues")
                #print(np.array(plot_data).shape)
                hmap = sns.heatmap( np.rot90(np.fliplr(plot_data)), 
                                   linewidth=0.5, 
                                   xticklabels=x_axis_labels, 
                                   yticklabels=y_axis_labels, 
                                   cmap=pal,
                                   cbar_kws={"shrink": 0.5},
                                   square=True)
                hmap.axes.set_title("%s haplotype map %s %s" % (chrom, '' if mapper is None else ' [%s]' % mapper, '[Only %i of %i HT shown]'%(max_haplotypes,shown_ht) if more_than_max_hp else ''),fontsize=8)
                hmap.tick_params(labelsize=5)
                hmap.figure.subplots_adjust(left = 0.2)
                hmap.figure.subplots_adjust(bottom = 0.1)
                for ytick, color in zip(hmap.axes.get_yticklabels(), y_axis_colors):
                    ytick.set_color(color)
                
                #hmap.figure.savefig(hap_file_prefix+'.' + chrom + '.pdf')  
                plt.close()
                plots+=[hmap]
        # write combined PDF
        with PdfPages(pdf_out_file) as pdf_pages:
            for hmap in plots:
                pdf_pages.savefig(hmap.figure)
            
    logging.info("Done calculating haplotypes.")

            
       
# ###############################################################################
#                Main 
# ###############################################################################
def nanopanel2_pipeline(config, outdir):
    """ Main nanopanel2 pipeline """
    startTime = time.time()
    
    # init logging
    if 'logfile' in config:
        logging.basicConfig(filename=outdir + config['logfile'],level=logging.DEBUG)
    logging.info(LOGO)
    logging.info( "Started: %s" % datetime.datetime.fromtimestamp(startTime).strftime('%Y-%m-%d %H:%M:%S') )
    
    if 'dataset_name' in config:
        logging.info('Analyzing %s ' % (config['dataset_name']) )

    if 'random_seed' in config:
        random.seed(config['random_seed'])

    # demultiplex ?
    demux_index = demultiplex_reads(config, outdir) if 'demultiplex' in config else None

    # extract FASTQ and create index file
    samples = extract_fastq(config, demux_index, outdir)
    
    # map reads with minimap2, last and/or ngmlr
    map_reads(config, samples)
    
    # decorate reads
    decorate_reads(config, samples)
    
    # downsample reads 
    downsample_reads(config, samples)
    
    # call variants
    call_variants(config, samples, outdir)

    # calculate consensus
    if 'consensus' in config:
        calc_consensus(config, samples, outdir)

    # calculate haplotypes
    calc_all_haplotypes( config, samples)
    
    logging.info("Wrote %s" % ( [a.vcf_file for a in samples] ) )
    logging.info("Finished in " + str(datetime.timedelta(seconds=time.time()-startTime)))
    print("All Done.")

if __name__ == "__main__":
    MODES = ['call', 'post_filter', 'build_demux_idx', 'calc_haplotypes']
    #============================================================================ 
    usage = LOGO.replace("VERS", VERSION) + '''                           
    
Nanopanel2 calls somatic variants in Nanopore panel sequencing data. 
    
USAGE: nanopanel2.py MODE'''
    if len(sys.argv) <= 1 or sys.argv[1] in ['-h', '--help']:
        print(usage.replace("MODE", "[" + ",".join(MODES) + "]"))
        sys.exit(1)
    mode = sys.argv[1]  
    if mode not in MODES:
        print("No/invalid mode provided. Please use one of " + ",".join(MODES))
        print(usage)
        sys.exit(1)
    usage = usage.replace("MODE", mode)

    parser = {}
    parser["call"] = MyArgumentParser(usage="nanopanel2.py call -c [config.json] -o [outdir]", description=usage, formatter_class=RawDescriptionHelpFormatter)
    parser["call"].add_argument("-c","--config", type=existing_file, required=True, dest="confF", help="JSON config file")
    parser["call"].add_argument("-o", "--out", required=True, dest="outF", help="Output dir", type=str)

    parser["post_filter"] = MyArgumentParser(usage="nanopanel2.py post_filter -c [infile.vcf] -f  [filters.json] -o [outdir]",description=usage, formatter_class=RawDescriptionHelpFormatter)
    parser["post_filter"].add_argument("-i","--in", type=existing_file, required=True, dest="inF", help="nanopanel2 VCF config file")
    parser["post_filter"].add_argument("-f","--filter", type=existing_file, required=True, dest="filterF", help="filter JSON config file")
    parser["post_filter"].add_argument("-o", "--out", required=True, dest="outF", help="output VCF", type=str)

    parser["build_demux_idx"] = MyArgumentParser(usage="nanopanel2.py build_demux_idx -i [input_demux_dir] -o [output_index_file]", description=usage, formatter_class=RawDescriptionHelpFormatter)
    parser["build_demux_idx"].add_argument("-i","--in", type=existing_file, required=True, dest="inF", help="demux dir (from porechop; will be scanned for BCXX.fasta files)")
    parser["build_demux_idx"].add_argument("-o", "--out", required=True, dest="outF", help="output index file", type=str)
        
    parser["calc_haplotypes"] = MyArgumentParser(usage="nanopanel2.py calc_haplotypes -i [input.vcf] -b [input.bam] -o [output_prefix]", description=usage, formatter_class=RawDescriptionHelpFormatter)
    parser["calc_haplotypes"].add_argument("-i","--in", type=existing_file, required=True, dest="inF", help="vcf file")
    parser["calc_haplotypes"].add_argument("-b","--bam", type=existing_file, required=True, dest="bamF", help="bam file")
    parser["calc_haplotypes"].add_argument("-o", "--out", required=True, dest="outPrefix", help="haplotype output file prefix", type=str)
    

    args = parser[mode].parse_args(sys.argv[2:])
    #============================================================================
    if mode == "call":    
        config = commentjson.load(open(args.confF), object_pairs_hook=OrderedDict)    
        # ensure dirs
        outdir = args.outF
        if not outdir.endswith("/"):
            outdir += "/"
        outdir += config["dataset_name"] + "/"
        if not os.path.exists(outdir):
            print("Creating dir " + outdir)
            os.makedirs(outdir)
        
        nanopanel2_pipeline(config, outdir)
    #============================================================================
    if mode == "post_filter":
        filter_config = commentjson.load(open(args.filterF), object_pairs_hook=OrderedDict)    
        post_filter( args.inF, filter_config, args.outF )
    #============================================================================
    if mode == "build_demux_idx":
        build_demux_idx( args.inF, args.outF )
    #============================================================================
    if mode == "calc_haplotypes":
        calc_haplotypes( args.inF, args.bamF, args.outPrefix )

        
        
        
        