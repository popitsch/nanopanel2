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
    nanopanel2 debug tools
'''


from argparse import ArgumentParser, RawDescriptionHelpFormatter
from collections import OrderedDict, Counter
import logging
import sys, os, commentjson
import pandas as pd
import binascii
from utils import existing_file
import pysam


# Necessary for including python modules from a parent directory
sys.path.append(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))



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

def ascii2prob(s):
    """ Converts an ASCII representation to a probability value """
    return binascii.unhexlify(s)

# .........................................................................
#                extract_guppy_distributions
# ......................................................................... 
def extract_guppy_distributions(bam_file, pos_file, out_file):
    """ call variants """
    logging.info("-----------------------------------------------")
    logging.info("extract_guppy_distributions")
    logging.info("-----------------------------------------------")
    config={}
    config['max_dp']=100000
    config['max_af']=0.1
    config['min_base_quality']=0
    config["tag_pA"]="pA" # BAM TAG used to store pA values
    config["tag_pC"]="pC" # BAM TAG used to store pC values
    config["tag_pG"]="pG" # BAM TAG used to store pG values
    config["tag_pT"]="pT" # BAM TAG used to store pT values


    allpos = pd.read_csv(pos_file, sep='\t')
    allpos=allpos[(allpos['type']=='SNV') & (allpos['vaf_raw']<=config['max_af'])]
    pos={}
    for status in ['TP','TN','FP','FN']:
        pos[status]=allpos[allpos['truth_status']==status]
    n_pos=len(pos['TP'].index)
    for status in ['TN','FP','FN']:
        if len(pos[status])>n_pos:
            pos[status]=pos[status].sample(n=n_pos)
    pos=pd.concat(pos)
    
    with open(out_file, 'w') as out:
        print('\t'.join([
                    'chrom', 
                    'ref_pos1',
                    'ref_base',
                    'alt_base',
                    'truth_status',
                    'vaf_raw',
                    'dp',
                    'p_ref',
                    'p_alt',
                    'p_mfa1',
                    'p_mfa2',
                    'p_mfa3',
                    'p_mfa4']), file=out)
            
        bam = pysam.AlignmentFile(bam_file, "rb") # @UndefinedVariable
        # @see https://pysam.readthedocs.io/en/latest/api.html
        flag_filter = BAM_FUNMAP | BAM_FSECONDARY | BAM_FQCFAIL | BAM_FDUP | BAM_SUPPLEMENTARY
        biter = bam.pileup(flag_filter=flag_filter, max_depth=config['max_dp'], truncate=True, mark_matches=False, mark_ends=False, add_indels=False, min_base_quality=config['min_base_quality'])
        for pileup_column in biter:
            chrom = pileup_column.reference_name
            ref_pos1=pileup_column.reference_pos+1
            p = pos.loc[(pos['chr'] == chrom) & (pos['pos'] == ref_pos1)]
            if p.empty:
                continue
            p=p.iloc[0]
            ref_base = str(p['ref']) # this is the refc at the current position
            alt_base = str(p['alt'])
            truth_status = str(p['truth_status'])
            vaf_raw= str(p['vaf_raw'])
            dp = pileup_column.get_num_aligned()
            #print(chrom, ref_pos1, ref_base, alt_base, dp, truth_status)
            # iterate reads and collect stats\
            read_alleles=Counter()
            
            data={}
            
            for _, rec in enumerate(pileup_column.pileups):
                read_is_del = rec.is_del or rec.is_refskip
                if read_is_del:
                    continue
                read_name = rec.alignment.query_name
                pos_in_read = rec.query_position
                read_base = rec.alignment.query_sequence[pos_in_read]
                if read_base != alt_base:
                    continue
                read_alleles[read_base]+=1
                read_strand = "-" if rec.alignment.is_reverse else "+"
                prob={}
                prob['A'] =ascii2prob(rec.alignment.get_tag(config['tag_pA']))[pos_in_read]
                prob['C'] =ascii2prob(rec.alignment.get_tag(config['tag_pC']))[pos_in_read]
                prob['G'] =ascii2prob(rec.alignment.get_tag(config['tag_pG']))[pos_in_read]
                prob['T'] =ascii2prob(rec.alignment.get_tag(config['tag_pT']))[pos_in_read]
                data[read_name]=prob
            mfaa = [x for x,_ in read_alleles.most_common()]
            for a in ['A','C','G','T']: # add missing alleles
                if a not in mfaa:
                    mfaa+=[a]
            for _,prob in data.items():
                
                print('\t'.join([str(x) for x in [
                    chrom, 
                    ref_pos1,
                    ref_base,
                    alt_base,
                    truth_status,
                    vaf_raw,
                    dp,
                    prob[ref_base],
                    prob[alt_base]] + [prob[x] for x in mfaa]]), file=out)
    print("All done.")
    
if __name__ == "__main__":
    MODES = ['extract_guppy_distributions']
    #============================================================================ 
    usage = '''                           
    
Nanopanel2 debug tools. 
    
USAGE: debug_tools.py MODE'''
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
    parser["extract_guppy_distributions"] = ArgumentParser(usage="nanopanel2.py call -c [config.json] -o [outdir]", description=usage, formatter_class=RawDescriptionHelpFormatter)
    parser["extract_guppy_distributions"].add_argument("-b","--bam", type=existing_file, required=True, dest="bam_file", help="BAM file")
    parser["extract_guppy_distributions"].add_argument("-p","--pos", type=existing_file, required=True, dest="pos_file", help="POS file")
    parser["extract_guppy_distributions"].add_argument("-o","--out", type=str, required=True, dest="out_file", help="OUT file")
    
    args = parser[mode].parse_args(sys.argv[2:])
    #============================================================================
    if mode == "extract_guppy_distributions":    
        extract_guppy_distributions(args.bam_file, args.pos_file, args.out_file)
        
