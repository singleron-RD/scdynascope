#!/bin/env python

import os
import pysam
import re
import pandas as pd
import argparse
import utils
from __init__ import ASSAY


class Substitution:
    """
    ## Features
    - Computes the overall conversion rates in reads and plots a barplot.

    ## Output
    - `{sample}.substitution.txt` Tab-separated table of the overall conversion rates.
    """

    def __init__(self, args):
        self.args = args
        # input
        self.sample = args.sample
        self.outdir = args.outdir
        self.in_bam = args.bam
        # set
        self.index_dict = {}
        self.report_dict = {}
        self.df = pd.DataFrame()
        # output files
        self.outstat = os.path.join(self.outdir, self.sample+'.substitution.tsv')


    def run(self):
        for_base, rev_base, is_forward, is_reverse = self.get_sub_tag()
        self.sub_stat(for_base, rev_base, is_forward, is_reverse)
        utils.write_multiqc(self.report_dict, self.sample, ASSAY, "substitution")


    def get_sub_tag(self):
        save = pysam.set_verbosity(0)
        bamfile = pysam.AlignmentFile(self.in_bam, 'rb')
        pysam.set_verbosity(save)
        is_reverse = {'cA': 0, 'gA': 0, 'tA': 0, 'aC': 0, 'gC': 0,
                      'tC': 0, 'aG': 0, 'cG': 0, 'tG': 0, 'aT': 0, 'cT': 0, 'gT': 0}
        is_forward = {'cA': 0, 'gA': 0, 'tA': 0, 'aC': 0, 'gC': 0,
                      'tC': 0, 'aG': 0, 'cG': 0, 'tG': 0, 'aT': 0, 'cT': 0, 'gT': 0}
        for_base = {'a': 0, 'c': 0, 'g': 0, 't': 0}
        rev_base = {'a': 0, 'c': 0, 'g': 0, 't': 0}
        snp_tags = ['', 'cA', 'gA', 'tA', 'aC', 'gC', 'tC', 'aG', 'cG', 'tG', 'aT', 'cT', 'gT']
        ref_tags = ['', 'a', 'c', 'g', 't']
        for read in bamfile.fetch(until_eof=True):
            try:
                if (not read.has_tag('TC')):
                    continue
                snpmatch = re.match(
                    r'cA(\d+);gA(\d+);tA(\d+);aC(\d+);gC(\d+);tC(\d+);aG(\d+);cG(\d+);tG(\d+);aT(\d+);cT(\d+);gT(\d+);', read.get_tag('SC'), re.M)
                totmatch = re.match(r'a(\d+);c(\d+);g(\d+);t(\d+)', read.get_tag('TC'), re.M)
                if snpmatch and totmatch:
                    if read.is_reverse:
                        for j in range(1, len(ref_tags)):
                            rev_base[ref_tags[j]] += int(totmatch.group(j))
                        for i in range(1, len(snp_tags)):
                            is_reverse[snp_tags[i]] += int(snpmatch.group(i))
                    else:
                        for j in range(1, len(ref_tags)):
                            for_base[ref_tags[j]] += int(totmatch.group(j))
                        for i in range(1, len(snp_tags)):
                            is_forward[snp_tags[i]] += int(snpmatch.group(i))
            except (ValueError, KeyError):
                continue
        bamfile.close()

        return for_base, rev_base, is_forward, is_reverse

    def sub_stat(self, for_base, rev_base, is_forward, is_reverse):
        convertdict = {'a': ['aC', 'aG', 'aT'],
                       'c': ['cA', 'cG', 'cT'],
                       'g': ['gA', 'gC', 'gT'],
                       't': ['tA', 'tC', 'tG']}
        subdict = {'a': 't', 't': 'a', 'c': 'g', 'g': 'c',
                   'aC': 'tG', 'aG': 'tC', 'aT': 'tA',
                   'cA': 'gT', 'cG': 'gC', 'cT': 'gA',
                   'gA': 'cT', 'gC': 'cG', 'gT': 'cA',
                   'tA': 'aT', 'tC': 'aG', 'tG': 'aC'}
        outdict = {'aC': 'A_to_C', 'aG': 'A_to_G', 'aT': 'A_to_T',
                   'cA': 'C_to_A', 'cG': 'C_to_G', 'cT': 'C_to_T',
                   'gA': 'G_to_A', 'gC': 'G_to_C', 'gT': 'G_to_T',
                   'tA': 'T_to_A', 'tC': 'T_to_C', 'tG': 'T_to_G'}

        outw = open(self.outstat, 'w')
        for x in ['a', 'c', 'g', 't']:
            fbase = for_base[x]
            rbase = rev_base[subdict[x]]
            for y in convertdict[x]:
                fcov = is_forward[y]*100 / float(fbase) if float(fbase)>0 else 0
                rcov = is_reverse[subdict[y]]*100 / float(rbase) if float(rbase)>0 else 0
                self.report_dict[outdict[y]] = fcov + rcov
                #outw.write(outdict[y]+'\t'+"%.3f" % fcov+'\t'+"%.3f" % rcov+'\n')
                outw.write(f"{outdict[y]}\t{fcov:.3f}\t{rcov:.3f}\n")

        outw.close()



def get_opts_substitution():
    parser = argparse.ArgumentParser(description='')
    parser.add_argument('--outdir', help='Output diretory.', required=True)
    parser.add_argument('--sample', help='Sample name.', required=True)
    parser.add_argument('--bam', help='Required. bam file from step conversion.', required=True)
    args = parser.parse_args()
    return args


if __name__ == '__main__':
    args=get_opts_substitution()
    Substitution(args).run()
