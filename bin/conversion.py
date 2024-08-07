#!/usr/bin/env python

import argparse
import os
import subprocess
from collections import defaultdict

import filter_gtf
import pandas as pd
import pysam
import utils


class Conversion:
    """
    Get conversion for each read and add tags:
    - get a list of sites of the specified conversion type
    - get statistics for all conversion types
    This function derived from the procedure implemented in the NASC-seq pipeline (DOI: 10.1038/s41467-019-11028-9).
    """

    def __init__(self, args):
        self.args = args
        # input files
        self.bam = args.bam
        self.gtf_file = args.gtf
        self.bclist = args.bclist
        self.outdir = args.outdir
        self.qual = args.basequalilty
        self.conversion_type = args.conversion_type

        # set
        gp = filter_gtf.GtfParser(self.gtf_file)
        self.strand = gp.get_id_strand()
        self.cells = utils.read_one_col(self.bclist)

        # output files
        outprefix = os.path.basename(self.bam)[:-4]
        self.outbam = f"{self.outdir}/{outprefix}.PosTag.bam"
        self.outcsv = f"{self.outdir}/{outprefix}.PosTag.csv"

    def run(self):
        # Adding tags and get conversion positions
        self.df_conv = self.add_tags()
        self.df_cover = self.count_read_cover_per_conv_pos()
        self.conv_candidate()

    def run_cmd(self, cmd):
        subprocess.call(" ".join(cmd), shell=True)

    def create_tag(self, d):
        return "".join(["".join(key) + str(d[key]) + ";" for key in d.keys()])[:-1]

    def get_types(self, contype):
        rdict = {"A": "T", "T": "A", "G": "C", "C": "G"}
        ftype = (contype[0].lower(), contype[1].upper())
        rtype = (rdict[contype[0]].lower(), rdict[contype[1]].upper())
        return ftype, rtype

    def check_md(self, md_tag):
        for x in ["A", "C", "G", "T"]:
            if x in md_tag:
                return True
        return False

    def conv_in_read(self, read):
        tc_loc, ag_loc = [], []
        total_content = {"a": 0, "c": 0, "g": 0, "t": 0}
        specific_conversions = {}
        specific_conversions[("c", "A")] = 0
        specific_conversions[("g", "A")] = 0
        specific_conversions[("t", "A")] = 0
        specific_conversions[("a", "C")] = 0
        specific_conversions[("g", "C")] = 0
        specific_conversions[("t", "C")] = 0
        specific_conversions[("a", "G")] = 0
        specific_conversions[("c", "G")] = 0
        specific_conversions[("t", "G")] = 0
        specific_conversions[("a", "T")] = 0
        specific_conversions[("c", "T")] = 0
        specific_conversions[("g", "T")] = 0
        specific_conversions[("a", "N")] = 0
        specific_conversions[("c", "N")] = 0
        specific_conversions[("g", "N")] = 0
        specific_conversions[("t", "N")] = 0

        try:
            refseq = read.get_reference_sequence().lower()
        except UnicodeDecodeError:
            return 0
        except AssertionError:
            return 0

        for base in total_content.keys():
            total_content[base] += refseq.count(base)

        if self.check_md(read.get_tag("MD")):
            ftype, rtype = self.get_types(self.conversion_type)
            for pair in read.get_aligned_pairs(with_seq=True):
                try:
                    if pair[0] is not None and pair[1] is not None and pair[2] is not None:
                        if str(pair[2]).islower() and not read.query_qualities[pair[0]] < self.qual:
                            specific_conversions[(pair[2], read.seq[pair[0]])] += 1
                            if (pair[2], read.seq[pair[0]]) == ftype:
                                tc_loc.append(pair[1])
                            if (pair[2], read.seq[pair[0]]) == rtype:
                                ag_loc.append(pair[1])

                except (UnicodeDecodeError, KeyError):
                    continue

        sc_tag = self.create_tag(specific_conversions)
        tc_tag = self.create_tag(total_content)

        if len(tc_loc) == 0:
            tc_loc.append(0)
        if len(ag_loc) == 0:
            ag_loc.append(0)

        return sc_tag, tc_tag, tc_loc, ag_loc

    def add_tags(self):
        site_depth = defaultdict(int)  ## conversion depth for each site
        save = pysam.set_verbosity(0)
        bamfile = pysam.AlignmentFile(self.bam, "rb")
        header = bamfile.header
        mod_bamfile = pysam.AlignmentFile(self.outbam, mode="wb", header=header, check_sq=False)
        pysam.set_verbosity(save)

        class GeneError(Exception):
            pass

        for read in bamfile.fetch(until_eof=True):
            try:
                ## check read info
                if (not read.has_tag("GX")) or read.get_tag("GX") == "-":
                    continue
                if read.get_tag("GX") not in self.strand:
                    raise GeneError
                if read.get_tag("CB") not in self.cells:
                    continue

                tags = self.conv_in_read(read)
                if tags == 0:
                    mod_bamfile.write(read)
                    continue

                read.set_tag("SC", tags[0], "Z")
                read.set_tag("TC", tags[1], "Z")
                read.set_tag("TL", tags[2])
                read.set_tag("AL", tags[3])
                read.set_tag("ST", self.strand[read.get_tag("GX")])

                if self.strand[read.get_tag("GX")] == "+":
                    locs = tags[2]
                else:
                    locs = tags[3]

                if locs[0] != 0:
                    for _, loc in enumerate(locs):
                        site = f"{read.reference_name}+{loc}"
                        site_depth[site] += 1

                mod_bamfile.write(read)

            except (ValueError, KeyError):
                continue
            except GeneError:
                print("{} is not in gtf file, please check your files.".format(read.get_tag("GX")))
                continue

        bamfile.close()
        mod_bamfile.close()

        if len(site_depth) == 0:  ## if no conversion site detected
            return pd.DataFrame()

        df = pd.DataFrame.from_dict(site_depth, orient="index")
        df.columns = ["convs"]
        return df

    def conv_candidate(self):
        df = self.df_conv
        if df.shape[0] == 0:
            df = pd.DataFrame(columns=["chrom", "pos", "convs", "covers"])
            df = df.set_index(["chrom", "pos"])
        else:
            df = df.reset_index()
            df[["chrom", "pos"]] = df["index"].str.split("+", expand=True)
            df["pos"] = df["pos"].astype(int)
            df = df.set_index(["chrom", "pos"])
            df.drop(["index"], axis=1, inplace=True)
            dep = pd.DataFrame.from_dict(self.df_cover, orient="index")
            dep = dep.reset_index()
            dep = dep.melt(id_vars="index").dropna()
            dep.columns = ["chrom", "pos", "covers"]
            dep = dep.set_index(["chrom", "pos"])
            df = pd.concat([df, dep], axis=1)
        # output
        df.to_csv(self.outcsv)

    def count_read_cover_per_conv_pos(self):
        bamfile = self.outbam
        df = self.df_conv
        cover_of_pos_with_convs = {}
        if df.shape[0] == 0:
            return cover_of_pos_with_convs
        df = df.reset_index()
        df[["chrom", "pos"]] = df["index"].str.split("+", expand=True)
        df["pos"] = df["pos"].astype(int)
        try:
            cmd = f"samtools index {bamfile} 2>&1 "
            subprocess.check_call(cmd, shell=True)
        except subprocess.CalledProcessError:
            cmd = f"samtools index -c {bamfile} 2>&1 "
            subprocess.check_call(cmd, shell=True)
        save = pysam.set_verbosity(0)
        bam = pysam.AlignmentFile(bamfile, "rb")
        pysam.set_verbosity(save)
        contig_locs = df.groupby("chrom")["pos"].apply(list).to_dict()
        for key in contig_locs.keys():
            contig_locs[key] = sorted(contig_locs[key])
            cover_of_pos_with_convs[key] = {}
            for key2 in contig_locs[key]:
                try:
                    cover_of_pos_with_convs[key][key2] = bam.count(key, key2, key2 + 1)
                except ValueError:
                    continue
        bam.close()
        return cover_of_pos_with_convs


def get_opts_conversion():
    parser = argparse.ArgumentParser(description="")
    parser.add_argument("--bam", required=True)
    parser.add_argument("--gtf", required=True)
    parser.add_argument("--bclist", required=True)
    parser.add_argument("--outdir", required=True)
    parser.add_argument(
        "--conversion_type",
        type=str,
        default="TC",
        help="conversion type, TC for dynaseq or CT for m6a",
        required=False,
    )
    parser.add_argument("--basequalilty", type=int, help="min base quality of the read sequence", required=False)

    args = parser.parse_args()
    return args


if __name__ == "__main__":
    args = get_opts_conversion()
    Conversion(args).run()
