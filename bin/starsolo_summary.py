#!/usr/bin/env python

import argparse
<<<<<<< HEAD
import json
=======
import os
>>>>>>> upstream/master
from collections import defaultdict

import pandas as pd
import utils
<<<<<<< HEAD

MAX_CELL = 2 * 10**5


def parse_read_stats(read_stats):
    dtypes = defaultdict(lambda: "int")
    dtypes["CB"] = "object"
    df = pd.read_csv(
        read_stats, sep="\t", header=0, index_col=0, skiprows=[1], dtype=dtypes
    )  # skip first line cb not pass whitelist
    umi_count = df["nUMIunique"]
    df = df.loc[
        :, ["cbMatch", "cbPerfect", "genomeU", "genomeM", "exonic", "intronic", "exonicAS", "intronicAS", "countedU"]
    ]
    s = df.sum()
    # json does not recognize NumPy data types. TypeError: Object of type int64 is not JSON serializable
    valid = int(s["cbMatch"])
    perfect = int(s["cbPerfect"])
    corrected = valid - perfect
    genome_uniq = int(s["genomeU"])
    genome_multi = int(s["genomeM"])
    mapped = genome_uniq + genome_multi
    exonic = int(s["exonic"])
    intronic = int(s["intronic"])
    antisense = int(s["exonicAS"] + s["intronicAS"])
    intergenic = mapped - exonic - intronic - antisense
    counted_uniq = int(s["countedU"])
    data_dict = {
        "Corrected Barcodes": corrected / valid,
        "Reads Mapped To Unique Loci": genome_uniq / valid,
        "Reads Mapped To Multiple Loci": genome_multi / valid,
        "Reads Mapped Uniquely To Transcriptome": counted_uniq / valid,
        "Mapped Reads Assigned To Exonic Regions": exonic / mapped,
        "Mapped Reads Assigned To Intronic Regions": intronic / mapped,
        "Mapped Reads Assigned To Intergenic Regions": intergenic / mapped,
        "Mapped Reads Assigned Antisense To Gene": antisense / mapped,
    }
    for k in data_dict:
        data_dict[k] = utils.get_frac(data_dict[k])

    return umi_count, data_dict


def parse_summary(summary):
    data = utils.csv2dict(summary)
    origin_new = {
        "Number of Reads": "Raw Reads",
        "Reads With Valid Barcodes": "Valid Reads",
        "Sequencing Saturation": "Saturation",
        "Estimated Number of Cells": "Estimated Number of Cells",
        "Fraction of Unique Reads in Cells": "Fraction Reads in Cells",
        "Mean Reads per Cell": "Mean Used Reads per Cell",
        "Median UMI per Cell": "Median UMI per Cell",
        "Median GeneFull_Ex50pAS per Cell": "Median Genes per Cell",
        "Total GeneFull_Ex50pAS Detected": "Total Genes",
    }
    parsed_data = {}
    for origin, new in origin_new.items():
        parsed_data[new] = data[origin]
    frac_names = {"Valid Reads", "Saturation", "Fraction Reads in Cells"}
    for k in frac_names:
        parsed_data[k] = utils.get_frac(parsed_data[k])
    for k in set(origin_new.values()) - frac_names:
        parsed_data[k] = int(parsed_data[k])
    return parsed_data
=======
from __init__ import ASSAY

MAX_CELL = 2 * 10**5
FEATURE_FILE_NAME = "features.tsv.gz"
BARCODE_FILE_NAME = "barcodes.tsv.gz"
MATRIX_FILE_NAME = "matrix.mtx.gz"


class StarsoloSummary:
    def __init__(self, args):
        self.args = args
        barcodes_file = os.path.join(args.filtered_matrix, BARCODE_FILE_NAME)
        self.matrix_file = os.path.join(args.filtered_matrix, MATRIX_FILE_NAME)
        self.cbs = utils.read_one_col(barcodes_file)
        self.stats = {}

    def add_total_genes(self):
        n = 0
        genes = set()
        with utils.openfile(self.matrix_file) as f:
            for line in f:
                n += 1
                if n <= 3:
                    continue
                gene = line.split(" ")[0]
                genes.add(gene)
        self.stats["Total Genes"] = len(genes)

    def parse_read_stats(self):
        dtypes = defaultdict(lambda: "int")
        dtypes["CB"] = "object"
        df = pd.read_csv(
            self.args.read_stats, sep="\t", header=0, index_col=0, skiprows=[1], dtype=dtypes
        )  # skip first line cb not pass whitelist
        rbs = list(df.index)
        umi_count = list(df["nUMIunique"])
        df = df.loc[
            :,
            [
                "cbMatch",
                "cbPerfect",
                "genomeU",
                "genomeM",
                "exonic",
                "intronic",
                "exonicAS",
                "intronicAS",
                "countedU",
                "nUMIunique",
                "nGenesUnique",
            ],
        ]
        s = df.sum()
        # json does not recognize NumPy data types. TypeError: Object of type int64 is not JSON serializable
        valid = int(s["cbMatch"])
        perfect = int(s["cbPerfect"])
        corrected = valid - perfect
        genome_uniq = int(s["genomeU"])
        genome_multi = int(s["genomeM"])
        mapped = genome_uniq + genome_multi
        exonic = int(s["exonic"])
        intronic = int(s["intronic"])
        antisense = int(s["exonicAS"] + s["intronicAS"])
        intergenic = mapped - exonic - intronic - antisense
        counted_uniq = int(s["countedU"])
        data_dict = {
            "Corrected Barcodes": corrected / valid,
            "Reads Mapped To Unique Loci": genome_uniq / valid,
            "Reads Mapped To Multiple Loci": genome_multi / valid,
            "Reads Mapped Uniquely To Transcriptome": counted_uniq / valid,
            "Mapped Reads Assigned To Exonic Regions": exonic / mapped,
            "Mapped Reads Assigned To Intronic Regions": intronic / mapped,
            "Mapped Reads Assigned To Intergenic Regions": intergenic / mapped,
            "Mapped Reads Assigned Antisense To Gene": antisense / mapped,
        }
        for k in data_dict:
            data_dict[k] = utils.get_frac(data_dict[k])
        self.stats.update(data_dict)

        n_cells = len(self.cbs)
        reads_cell = df.loc[self.cbs, "countedU"].sum()
        fraction_reads_in_cells = utils.get_frac(float(reads_cell / counted_uniq))
        mean_used_reads_per_cell = int(reads_cell / n_cells)
        median_umi_per_cell = int(df.loc[self.cbs, "nUMIunique"].median())
        median_genes_per_cell = int(df.loc[self.cbs, "nGenesUnique"].median())
        data_dict = {
            "Estimated Number of Cells": n_cells,
            "Fraction Reads in Cells": fraction_reads_in_cells,
            "Mean Used Reads per Cell": mean_used_reads_per_cell,
            "Median UMI per Cell": median_umi_per_cell,
            "Median Genes per Cell": median_genes_per_cell,
        }
        self.stats.update(data_dict)
        return rbs, umi_count

    def parse_summary(self):
        data = utils.csv2dict(self.args.summary)
        origin_new = {
            "Number of Reads": "Raw Reads",
            "Reads With Valid Barcodes": "Valid Reads",
            "Sequencing Saturation": "Saturation",
        }
        parsed_data = {}
        for origin, new in origin_new.items():
            parsed_data[new] = data[origin]
        frac_names = {"Valid Reads", "Saturation"}
        for k in frac_names:
            parsed_data[k] = utils.get_frac(parsed_data[k])
        for k in set(origin_new.values()) - frac_names:
            parsed_data[k] = int(parsed_data[k])
        self.stats.update(parsed_data)

    def run(self):
        rbs, umis = self.parse_read_stats()
        self.add_total_genes()
        self.parse_summary()
        plot_data = utils.get_umi_count(rbs, umis, self.cbs, self.args.sample)
        utils.write_multiqc(plot_data, args.sample, ASSAY, "umi_count")
        utils.write_multiqc(self.stats, args.sample, ASSAY, "starsolo_summary.stats")
>>>>>>> upstream/master


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Starsolo summary")
    parser.add_argument("--read_stats", help="cellReadsStats file")
<<<<<<< HEAD
    parser.add_argument("--barcodes", help="barcode file")
=======
    parser.add_argument("--filtered_matrix", help="filtered_matrix")
>>>>>>> upstream/master
    parser.add_argument("--summary", help="summary file")
    parser.add_argument("--sample", help="sample name")
    args = parser.parse_args()

<<<<<<< HEAD
    umi_count, data_dict = parse_read_stats(args.read_stats)
    read_stats_file = args.sample + ".scrna.read.stats.json"
    utils.write_json(data_dict, read_stats_file)

    # summary
    data_dict = parse_summary(args.summary)
    summary_file = args.sample + ".scrna.starsolo.stats.json"
    utils.write_json(data_dict, summary_file)

    # UMI count
    umi_count.loc[lambda x: x > 0]
    umi_count = umi_count.sort_values(ascending=False)
    cbs = set(utils.read_one_col(args.barcodes))
    plot_data = {}
    n = len(umi_count)
    first_noncell = n - 1
    for i, bc in enumerate(umi_count.index):
        if bc not in cbs:
            first_noncell = i
            break
    last_cell = 0
    for i in range(min(n - 1, MAX_CELL), -1, -1):
        bc = umi_count.index[i]
        if bc in cbs:
            last_cell = i
            break
    pure = args.sample + ".cells.pure" + f"({first_noncell}/{first_noncell}, 100%)"
    bg = args.sample + ".cells.background"
    plot_data[pure] = {}
    plot_data[bg] = {}
    for i in range(first_noncell):
        plot_data[pure][i + 1] = int(umi_count.iloc[i])

    n_mix = last_cell - first_noncell + 1
    if n_mix != 0:
        n_total = len(cbs)
        n_mix_cell = n_total - first_noncell
        mix_rate = round(n_mix_cell / n_mix * 100, 2)
        mix = args.sample + ".cells.mix" + f"({n_mix_cell}/{n_mix}, {mix_rate}%)"
        plot_data[mix] = {}
        for i in range(first_noncell, last_cell + 1):
            plot_data[mix][i + 1] = int(umi_count.iloc[i])

    for i in range(last_cell + 1, min(MAX_CELL, n), 10):
        plot_data[bg][i + 1] = int(umi_count.iloc[i])
    # do not record every umi count
    for i in range(MAX_CELL, n, 1000):
        plot_data[bg][i + 1] = int(umi_count.iloc[i])

    umi_file = args.sample + ".scrna.umi_count.json"
    with open(umi_file, "w") as f:
        json.dump(plot_data, f)
=======
    StarsoloSummary(args).run()
>>>>>>> upstream/master
