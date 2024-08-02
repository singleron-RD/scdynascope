#!/usr/bin/env python

import pandas as pd
import argparse


class Conversion_merge:
    """
    Merge conversion csv files and select snp candidates:
    """

    def __init__(self, args):
        self.args = args
        # input files    
        self.sample = args.sample
        self.outdir = args.outdir
        self.snp_min_depth = args.snp_min_depth
        self.snp_threshold = args.snp_threshold

        # set
        self.csvlist = args.csvlist.strip().split(',')
        self.df_list = []
        self.df_conv = pd.DataFrame()

        # output files 
        self.outcsv = f'{self.outdir}/{self.sample}.PosTag.csv'
        self.outsnp = f'{self.outdir}/{self.sample}.snp.csv'

    
    def run(self):
        # Adding tags and get conversion positions
        self.read_csv_list()
        self.merge_csv()
        self.snp_candidate()


    def read_csv_list(self):
        for csvfile in self.csvlist:
            df = pd.read_csv(csvfile)
            self.df_list.append(df)

    def merge_csv(self):
        df = pd.concat(self.df_list)
        self.df_conv = df.groupby(['chrom', 'pos']).agg({
                    'convs': 'sum',
                    'covers': 'sum'
                    }).reset_index()
        self.df_conv['posratio'] = self.df_conv['convs'] / self.df_conv['covers']
        self.df_conv.to_csv(self.outcsv, index=False)                

    def snp_candidate(self):
        snp_df = self.df_conv[ self.df_conv['covers'] >= self.snp_min_depth ]
        snp_df = snp_df[ snp_df['posratio'] >= self.snp_threshold]
        snp_df.to_csv(self.outsnp, index=False)



def get_opts_conversion():
    parser = argparse.ArgumentParser(description='')
    parser.add_argument("--csvlist", required=True)
    parser.add_argument("--sample", required=True)
    parser.add_argument('--outdir',  required=True)
    parser.add_argument('--snp_threshold', type=float, 
                        help='snp threshold filter, greater than snp_threshold will be recognized as snp')
    parser.add_argument('--snp_min_depth', type=int,
                        help='Minimum depth to call a variant')
    
    args = parser.parse_args()
    return args



if __name__ == '__main__':
    args=get_opts_conversion()
    Conversion_merge(args).run()