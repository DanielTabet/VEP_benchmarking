# Authors: Roujia Li, Kevin Kuang
# email: Roujia.li@mail.utoronto.ca, kvn.kuang@mail.utoronto.ca
import glob
import os
import gzip
import pandas as pd
import numpy as np
import math
import argparse

# Go through all the output files generated by VEP (*.txt), also go through the raw vcf files
# make two files for each input, *_all_weights and *_filtered_mut
# the purpose is to clean up VEP output files, also clean vcf files to see what variants does each participant carry

# global variables
total_n = 454801
VARITY = "/mnt/project/Caches/varity_all_predictions.txt"
pvcf_block = "/mnt/project/Caches/pvcf_blocks.txt"


class parseVEP(object):

    def __init__(self, vep_output, raw_vcf):
        """
        Read VEP output and vcf file together, to get how many people carrying exon variant
        """
        self._vep_output = vep_output
        self._vcf = raw_vcf
        self._chr = os.path.basename(vep_output).split("_")[1][1:]
        self._block = os.path.basename(vep_output).split("_")[2][1:]

        if self._chr == "X":
            self._chr = "23"
        elif self._chr == "Y":
            self._chr = "24"
        # load pvcf df
        block_df = pd.read_csv(pvcf_block, sep="\t", header=None)
        block_df.columns = ["index", "chr", "block", "block_start", "block_end"]

        self._coord = block_df[(block_df["chr"].astype(str) == self._chr) & (block_df["block"].astype(str) == self._block)][["block_start", "block_end"]].values.tolist()[0]

    def _load_vep_output(self):
        """
        Read through output file from VEP
        Only keep missense_variant + syn variants
        """
        output_file = []
        with open(self._vep_output, "r") as vepIn:
            for line in vepIn:
                # skip header line intronic variant
                if "intron_variant" in line or "## " in line:
                    continue
                if "#Uploaded_variation" in line:  # header line
                    header_line = line.strip().split()
                else:
                    if not "CANONICAL=YES" in line:
                        continue

                    line = line.strip().split()
                    info = line[-1]
                    info_dict = dict( (n,v) for n,v in (a.split('=') for a in info.split(";") ) )
                    #info_df = pd.DataFrame(info_dict)
                    line[-1] = info_dict
                    output_file.append(line)
        #variants_to_keep = ["missense_variant", "synonymous_variant", "frameshift_variant", "inframe_deletion", "stop_gained", "missense_variant,splice_region_variant", "splice_region_variant,synonymous_variant", "inframe_insertion", "start_lost"]
        #pd.set_option('display.max_columns', None) 
        output_df = pd.DataFrame(output_file)
        output_df.columns = header_line
        output_df = output_df[output_df["CDS_position"] != "-"].reset_index(drop=True)
        # select every column but the last
        self._vep_df = output_df.iloc[:,:-1]
        # build all v with weights
        info_df = output_df['Extra'].apply(pd.Series)
        merged_weights = self._vep_df.join(info_df)
        merged_varity = self._add_varity_w(merged_weights, VARITY)
        return merged_varity

    def _add_evmut_w(self, merged_weights, evmut_file):
        """
        add EVMutation weightd df
        """
        pass


    def _add_varity_w(self, merged_weights, varity_file):
        """
        add VARITY to merged_weights df
        """
        all_genes = merged_weights["Uniprot_acc"].unique().tolist()
        print(all_genes)
        # convert ENSG to varity identifier
        variants_in_range = []
        # ["X", "Y", "M"] in chrs
        all_unique_change = merged_weights["#Uploaded_variation"].unique().tolist()
        with open(varity_file, "r") as varity_raw:
            for line in varity_raw:
                if "chr" in line: 
                    header = line.strip().split()
                    continue
                # select row with gene name 
                l = line.strip().split()
                if l[0] != self._chr:
                    continue
                if int(self._coord[0]) <= int(l[1]) <= int(self._coord[1]):
                    # check if the variant is in the range of this block 
                    loc = f"chr{l[0]}_{l[1]}_{l[2]}_{l[3]}"
                    if loc in all_unique_change:
                        variants_in_range.append([loc]+l)
        variants_in_range_df = pd.DataFrame(variants_in_range)
        if variants_in_range_df.empty:
            return merged_weights

        print(variants_in_range_df)
        print(["variant_name"] + header)
        variants_in_range_df.columns = ["variant_name"] + header
        print(variants_in_range_df[variants_in_range_df["p_vid"] == "P30304"])

        # merge variants with input merged_weights
        merged_weights = pd.merge(merged_weights, variants_in_range_df[["variant_name", "p_vid", "aa_pos", "aa_ref", "aa_alt", "VARITY_R", "VARITY_ER"]], how="left", 
                left_on="#Uploaded_variation", right_on="variant_name")
        print(merged_weights[merged_weights["p_vid"] == "P30304"])
        return merged_weights

    def _merge_vcf(self):
        """
        After getting filtered vep df, merge the data with vcf file
        """
        # open vcf file for reading
        # parsed_vcf_file = os.path.join(self._output, f"{os.path.basename(self._vcf).replace('.vcf.gz', '_filtered.vcf')}")
        list_loc = set(self._vep_df["#Uploaded_variation"].tolist())
        all_mut = []
        with gzip.open(self._vcf, "rb") as in_vcf:
            for line in in_vcf:
                l = line.decode("utf-8").strip()
                if l.startswith("##"):  # skip header
                    continue
                if l.startswith("#CHROM"): # header line contains eid info
                    eid_load = np.array(l.split())[9:]
                    continue

                l = np.array(l.split())
                # skip variants that are not in the list
                variant_name = l[2].split(";")
                if list_loc.intersection(set(variant_name)):

                    indices = [i for i, s in enumerate(l[9:]) if ("0/0" not in s)]
                    gt_with_v = l[9:][indices]
                    eid_with_v = eid_load[indices]
                    # output variant_name, quality and AF/AQ
                    mut = [l[2], l[5], l[7]] +eid_with_v.tolist() + gt_with_v.tolist()
                    # process this entry
                    processed_mut = self.parse_mut(mut)
                    if not processed_mut.empty:
                        all_mut.append((processed_mut))
                    # l = "\t".join(l)
                   # out_vcf.write(l+"\n")
        all_mut = pd.concat(all_mut)
        all_mut[["chr", "pos", "ref", "alt"]] = all_mut["variant_name"].str.split("_", expand=True)
        return all_mut

    @staticmethod
    def parse_mut(entry):
        """
        """
        # get list of mutation in this entry (each entry might contain multiple mutations
        mut = entry[0].split(";")
        # get AF for each mut
        af = entry[2].split(";")[0].replace("AF=", "").split(",")
        # all the eids and genotype
        eid_gt = entry[3:]
        # divide eid and gt into two lists
        half = len(eid_gt)/2

        # error check: if we get .5 means there is something wrong
        if half != math.floor(half):
            raise ValueError("eid GT list can't be divided in half")
        eid = np.array(eid_gt[:int(half)])
        gt = np.array(eid_gt[int(half):])
        mut_index = 1
        all_mut = []
        mut_df = pd.DataFrame({})
        for i in mut:
            # if more than 10% of people have the gt ./., skip this variant
            missing_ind = [i for i,s in enumerate(gt) if (f"./." in s) or (f"{mut_index}/." in s) or (f"./{mut_index}" in s)]
            missing_perc = len(missing_ind)/total_n
            if missing_perc > 0.1:
                continue

            # select gt contains this alt
            # output eid list and gt list for people with this mut
            ind = [i for i, s in enumerate(gt) if (f"0/{mut_index}" in s) or (f"{mut_index}/0" in s)
                   or (f"{mut_index}/{mut_index}" in s)]

            eid_with_mut = eid[ind]
            gt_with_mut = gt[ind]
            mut_df = pd.DataFrame({"eid": eid_with_mut, "gt": gt_with_mut})
            mut_df["variant_name"] = i
            mut_df["AF"] = af[mut_index-1]

            # split gt column into multiple columns, same order as INFO section
            # GT:DP:AD:GQ:PL:RNC
            if mut_df.empty:
                continue
            mut_df[["GT", "DP", "AD", "GQ", "PL", "RNC"]] = mut_df["gt"].str.split(":", expand=True)
            mut_index += 1
            all_mut.append(mut_df)
        if all_mut != []:
            all_mut = pd.concat(all_mut)
        else:
            return pd.DataFrame({})
        return all_mut

def main(arguments):
    """
    Given directory contains output txt files from VEP, and directory contains raw vcf files,
    Merge all the filtered missense/syn coding variants into one file per chr/block
    """
    # input vcf and vep file
    vcf_file = arguments.vcf
    vep_file = arguments.vep

    if (not os.path.isfile(vep_file)) or (not os.path.isfile(vcf_file)):
        raise FileNotFoundError("Check input VEP/VCF file!!")

    # vep and vcf files have the same base name
    base = os.path.basename(vcf_file).split(".")[0]
    if base != os.path.basename(vep_file).split(".")[0]:
        raise ValueError(f"Base name doesn't match! {base}, {os.path.basename(vep_file).split('.')[0]}")

    parser = parseVEP(vep_file, vcf_file)
    # first load VEP to object
    merged_weights = parser._load_vep_output()
    merged_weights_file = f"{base}_all_weights.csv"
    merged_weights.to_csv(merged_weights_file)
    print(merged_weights_file)
    all_mut_for_block = parser._merge_vcf()

    # make output file
    mut_file = f"{base}_filtered_mut.csv"
    all_mut_for_block.to_csv(mut_file)
    print(mut_file)


if __name__ == '__main__':
    # required arguments
    parser = argparse.ArgumentParser(description="Parse files from vep analysis (per chr/block)")
    parser.add_argument("-vep", help="VEP txt file", required=True)
    parser.add_argument("-vcf", help="VCF file", required=True)
    # parser.add_argument("-cb", help="chromosome and block number, each time process one chromosome/block combination, eg. -cb 2 3 means chr2_block3", nargs=2, required=True)
    # parser.add_argument("-gInfo", help="gene info json file")
    # parser.add_argument("-o", help="Output folder", required=True)
    args = parser.parse_args()
    main(args)
