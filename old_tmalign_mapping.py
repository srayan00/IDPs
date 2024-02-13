import pandas as pd
import numpy as np
import os
import sys
import gzip
import shutil
import json
import re
from tmalign_mapping import *

def read_zipped_status(filepath = "/nfs/turbo/lsa-tewaria/zipped_status.csv"):
    if not os.path.exists(filepath):
        list_mmcifs = os.listdir("/nfs/turbo/lsa-tewaria/mmCIF/")
        zipped_status = pd.DataFrame({"status" : [True for i in range(len(list_mmcifs))]}, index = list_mmcifs)
        zipped_status.to_csv(filepath)
    else:
        zipped_status = pd.read_csv(filepath, index_col = 0)
    return zipped_status

def uniprot_to_pdb(uniprot_id, uniprot_df, zip_status, write = True):
    uniprot_path = os.path.join("/nfs/turbo/lsa-tewaria/uniprot/", uniprot_id)
    mmcif_path = os.path.join("/nfs/turbo/lsa-tewaria/", "mmCIF")
    if not os.path.exists(uniprot_path):
        os.mkdir(uniprot_path)
    uniprot_pdb_ids = uniprot_df[uniprot_df["uniprot_id"] == uniprot_id]["from"]
    for pdb_id in uniprot_pdb_ids:
        pdb_dir_path = os.path.join(mmcif_path, pdb_id[1:3])
        if zip_status["status"][pdb_id[1:3]]:
            unzip_cif_folder(pdb_dir_path)
            zip_status["status"][pdb_id[1:3]] = False
        pdb_path = os.path.join(pdb_dir_path, pdb_id + ".cif")
        shutil.copy(pdb_path, uniprot_path)
    if write:
        write_chain_list(uniprot_path)
    return uniprot_path, zip_status

def compare_proteins_domain(directory, domain_info):
    list_of_domains = list(domain_info.columns[3:])
    unique_combos = domain_info.groupby(list_of_domains)
    i = 0
    for _, group in unique_combos:
        if len(group) > 1:
            write_chain_list(directory, chain_list = [x.lower() + ".cif" for x in list(group.index)], filename = "chain_list_{}".format(i))
            compare_proteins_dir(directory, chain_list = "chain_list_{}".format(i), output_file = "TMalign_raw_output_{}.txt".format(i))
            clean_TMalign_output(directory, raw_file = "TMalign_raw_output_{}.txt".format(i), clean_file = "TMalign_output_{}.csv".format(i))
            i += 1

def get_domain_info_for_pdbs(json_file, scop_file):
    domain_boundaries = get_domain_boundaries(json_file, scop_file)
    residue_boundaries = get_residue_boundaries(json_file)
    for x in domain_boundaries["domain_id"].unique():
        dom_start = domain_boundaries[domain_boundaries["domain_id"] == x]["domain_start"].values
        dom_end = domain_boundaries[domain_boundaries["domain_id"] == x]["domain_end"].values
        residue_boundaries["d_{}".format(x)] = ((residue_boundaries["start"] <= min(dom_start)) & (residue_boundaries["end"] >= max(dom_end)))
    return residue_boundaries

def get_residue_boundaries(json_file):
    sample = {x["id"] : {y["key"] : y["value"] for y in x["properties"]} for x in json_file["uniProtKBCrossReferences"] if x["database"] == "PDB"}
    domain_boundaries = {idx : item["Chains"].split("=")[-1] for idx, item in sample.items() if "Chains" in item.keys()}
    domain_boundaries = pd.DataFrame.from_dict(domain_boundaries, orient = "index", columns = ["Domain_Boundaries"])
    domain_boundaries[["start", "end"]] = domain_boundaries["Domain_Boundaries"].str.split("-", expand = True).astype(int)
    domain_boundaries = domain_boundaries.drop("Domain_Boundaries", axis = 1).reset_index().rename(columns = {"index" : "name"})
    domain_boundaries["name"] = domain_boundaries["name"].str.lower() + ".cif"
    return domain_boundaries