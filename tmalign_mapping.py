import pandas as pd
import numpy as np
import os
import sys
import gzip
import shutil
import json
import re

"""
Okay! 
1) So, we want to loop through uniprot_df["uniprot_id"] and find the corresponding pdb files for that protein. (DONE)
2) We want to create a folder for each protein, and then put all the pdb files for that protein in that folder. (DONE)
3) Unzip all the mmcif files in the folder (DONE)
4) Use the following command to run the TM align program: ./TMalign -dir chain_folder/ chain_list (DONE)
5) Clean the TMalign output (DONE)

6) We want to make sure the protein structures are aligned correctly. Two proteins can only be compared if they have the same
residue boundaries. We can check the domain stuff here too using Scop?? (trying this)
-> Okay here is the plan for this part:
    -> We want to find the domain boundaries for each protein. We can do this using the uniprot_df "to" column. the uniProtKBCrossReferences json part has it all
    -> remove partial domain crap and proteins with no domain boundary info
    -> And then do tmalign

ENLYS_BPT4 has been the example protein so far. 
"""

def read_scop_file(filepath = "/nfs/turbo/lsa-tewaria/scop-cla-latest.txt"):
    scop_df = pd.read_csv(filepath, sep = " ", skiprows = 5)
    new_col_names = list(scop_df.columns)[1:]
    new_col_names.append("def")
    scop_df.columns = new_col_names
    scop_df.drop("def", axis = 1, inplace = True)
    return scop_df

def find_uniprot(df, uniprot_id):
    test = df[df["uniprot_id"] == uniprot_id]["to"]
    s = test.values[0]

    split_s = s.split("\"")
    split_s[0::2] = [x.replace("\'", '\"') for x in split_s[0::2]]
    # split_s[1::2] = [x.replace("\'", '') for x in split_s[1::2]]
    s = "\"".join(split_s)
    s = s.replace(": False", ": \"False\"")
    s = s.replace(": True", ": \"True\"")

    y = json.loads(r"{}".format(s))
    return y

def unzip_cif_folder(directory):
    for filename in os.listdir(directory):
        if filename.endswith(".gz"):
            filename_new = os.path.join(directory, filename)
            with gzip.open(filename_new, 'rb') as f_in:
                with open(filename_new[:-3], 'wb') as f_out:
                    shutil.copyfileobj(f_in, f_out)


def write_chain_list(directory):
    with open(directory + "/chain_list", "w") as f:
        for filename in os.listdir(directory):
            if filename.endswith(".cif"):
                f.write(filename + "\n")

def uniprot_to_pdb(uniprot_id, uniprot_df, zip_status):
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
    write_chain_list(uniprot_path)
    return uniprot_path, zip_status


def compare_two_proteins(directory, protein_1, protein_2):
    terminal_command = "./TMalign " + directory + "/" + protein_1 + ".cif " + directory + "/" + protein_2 + ".cif"
    os.system(terminal_command)


def compare_proteins_dir(directory, chain_list = "chain_list"):
    terminal_command = "./TMalign -dir " + directory + "/ " + directory + "/" + chain_list + " -outfmt 2 > " + directory + "/TMalign_raw_output.txt"
    os.system(terminal_command)

def clean_TMalign_output(directory):
    tmalign_output = pd.read_csv(os.path.join(directory, "TMalign_raw_output.txt"), sep = "\t", skipfooter = 1)
    tmalign_output["Chain_1"] = tmalign_output["#PDBchain1"].apply(lambda x: x.split("/")[-1].split(".")[0])
    tmalign_output["Chain_2"] = tmalign_output["PDBchain2"].apply(lambda x: x.split("/")[-1].split(".")[0])
    uniprot_id = directory.split("/")[-1]
    tmalign_output["Uniprot_ID"] = uniprot_id
    tmalign_output.to_csv(os.path.join(directory, uniprot_id + "_TMalign_output.csv"))

def get_domain_boundaries(json_file, scop_file):
    uniprot_id = json_file["primaryAccession"]
    uniprot_scop_file = scop_file[scop_file["FA-UNIID"] == uniprot_id]
    return uniprot_scop_file["FA-UNIREG"]

def get_residue_boundaries(json_file, scop_file):
    pass
    


if __name__ == "__main__":
    # Zip status
    zipped_status_path = "/nfs/turbo/lsa-tewaria/zipped_status.csv"
    if not os.path.exists(zipped_status_path):
        list_mmcifs = os.listdir("/nfs/turbo/lsa-tewaria/mmCIF/")
        zipped_status = pd.DataFrame({"status" : [True for i in range(len(list_mmcifs))]}, index = list_mmcifs)
        zipped_status.to_csv(zipped_status_path)
    else:
        zipped_status = pd.read_csv(zipped_status_path, index_col = 0)
    
    # Load dataframe PDBs to Uniprot
    uniprot_pdb_path = "/nfs/turbo/lsa-tewaria/uniprot_df_small.csv"
    uniprot_df = pd.read_csv(uniprot_pdb_path)

    # Load Scop data
    scop_df = read_scop_file()
    print(scop_df.head())

    # Get Unique Uniprot IDs
    uniprot_ids = uniprot_df["uniprot_id"].unique()


    
    # for uniprot_id in uniprot_ids:
    #     # Find Uniprot data
    #     uniprot_json = find_uniprot(uniprot_df, uniprot_id)

    #     # Move PDBs to Uniprot folders
    #     uniprot_path, zipped_status = uniprot_to_pdb(uniprot_id, uniprot_df, zipped_status)

    #     # Compare proteins in Uniprot folder
    #     compare_proteins_dir(uniprot_path)

    #     # Clean TMalign output
    #     clean_TMalign_output(uniprot_path)
    

    # Find Uniprot data
    uniprot_json = find_uniprot(uniprot_df, "CHEY_ECOLI")
    print(get_domain_boundaries(uniprot_json, scop_df))
    sample = [x for x in uniprot_json["uniProtKBCrossReferences"] if x["database"] == "PDB"]
    print(sample[:][["id", "properties"]])

    # # Move PDBs to Uniprot folders
    # uniprot_path, zipped_status = uniprot_to_pdb("CHEY_ECOLI", uniprot_df, zipped_status)

    # # Compare proteins in Uniprot folder
    # compare_proteins_dir(uniprot_path)

    # # Clean TMalign output
    # clean_TMalign_output(uniprot_path)

    # # Save zipped status
    # zipped_status.to_csv(zipped_status_path)
    
