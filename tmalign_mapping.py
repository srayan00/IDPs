import pandas as pd
import numpy as np
import os
import sys
import gzip
import shutil
import json
import re
import argparse
import requests
import Bio.PDB as bpdb
from Bio import Align

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

"""
POST MEETING STUFF
TODOS
1) When I compare two proteins, I want to first make sure the sequences are aligned correctly using NWAlign
2) Then I want to delete the unnecessary/misaligned parts of the structure
3) Then I want to compare the structures using TMScore 
"""

# Parser
parser = argparse.ArgumentParser()
parser.add_argument('--uniprotID', default='BIOD_MYCTU', type=str,
                    help='Uniprot ID for protein you want to compare')
parser.add_argument('--clustering', default= 1, type=int,
                    help='0 - all proteins with the same uniprot id, 1 - all proteins containing the same domains, 2 - all proteins in the same sequence cluster')
parser.add_argument("--testing_phase", default = 1, type = int,
                    help = "Testing phase or not")
parser.add_argument("--playing", default = 1, type = int,
                    help = "Playing around or not")
parser.add_argument("--similarity_score", default = 1 , type = int,
                    help = "0 - TM score, 1 - TM align")

# Aligner
aligner = Align.PairwiseAligner()
aligner.mode = 'local'
aligner.open_gap_score = -11
aligner.extend_gap_score = -1
aligner.substitution_matrix = Align.substitution_matrices.load("BLOSUM62")

# Functions 

def read_zipped_status(filepath = "/nfs/turbo/lsa-tewaria/zipped_status_efficient.csv", uniprot_ids = None):
    if not os.path.exists(filepath):
        zipped_status = pd.DataFrame({"status" : [True for i in range(len(uniprot_ids))]}, index = uniprot_ids)
        zipped_status.to_csv(filepath)
    else:
        zipped_status = pd.read_csv(filepath, index_col = 0)
    return zipped_status

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

def extract_sequence(json_file):
    return json_file["sequence"]["value"]

def extract_uniprot_id(json_file):
    return json_file["primaryAccession"]

def unzip_cif_folder(directory):
    for filename in os.listdir(directory):
        if filename.endswith(".gz"):
            filename_new = os.path.join(directory, filename)
            with gzip.open(filename_new, 'rb') as f_in:
                with open(filename_new[:-3], 'wb') as f_out:
                    shutil.copyfileobj(f_in, f_out)



def write_chain_list(directory, chain_list = None, filename = "chain_list"):
    with open(directory + "/" + filename, "w") as f:
        if chain_list is None:
            for filename in os.listdir(directory):
                if filename.endswith(".cif"):
                    f.write(filename + "\n")
        else:
            for chain in chain_list:
                f.write(chain + "\n")



def uniprot_to_pdb(uniprot_id, uniprot_df, write = True):
    uniprot_path = os.path.join("/nfs/turbo/lsa-tewaria/uniprot/", uniprot_id)
    mmcif_path = os.path.join("/nfs/turbo/lsa-tewaria/", "mmCIF")
    if not os.path.exists(uniprot_path):
        os.mkdir(uniprot_path)
    uniprot_pdb_ids = uniprot_df[uniprot_df["uniprot_id"] == uniprot_id]["from"]
    for pdb_id in uniprot_pdb_ids:
        pdb_dir_path = os.path.join(mmcif_path, pdb_id[1:3])
        pdb_path = os.path.join(pdb_dir_path, pdb_id + ".cif.gz")
        shutil.copy(pdb_path, uniprot_path)
    if write:
        write_chain_list(uniprot_path)
    unzip_cif_folder(uniprot_path)
    return uniprot_path  

def get_true_residue_boundaries(seq1, seq2):
    alignment = aligner.align(seq1, seq2)
    start_1, end_1 = alignment[0].aligned[0][0][0], alignment[0].aligned[0][0][1]
    start_2, end_2 = alignment[0].aligned[1][0][0], alignment[0].aligned[1][0][1]
    return start_1, end_1, start_2, end_2

def get_interpro_domain_boundaries(uniprot_id):
    response = requests.get("https://www.ebi.ac.uk/interpro/api/entry/all/protein/reviewed/" + uniprot_id + "/")
    if response.status_code == 200:
        json_file = response.json()["results"]
        domain_stuff = {x["metadata"]["accession"]: x for x in json_file if x["metadata"]["type"] == "domain"}
        domain_boundaries = pd.Series({k: [(x["fragments"][0]["start"],x["fragments"][0]["end"]) 
                                   for x in v["proteins"][0]["entry_protein_locations"]] for k, v in domain_stuff.items() 
                                   if v["proteins"][0]["entry_protein_locations"] is not None})
        domain_db = pd.Series({k: v["metadata"]["source_database"] for k, v in domain_stuff.items() if v["proteins"][0]["entry_protein_locations"] is not None}).reset_index().rename(columns = {"index" : "domain_id", 0 : "database"})
        domain_boundaries = domain_boundaries.explode().reset_index().rename(columns = {"index" : "domain_id", 0 : "domain_start_end"})
        domain_boundaries["domain_start"] = domain_boundaries["domain_start_end"].apply(lambda x: x[0])
        domain_boundaries["domain_end"] = domain_boundaries["domain_start_end"].apply(lambda x: x[1])
        domain_boundaries = domain_boundaries.join(domain_db.set_index("domain_id"), on = "domain_id").drop("domain_start_end", axis = 1)
        return domain_boundaries
    else:
        print("Could not find domain boundaries for {}".format(uniprot_id))
        return None

def remove_unwanted_residues(chain, start, end):
    ids_to_remove = [res.get_id() for i, res in enumerate(chain) if i < start or i > end]
    for res_id in ids_to_remove:
        chain.detach_child(res_id)
    return chain

def split_pdb_structure(uniprot_id, json_file):
    uniprot_path = os.path.join("/nfs/turbo/lsa-tewaria/uniprot/", uniprot_id)
    main_seq = extract_sequence(json_file)
    residue_boundaries = {"name" : [], "start" : [], "end" : []}
    for pdb_id in os.listdir(uniprot_path):
        if not pdb_id.endswith(".cif"):
            continue
        pdb_path = os.path.join(uniprot_path, pdb_id)
        pdb = bpdb.MMCIFParser().get_structure(pdb_id, pdb_path)
        for model in pdb:
            polypeptides = bpdb.PPBuilder().build_peptides(model)
            for chain, pp in zip(model, polypeptides):
                model_id, chain_id = str(model.get_id()), str(chain.get_id())
                name_of_chain = pdb_id[:-4] + "_"+ model_id + chain_id
                chain_path = os.path.join(uniprot_path, name_of_chain + ".pdb")
                start_pdb, end_pdb, start_main, end_main = get_true_residue_boundaries(pp.get_sequence(), main_seq)
                # if start_1 > 0:
                chain = remove_unwanted_residues(chain, start_pdb, end_pdb)
                io = bpdb.PDBIO()
                io.set_structure(chain)
                io.save(chain_path)
                residue_boundaries["name"].append(name_of_chain + ".pdb")
                residue_boundaries["start"].append(start_main)
                residue_boundaries["end"].append(end_main)
    return pd.DataFrame(residue_boundaries)

def TMalign_two_proteins(directory, protein_1, protein_2):
    terminal_command = "./TMalign " + directory + "/" + protein_1 + ".cif " + directory + "/" + protein_2 + ".cif"
    os.system(terminal_command)


def TMalign_dir(directory, chain_list = "chain_list", output_file = "TMalign_output.txt"):
    terminal_command = "./TMalign -dir " + directory + "/ " + directory + "/" + chain_list + " -outfmt 2 > " + directory + "/" + output_file
    os.system(terminal_command)

def TMscore_dir(directory, chain_list = "chain_list", output_file = "TMscore_output.txt"):
    terminal_command = "./TMscore -seq -dir " + directory + "/ " + directory + "/" + chain_list + " -outfmt 2 > " + directory + "/" + output_file
    os.system(terminal_command)

def clean_TMalign_output(directory, raw_file = "TMalign_raw_output.txt", clean_file = None):
    tmalign_output = pd.read_csv(os.path.join(directory, raw_file), sep = "\t", skipfooter = 1)
    tmalign_output["Chain_1"] = tmalign_output["#PDBchain1"].apply(lambda x: x.split("/")[-1].split(".")[0])
    tmalign_output["Chain_2"] = tmalign_output["PDBchain2"].apply(lambda x: x.split("/")[-1].split(".")[0])
    uniprot_id = directory.split("/")[-1]
    tmalign_output["Uniprot_ID"] = uniprot_id
    if clean_file is None:
        tmalign_output.to_csv(os.path.join(directory, uniprot_id + "_TMalign_output.csv"))
    else:
        tmalign_output.to_csv(os.path.join(directory, clean_file))

def clean_TMscore_output(directory, raw_file = "TMscore_raw_output.txt", clean_file = None):
    tmscore_output = pd.read_csv(os.path.join(directory, raw_file), sep = "\t", skipfooter = 1)

def get_domain_boundaries(json_file, scop_file):
    uniprot_id = extract_uniprot_id(json_file)
    uniprot_scop_file = scop_file[scop_file["FA-UNIID"] == uniprot_id][["FA-DOMID", "FA-UNIID", "FA-UNIREG"]]
    uniprot_scop_file["FA-UNIREG"] = uniprot_scop_file["FA-UNIREG"].str.split(",")
    uniprot_scop_file = uniprot_scop_file.explode("FA-UNIREG")
    uniprot_scop_file[["domain_start", "domain_end"]] = uniprot_scop_file["FA-UNIREG"].str.split("-", expand = True).astype(int)
    uniprot_scop_file = uniprot_scop_file.drop(["FA-UNIREG", "FA-UNIID"], axis = 1).rename(columns = {"FA-DOMID" : "domain_id"})
    uniprot_scop_file["database"] = "SCOP"
    other_domain_boundaries = get_interpro_domain_boundaries(uniprot_id)
    return pd.concat([uniprot_scop_file, other_domain_boundaries]).reset_index(drop = True)


def get_domain_info_for_pdbs(json_file, scop_file, uniprot_id):
    domain_boundaries = get_domain_boundaries(json_file, scop_file)
    residue_boundaries = split_pdb_structure(uniprot_id, json_file)
    for x in domain_boundaries["domain_id"].unique():
        dom_start = domain_boundaries[domain_boundaries["domain_id"] == x]["domain_start"].values
        dom_end = domain_boundaries[domain_boundaries["domain_id"] == x]["domain_end"].values
        residue_boundaries["d_{}".format(x)] = ((residue_boundaries["start"] <= min(dom_start)) & (residue_boundaries["end"] >= max(dom_end)))
    return residue_boundaries

def get_sequence_clusters(directory, domain_info, uniprot_json):
    fasta_path = os.path.join(directory, "fasta_file.fasta")
    main_sequence = extract_sequence(uniprot_json)
    with open(fasta_path, "w") as f:
        for x in range(len(domain_info["name"])):
            start = domain_info.iloc[x]["start"]
            end = domain_info.iloc[x]["end"]
            f.write(f">{domain_info.iloc[x]['name']}\n")
            f.write(main_sequence[start:end] + "\n")
    terminal_command = "mmseqs easy-cluster " + fasta_path + " " + os.path.join(directory, "clusterRes") + " tmp --min-seq-id 0.8 -c 0.8 --cov-mode 1"
    os.system(terminal_command)
    cluster_put = pd.read_csv(os.path.join(directory, "clusterRes_cluster.tsv"), sep="\t", header=None, names=["cluster", "pdb"])
    return cluster_put

def compare_proteins_domain(directory, domain_info):
    list_of_domains = list(domain_info.columns[3:])
    unique_combos = domain_info.groupby(list_of_domains)
    i = 0
    for _, group in unique_combos:
        if len(group) > 1:
            print(list(group["name"]))
            write_chain_list(directory, chain_list = list(group["name"]), filename = "chain_list_{}".format(i))
            TMalign_dir(directory, chain_list = "chain_list_{}".format(i), output_file = "TMalign_raw_output_{}.txt".format(i))
            clean_TMalign_output(directory, raw_file = "TMalign_raw_output_{}.txt".format(i), clean_file = "TMalign_output_{}.csv".format(i))
            i += 1

def compare_proteins_seq_clusters(directory, domain_info, uniprot_json):
    seq_clusters = get_sequence_clusters(directory, domain_info, uniprot_json)
    unique_combos = seq_clusters.groupby("cluster")
    i = 0
    for _, group in unique_combos:
        if len(group) > 1:
            write_chain_list(directory, chain_list = list(group["pdb"]), filename = "chain_listseq_{}".format(i))
            TMalign_dir(directory, chain_list = "chain_listseq_{}".format(i), output_file = "TMalign_raw_outputseq_{}.txt".format(i))
            clean_TMalign_output(directory, raw_file = "TMalign_raw_outputseq_{}.txt".format(i), clean_file = "TMalign_outputseq_{}.csv".format(i))
            i += 1




if __name__ == "__main__":

    args, unknown = parser.parse_known_args()
    uniprot_id = args.uniprotID
    clustering = args.clustering
    testing_phase = args.testing_phase
    playing = args.playing
    similarity_score = args.similarity_score

    # Load dataframe PDBs to Uniprot
    if testing_phase:
        uniprot_pdb_path = "/nfs/turbo/lsa-tewaria/uniprot_df_small.csv"
    else:
        uniprot_pdb_path = "/nfs/turbo/lsa-tewaria/uniprot_df_2.csv"
    uniprot_df = pd.read_csv(uniprot_pdb_path)

    # Load Scop data
    scop_df = read_scop_file()

    if playing:
        uniprot_json = find_uniprot(uniprot_df, uniprot_id)
        uniprot_path = uniprot_to_pdb(uniprot_id, uniprot_df)
        domain_info = get_domain_info_for_pdbs(uniprot_json, scop_df, uniprot_id)
        domain_info.to_csv(os.path.join(uniprot_path,"domain_info.csv"))
        if clustering == 1:
            compare_proteins_domain(uniprot_path, domain_info)
        elif clustering == 0:
            TMalign_dir(uniprot_path)
        else:
            compare_proteins_seq_clusters(uniprot_path, domain_info, uniprot_json)
        exit()

    if uniprot_id == "all":
        # Get Unique Uniprot IDs
        uniprot_ids = uniprot_df["uniprot_id"].unique()
    
        for uniprot_id in uniprot_ids:
            # Find Uniprot data
            uniprot_json = find_uniprot(uniprot_df, uniprot_id)

            # Move zipped cif files from pdb to uniprot folder
            uniprot_path = uniprot_to_pdb(uniprot_id, uniprot_df)

             # Get the domain info for the proteins
            domain_info = get_domain_info_for_pdbs(uniprot_json, scop_df)

            # Compare proteins in Uniprot folder
            if clustering == 1:
                compare_proteins_domain(uniprot_path, domain_info)
            elif clustering == 0:
                TMalign_dir(uniprot_path)
            else:
                compare_proteins_seq_clusters(uniprot_path, domain_info, uniprot_json)
    else:
        uniprot_json = find_uniprot(uniprot_df, uniprot_id)
        uniprot_path = uniprot_to_pdb(uniprot_id, uniprot_df)
        domain_info = get_domain_info_for_pdbs(uniprot_json, scop_df)
        if clustering == 1:
            compare_proteins_domain(uniprot_path, domain_info)
        elif clustering == 0:
            TMalign_dir(uniprot_path)
        else:
            compare_proteins_seq_clusters(uniprot_path, domain_info, uniprot_json)

