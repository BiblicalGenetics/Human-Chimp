# This program is designed to test the similarity between the human and chimp
# genomes using BLASTn. It contains three subroutines.
# main() sets up a series of BLAST searches, parameters can be set as desired.
# The length of the query strings are set with the {sample_length} variable.
# The {input_fasta_file} contains a set of sequences in FASTA format.
# main() will choose a random starting location on each sequence and create a
# subsequence that is {sample_length} nucleotides long.
# These are sent to blast_sequence(). This subroutine performs the searches and
# saves the results in {blast_results_directory}.
# The combine_blast_results() subroutine takes the reports in {blast_results_directory} 
# and combines all reports whose names include the test string {model} into a
# single file and saves this in the {save_directory}.
# If the reports contain more than one line, combine_blast_results() will only 
# retrieve the first hit.
# This is a working program, meaning it is subject to modification and change. I 
# make no promises as to readability or efficiency, but I can attest that it works
# as designed.
# Robert Carter, 28 Sep 2023.

import subprocess
import os
import io
from io import StringIO
import pandas as pd
import random
from Bio import SeqIO

# Pick one:
#input_fasta_file = "E:/PT3.1.1/Sequences/PT3.fna"                              # Chimp PanTro3.1.1 PT3
input_fasta_file = "E:/T2T/Sequence/chm13v2.0.fa"                               # Human Telomere-to-Telomere T2T

# Pick one:
blast_db = "E:\PT3.1.1\BlastDB\PT3full_db"                                      # PT3, \ is required as are all caps and lc on Windows machines
#blast_db = "E:\T2T\BLASTDB\chm13_db"                                           # T2T, \ is required as are all caps and lc on Windows machines

# Set search parameters
sample_length = 300                                                             # One sample of this length will be selected from each sequence in the input file
loops = 100                                                                     # Use to repeat the loop X times, each sequence in the input file will be processed this many times

# Change to suit, carefully
model = "T2TonPT3-300u"                                                         # Use this to save the results from different model runs separately

# Set the save directories
blast_results_directory = "e:/T2T/BlastResults/"                                # Where the blast reports will be saved
log_file = "E:/T2T/Blast_log_file.csv"                                          # Where the log file will be saved
save_directory = "E:/T2T/"                                                      # Where the summary file will be saved

print ("Blasting...")

# * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * 

def main():

    processed_sequence_count = 1

    # Setup BLAST parameters, most of these are default setings, additional parameters could be added if necessary
    eval = '0.1'          # e-value, how likely is the query expected to appear in the target at random?
    mts = '1'             # Maximum number of target sequences in report
    mhsps = '1'           # Maximum number of high-scoring sequence pairs among the target sequences in report
    perc = '50'           # Percent identity threshold for reporting
    word = '11'           # Word size
    threads = '8'         # Number of threads, machine specific
    gapo = '3'            # Penalty for opening a gap
    gape = '3'            # Penalty for extending a gap
    reward = '2'          # Reward for a matching base pair
    penalty = '-3'        # Penalty for a mismatching base pair
    xdru = '20'           # X-drop for ungapped searches
    xdrg = '30'           # X-drop for gapped searches
    xdrgf = '100'         # Final x-drop for gapped searches
    dust = 'no'           # Dust filtering yes or no
    soft = 'false'        # Soft masking true or false
    gapped = 'ungapped'   # Gapping parameter

    # Loop X times
    for loop in range(loops):
       print ("Input FASTA:", input_fasta_file)
       print (f"*** LOOP {loop + 1} ***")
       for record in SeqIO.parse(input_fasta_file, "fasta"):
           sequence = str(record.seq)
           id = record.id
           sequence_length = len(sequence)
           start_position = random.randint(0, sequence_length - sample_length)
           end_position = start_position + sample_length
           if end_position > sequence_length:
               end_position = sequence_length
           subsequence = sequence[start_position:end_position]
           subseq_name = f"{id}_{start_position}_{sample_length}"
           print(processed_sequence_count, id, start_position) # use with long subsequence searches
   
           output_file_control = os.path.join(blast_results_directory, f"{id}_{model}_{start_position}_{sample_length}_Control.csv")
           output_file_test = os.path.join(blast_results_directory, f"{id}_{model}_{start_position}_{sample_length}_Test.csv")
           if os.path.exists(output_file_control):
               print(f"{processed_sequence_count} {output_file_control} already exists. Skipping.")
               processed_sequence_count += 1
               continue
           else:
               processed_sequence_count += 1
   
           print("   Control...")
           # Control parameters go here:
           gapped = "ungapped"
           blast_sequence(subsequence, output_file_control, blast_db, eval, mts, mhsps, perc, word, threads, gapo, gape, reward, penalty, xdru, xdrg, xdrgf, dust, soft, gapped)

           print("   Test...")
           # Parameters to test go here:
           gapped = "gapped"
           blast_sequence(subsequence, output_file_test, blast_db, eval, mts, mhsps, perc, word, threads, gapo, gape, reward, penalty, xdru, xdrg, xdrgf, dust, soft, gapped)

    combine_blast_results()
    bell = '\a'
    print(bell)

# * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * 

def blast_sequence(query_sequence, output_file, blast_db, eval, mts, mhsps, perc, word, threads, gapo, gape, reward, penalty, xdru, xdrg, xdrgf, dust, soft, gapped):

    blastn_executable = "e:/BLAST/blast-2.14.1+/bin/blastn"
    command = [
        blastn_executable,
        "-db", blast_db,
        "-evalue", eval,
        "-outfmt", "10 qid qstart qend sseqid sstart send pident nident length mismatch gapopen gaps evalue bitscore",
        "-max_target_seqs", mts,
        "-max_hsps", mhsps,
        "-perc_identity", perc,
        "-word_size", word,
        "-num_threads", threads,
        "-gapopen", gapo,
        "-gapextend", gape,
        "-reward", reward,
        "-penalty", penalty,
        "-xdrop_ungap", xdru,
        "-xdrop_gap", xdrg,
        "-xdrop_gap_final", xdrgf,
        "-dust", dust,
        "-soft_masking", soft,
    ]
    if gapped == "ungapped":
        command.append("-ungapped")

#    print ("Blasting away with", command)

    try:
        with open(output_file, "w") as f_out:
            # Run BLAST and pipe the output to the output file
            p = subprocess.Popen(command, stdin=subprocess.PIPE, stdout=f_out, stderr=subprocess.PIPE, text=True)
            stdout, stderr = p.communicate(query_sequence)
            with open(log_file, "a") as log:
                log.write("Executed BLAST command:\n")
                log.write(" ".join(command) + "\n\n")
            if p.returncode != 0:
                print(f"Error occurred during BLAST search: {stderr}")
    except Exception as e:
        print(f"Error occurred: {str(e)}")

# * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * 

def combine_blast_results():
    print ("Combining output data...")
    results = {}
    data_dict = {}
 
    # Get chromosome names for the PT3 sequences
    mapping_file_path = "E:/PT3.1.1/Sequence names.csv"
    sequence_mapping = {}
    with open(mapping_file_path, "r") as mapping_file:
        lines = mapping_file.readlines()
        for line in lines[1:-1]:  # Skip header and footer
            columns = line.strip().split(",")
            new_name = columns[1]
            old_name = columns[0]
            sequence_mapping[old_name] = new_name
 
    # Iterate through the sequence files and extract data
    file_list = os.listdir(blast_results_directory)
    filtered_files = [filename for filename in file_list if model in filename]
    filenum = 0
    for filename in filtered_files:
        name_data = filename.split('_')
        seq_id = name_data[0]
        if seq_id not in data_dict:
            data_dict[seq_id] = {}
        start = name_data[2]
        if start not in data_dict[seq_id]:
            data_dict[seq_id][start] = {}
        length = name_data[3]
        if length not in data_dict[seq_id][start]:
            data_dict[seq_id][start][length] = {} 

        if name_data[4] == "Control.csv":
            print (seq_id, start, length, "control")
            file_path = os.path.join(blast_results_directory, filename)
            with open(file_path, 'r') as file:
                data = file.readline()
                control_data = data.rstrip("\n")
                data_dict[seq_id][start][length]['control'] = control_data

        if name_data[4] == "Test.csv":
            print (seq_id, start, length, "test")
            file_path = os.path.join(blast_results_directory, filename)
            with open(file_path, 'r') as file:
                data = file.readline()
                test_data = data.rstrip("\n")
                data_dict[seq_id][start][length]['test'] = test_data

    results_file_path = os.path.join(save_directory, f"{model} combined_blast_results.csv")
    with open(results_file_path, "w") as output_file:
        output_string = ",,,Control,,,,,,,,,,,,,,,Test,,,,,,,,,,,,,,\n"
        output_string += "QID,start,len,qstart,qlen,sseqid,sstart,send,pident,nident,length,mismatch,gapopen,gaps,evalue,bitscore,qstart,qlen,sseqid,sstart,send,pident,nident,length,mismatch,gapopen,gaps,evalue,bitscore\n"
        output_file.write(output_string)
        for seq_id, start_data in data_dict.items():
            for start, len_data in start_data.items():
                for length, data in len_data.items():
                    control_data = data.get("control", ",,,,,,,,,,,,,,")
                    test_data = data.get("test", ",,,,,,,,,,,,,,")
                    output_string = f"{seq_id},{start},{length},"
                    output_string += f"{control_data},"
                    output_string += f"{test_data}\n"
                    for old_name, new_name in sequence_mapping.items():
                        output_string = output_string.replace(old_name, new_name)
                    output_file.write(output_string)

# * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * 

if __name__ == "__main__":

    main()