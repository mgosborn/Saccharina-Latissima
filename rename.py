import os,re,sys
import gzip
import subprocess

# This script renames all files in a folder to preferred format.

## the files we want to rename are in the format:
## alpha_NanAmp_Sac1_species_sampleid_num_bases_Sac2_sequencer_lane_read.fastq
## examples:
### IZXU_NanoAmplified_Saccharina_latissima_SL-JS-19-FG-2_1_TCCAACGC_Saccharina_I997_L1_R1.fastq
### IZXU_NanoAmplified_Saccharina_latissima_SL-JS-19-FG-2_1_TCCAACGC_Saccharina_I997_L1_R2.fastq
### IZYN_NanoAmplified_Saccharina_angustissima_SA-CB-5-MG-3_1_GCAATGCA_Saccharina_I997_L1_R1.fastq
### IZYN_NanoAmplified_Saccharina_angustissima_SA-CB-5-MG-3_1_GCAATGCA_Saccharina_I997_L1_R2.fastq

## Some files (sequenced by machine I1018 or I1019 have additional _num_ before sampleid

### Renamed format:
## sampleid_sequencer_lane_read.fastq

print("Creating list of allreads fqs in folder...")
cmd = "ls *NanoAmplified*.fastq > files.txt"
os.system(cmd)

print("Renaming...")

with open('files.txt','r') as files:
    for file in files:
        file_stripped = file.strip()
        alts = ['I1018','I1019']
        if any(item in file_stripped for item in alts):
            sampleid = file_stripped.rsplit('_')[5]
            sequencer = file_stripped.rsplit('_')[9]
            lane = file_stripped.rsplit('_')[10]
            read_fq = file_stripped.rsplit('_')[11]
        else:
            sampleid = file_stripped.rsplit('_')[4]
            sequencer = file_stripped.rsplit('_')[8]
            lane = file_stripped.rsplit('_')[9]
            read_fq = file_stripped.rsplit('_')[10]
        newname = sampleid + '_' + sequencer + '_' + lane + '_' + read_fq
        #print(newname)
        cmd = 'mv ' + file_stripped + ' ' + newname
        #print(cmd)
        os.system(cmd)
    files.close()
