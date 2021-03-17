import os,re,sys
import gzip
import subprocess

## This script takes a list of sample IDs, searches two folders for fastq files with the matching ID,
## and moves all matching files to a third directory for further analysis.

## format of sample IDs in list: SA-CB-9-FG-1

## Parent directory: /project/noujdine_61/mgosborn/Latissima/

## Pulling fastq files from: /project/noujdine_61/mgosborn/Latissima/Plate1_and_2 and /project/noujdine_61/mgosborn/Latissima/Plate3
## Destination directory: /project/noujdine_61/mgosborn/Latissima/Biomass_Subset

subsetList = 'list_of_60_biomass_sampleids.txt'

subsetFolder = '/project/noujdine_61/mgosborn/Latissima/Biomass_Subset/'
parentFolder = '/project/noujdine_61/mgosborn/Latissima/'
plates = ['Plate1_and_2','Plate3']

with open(subsetList,'r') as files:
    for plate in plates:
        files.seek(0)
        searchFolder = parentFolder + plate
        print(searchFolder)
        for file in files:
            file_stripped = file.strip() #get rid of any trailing/leading white space
            cmd = 'cp ' + searchFolder + '/*' + file_stripped + '* ' + subsetFolder
            print(cmd)
            os.system(cmd)

cmd = 'sed -n -e \'1,20p\' ' + subsetList + ' > batch1.txt'
print(cmd)
os.system(cmd)
cmd = 'sed -n -e \'21,40p\' ' + subsetList + ' > batch2.txt'
print(cmd)
os.system(cmd)
cmd = 'sed -n -e \'41,60p\' ' + subsetList + ' > batch3.txt'
print(cmd)
os.system(cmd)

batches = ['batch1','batch2','batch3']
for batch in batches:
    cmd = 'for file in $(cat ' + subsetFolder + batch +'.txt); do mv "$file" /' + batch + '; done'
    print(cmd)
    os.system(cmd)
