#Code written by Darren Mallett, Biggins Lab @ Fred Hutchinson Cancer Center. This program has only been tested on Mac.
version = '1 [Published Version]'

print(f'Running: Version {version}')
print('Initializing... please wait...')

import os, sys, filereader as fr, urllib.parse, mass_spec as ms, re, proteins, compare_files, time

def print_function(output):
    print("--------> " + output)
    
start_over = "YES"
list_of_successful_files = {} # a dictionary containing a list of files that were uploaded successfully to analyze.

#user-set commands:
minimum_PSM = 1 # Does the user want to only show phosphorylations that have greater than X PSM counts? anything ABOVE this (not including this value)
minimum_confidence = float(75) # % confidence requried to report phosphorylation. This will filter all phosphorylation calls that are less than X percent.

while start_over == "YES":

    print('Place Excel (.xlsx) file(s) exported from ProteomeDiscoverer that you want to analyze and compare in folder "msfiles". The program will analyze and compare all files in this folder. Make sure to remove old files you do not want analyzed. Note that the Excel files cannot be currently opened in Excel.\nContinue? Type Y')
    cont = input()
    cont = cont.upper()
    
    if cont == 'Y':

        filepath = os.path.dirname(os.path.realpath(__file__)) + '/msfiles'
        accepted_file_exts = ["xlsx", "xls", "xlsb", "xlsm"]
        
        for filename in os.listdir(filepath):
            
            file_ext = filename.split('.')
            file_ext = file_ext[-1]
            file_ext = file_ext.lower()
                    
            #delete all .csv files from the previous run, so the folder doesn't accumulate files.
            
            if file_ext == 'csv':
                
                os.remove(filepath + '/' + filename)
                
        #now that they have been deleted, loop through each file again and create csv files:
        
        count = 0 #file number....
        
        for filename in os.listdir(filepath):
            
            file_ext = filename.split('.')
            file_ext = file_ext[-1]
            file_ext = file_ext.lower()
            
            if file_ext in accepted_file_exts:

                file_uploads = [] #a list of successful file uploads? True = successful upload; False = not successful upload
                
                #for each file in the input folder, loop through the FileReader to convert each file to a CSV.

                new_filename = filename.strip()
                read_file = fr.FileReader(new_filename)

                if read_file.Error == None:
                    # we had a successful file upload:
                    file_uploads.append(True)
                    list_of_successful_files[count] = read_file.csv_filename
                    print_function("Success: The file " + new_filename + " was read successfully...")
                    count += 1
                else:
                    file_uploads.append(False)
                    print_function("Error: " + read_file.Error)

    else:
        
        exit()

    #if there are successful and not successful uploads, ask the user if they want to continue or start over
    if len(file_uploads) >= 2:
        #there were more than 2 files that were attempted to be analyzed...
        if not all(file_uploads):
            #not all of the uploads were True
            #there was an unsuccessful upload. 
            #Ask the user if they want to continue or start over
            
            print_function("Continue the analysis on the file(s) that were successfully found? Type YES or NO")
            continue_analysis = input()
            
            if continue_analysis != "YES":
                print_function("Program will now quit...")
                exit()
            elif continue_analysis == "NO":
                print_function("Do you want to start over? Type YES or NO")
                
                start_over = input()
                if start_over == "NO":
                    #quit the program:
                
                    print_function("Program will now quit...")
                    exit()
                
            else:
                #the user wants to continue analyzing the files that uploaded successfully:
                #MassSpec.analyze
                start_over = "NO"
                break #break the while() loop and continue with analysis..
                
        else:
            #all the files are true... continue the analysis:
            start_over = "NO"
            break #break the while() loop and continue with analysis..
                
    elif 0 in file_uploads: 
        #the user uploaded one file and it was not read properly:
        print_function("Do you want to start over? Type YES or NO")

        start_over = input()
        if start_over == "NO":
            #quit the program:
            print_function("Program will now quit...")
            exit()
            
    elif 1 in file_uploads:
        #the user successfully uploaded a file:
        start_over = "NO"
        break #break the while() loop and continue with analysis..
        
print(f'\nFilter out peptides with less than or equal to {minimum_PSM} PSM(s)? Type \'OK\' or type your desired PSM cutoff (anything above this number will be included in results):')
minimum_PSM_input = input().strip().upper()

if minimum_PSM_input != 'OK':
    minimum_PSM = int(minimum_PSM_input)

print(f'\nPhosphorylation call confidence minimum: {minimum_confidence} (%). Type \'OK\' or type your desired confidence minimum (any phos with this confidence or above will be included in results):')
min_phos_conf_input = input().strip().upper()
if min_phos_conf_input != 'OK':
    minimum_confidence = float(min_phos_conf_input)

print_function("Analyzing... please wait")

ms_files = {}
    
for key, file in list_of_successful_files.items():
    
    print(f"- - - - >>> Analyzing: {file}")
    #analyze the MS data...
    ms_files[file] = ms.MassSpec(file).analyze(minimum_PSM, minimum_confidence) #analyze MS files with the user-set settings
    
compare = compare_files.compare_all_phosphorylation(ms_files)

if compare:
    
    time.sleep(1)
    print_function("Complete! Your file is now available in the 'output' folder.")