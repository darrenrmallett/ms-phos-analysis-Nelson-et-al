import sys, os, proteins, re

current_filepath = os.path.dirname(os.path.realpath(__file__)) + '/msfiles'
current_output_path = os.path.dirname(os.path.realpath(__file__)) + '/output'

class MassSpec:
    
    #this class parses the MS data that has been already converted to a tab-separated file with a .csv extension.
    
    def __init__(self, filename):
    
        self.filename = filename
        
    def mods(self, mod_column, protein, peptide, mods_in_master, minimum_confidence):
        #reads the Modification column on a file line containing a peptide and determines whether there are mods.
        #returns values to put in 'mods' key of the peptide dictionary.
        
        mod_column = mod_column
        mod_string = None #the list item. default to no modifications.
        protein = protein #the protein name. Only needed if mods_in_master is False.
        peptide = peptide #the peptide. Only needed if mods_in_master is False.
        mods_in_master = mods_in_master #this is a boolean True or False. True = the modification info refers to residues in the master protein. False = refers to residues in the PEPTIDE. If this is false, we have to match the residue in the peptide to the residue in the master protein.
        minimum_confidence = float(minimum_confidence) #the mimumum confidence (percentage) cutoff of the modification. Set by the user in phospho.py. Default = 75 (%)
        
        # # # PHOSPHORYLATION # # # 
        
        #match to see if there is a phosphorylation mod:
        how_many_phosphos = re.search(r'(\d)(?:x)(?:Phospho)', mod_column)
        
        if how_many_phosphos is not None:
            
            #there is a phosphorylated residue in this peptide
            
            mod_string = {'phos': {}} #what will go into the mods list at the very end. each item in the list will be a dictionary.
            
            #get the residues:
            
            phospho_residues_search = re.search(r'(?:\d)(?:x)(?:Phospho \[)([A-Z0-9\(\).; \/\-]+)(?:\])', mod_column, re.I) #case insensitive
            phospho_residues = phospho_residues_search.groups()[0]
            phospho_residues = phospho_residues.split('; ') #they are separated by a ; and space if there are more than one
            
            count = 0
            
            for residue in phospho_residues:
                
                #get the residue, and it's confidence level (if it knows which residue).
                
                #some residues are NOT known. Therefore the MS software outputs that as "S/T" or "S/T/Y" without a residue number and without a confidence.
                #first, figure out if the residue is known or not.
                
                residue_has_position = re.search(r'(^[A-Z])(?:\d+)', residue) #this matches S394 for instance (a residue with known position)
                
                if residue_has_position is not None:
                    
                    #this residue position is known.
                    amino_acid = re.search(r'(^[A-Z])', residue)
                    amino_acid = amino_acid.group()[0]
                    
                    position = re.search(r'(?:[A-Z])([\d]+)', residue)
                    position = int(position.groups()[0])
                    
                    confidence = re.search(r'(?:\()([0-9.]+)(?:\))', residue)

                    if confidence is not None:
                        confidence = float(confidence.groups()[0])
                    else:
                        confidence = None
                    
                else:
                    
                    #this must be a residue or residue group that does not have a known position.
                    amino_acid = residue #whatever the residue value is...
                    
                    confidence = None
                    
                    position = None
                    
                #IMPORTANT: if mods_in_master is False, we have to match the residue to the residue in the master protein.
                
                if mods_in_master is False:
                    
                    peptide_positions = self.get_peptide_position(protein, peptide)
                    starting_residue = peptide_positions[0]
                    ending_residue = peptide_positions[1]  
                    
                    if residue_has_position:
                        
                        #update the position to the position in the MASTER PROTEIN (not in the peptide).
                        position = starting_residue + (position-1)
                
                if minimum_confidence is not float(0): #if there is no user-set minimum confidence, put the phosphorylation mark in the mod_string
                    
                    if confidence is not None: #do not include phosphorylation of unknown position
                        if confidence >= minimum_confidence:
                        
                            mod_string['phos'][count] = {

                                'AA': amino_acid,
                                'conf': confidence,
                                'pos': position

                            }     

                else:
                    
                    #there is no minimum confidence, so include everything (low confidence as well as no known position without any confidence)                
                    mod_string['phos'][count] = {

                        'AA': amino_acid,
                        'conf': confidence,
                        'pos': position

                    }

                    count += 1
                
        #check to see whether there's any acetylation detected!

        how_many_acetyls = re.search(r'(\d)(?:x)(?:Acetyl)', mod_column)

        if how_many_acetyls is not None:

            #there's an acetylation.

            #if there is not any phosphorylation in the mod string already...
            if mod_string is None:

                mod_string = {'acetyl': {}}

            else:

                mod_string['acetyl'] = {}

            #get the residues, if known:

            acetyl_residues_search = re.search(r'(?:\d)(?:x)(?:Acetyl \[)([A-Z0-9\(\).; \/\-]+)(?:\])', mod_column, re.I)  #case insensitive      
            acetyl_residues = acetyl_residues_search.groups()[0]
            acetyl_residues = acetyl_residues.split('; ') #they are separated by a ; and space if there are more than one

            count = 0

            for residue in acetyl_residues:

                #get the residue, and it's confidence level (if it knows which residue).

                #some residues are NOT known. Therefore the MS software outputs that as "S/T" or "S/T/Y" without a residue number
                #if the N-term is acetylated, it will output 'N-term'
                #first, figure out if the residue is known or not.

                residue_has_position = re.search(r'(^[A-Z])(?:\d+)', residue) #this matches S394 for instance (a residue with known position)

                if residue_has_position is not None:

                    #this residue position is known.
                    amino_acid = re.search(r'(^[A-Z])', residue)
                    amino_acid = amino_acid.group()[0]

                    confidence = re.search(r'(?:\()([0-9.]+)(?:\))', residue)
                    confidence = float(confidence.groups()[0])

                    position = re.search(r'(?:[A-Z])([\d]+)', residue)
                    position = int(position.groups()[0])

                else:

                    #this must be a residue or residue group that does not have a known position.
                    amino_acid = residue #whatever the residue value is...

                    confidence = None

                    position = None


                mod_string['acetyl'][count] = {

                    'AA': amino_acid,
                    'conf': confidence,
                    'pos': position

                }

                count += 1
                
        #delete the 'phos' key if there is nothing in it now.
        if type(mod_string) is dict:
            if 'phos' in mod_string.keys():
                if len(mod_string['phos']) == 0:
                    del(mod_string['phos'])

        return mod_string              
    
    def analyze(self, minimum_PSM = 0, minimum_confidence = 75):
        
        minimum_PSM = minimum_PSM
        self.minimum_PSM = minimum_PSM
        minimum_confidence = minimum_confidence
        self.minimum_conf = minimum_confidence
        
        self.proteins = {}
        
        #we want to open the mass spec data (CSV file) that we previously created.
        self.filehandle = open(current_filepath + '/' + self.filename, 'r')
        lines = self.filehandle.readlines() #read line by line...
        
        count = 1
        
        is_row_protein_info = False #to keep track of which row we are on...
        
        #pre-set some variables.
        accession_row = None
        description_row = None
        peptide_row = None
        modifications_row = None
        peptide_PSMs_row = None
        
        last_protein = None
        contaminant_protein = False #set to true if the rows refer to a contaminant protein.
        
        for line in lines:
            
            line = line.strip()
            
            #find which line contains the accession/identifier and description
            line_columns = line.split('\t')
            
            if count == 1:
                #this is the first row of the document, with the header info.
            
                master_row          = line_columns.index('Master')
                accession_row       = line_columns.index('Accession')
                description_row     = line_columns.index('Description')
                #coverage_row        = line_columns.index('Coverage [%]')
                num_peptides_row    = line_columns.index('# Peptides')
                #isoforms_row        = line_columns.index('# Isoforms')
                PSMs_row            = line_columns.index('# PSMs')
            
            else:
                #we are past the first row.
                
                #are we on a row that contains protein info? The "Master" column should have "Master Protein" if it is referring to the protein info.
                
                if line_columns[master_row] == 'Master Protein':
                    
                    #this is a row containing protein info...
                
                    #The 'Accession' row will contain the Y(accession) number.
                    match = re.search(r'(^[Y]\w+[-]*\w+)', line_columns[accession_row]) #this pattern will not select keratin contaminants and non-SGD stuff.

                    if match is not None:

                        #this is a row with protein info in it. It is not a contaminant protein
                        contaminant_protein = False
                        is_row_protein_info = True

                        #get the protein/gene name...

                        gene_name_match = re.search(r'^(?:["\'])?([A-Z0-9\-]+){1}', line_columns[description_row])
                        

                        if gene_name_match is not None:
                            
                            #check to see if the protein/gene name is in the Master List (master_list).

                            if gene_name_match.groups()[0] not in proteins.master_list.keys():

                                #this gene might have a different name in the MS file than in SGD (alias). Fetch the common/display name.

                                for common_name, info in proteins.master_list.items():

                                    if info['identifier'] == match.groups()[0]:

                                        protein_name = common_name

                            else:

                                protein_name = gene_name_match.groups()[0]

                            last_protein = protein_name

                            #enter the protein into the master protein dictionary for this MS file.
                            #get other info too:
                            PSM_count       = line_columns[PSMs_row]
                            #coverage        = line_columns[coverage_row]
                            num_peptides    = line_columns[num_peptides_row]
                            #isoforms        = line_columns[isoforms_row]

                            self.proteins[protein_name] = { 

                                #'coverage': coverage,
                                'num_peptides': num_peptides,
                                #'isoforms': isoforms,
                                'PSMs': PSM_count,
                                'peptides': {} 

                            }
                            

                    else:
                        
                        #this is a contaminant protein
                        contaminant_protein = True
                                
                else:
                    
                    #this is either a row with peptide info or the header for peptide info. Determine which one...

                    is_row_protein_info = False

                    if 'Annotated Sequence' in line or 'Modifications' in line:
                        #this must be the header info. get the rows... oops I mean columns but I already wrote the script as rows so keep that in mind.
                        
                        #IMPORTANT: Some files do not have 'Annotated Sequence' but instead have 'Sequence':
                        
                        if 'Annotated Sequence' in line:
                            peptide_row = line_columns.index('Annotated Sequence')
                            
                        elif 'Sequence' in line:
                            peptide_row = line_columns.index('Sequence')
                            
                        else:
                            #there is no peptide information. Throw an error.
                            print(f"ERROR: the following file does not contain peptide information in the Peptides & Proteins tab: {self.filename} \n\nPlease fix and re-run the program. The program will now quit.")
                            sys.exit()
                            
                        
                        #IMPORTANT: Some files do not have 'Modifications in Master Proteins' but rather 'Modifications' which makes my life SO much harder.
                        # 'Modifications' contains modifications with positions on the peptide, not the protein. Therefore if the file only contains 'Modifications' row then we have to map the modified residue back to the matching residue in the whole protein.
                        
                        mods_in_master = None
                            
                        if 'Modifications in Master Proteins' in line:
                            #default and hopefully the case:
                            modifications_row = line_columns.index('Modifications in Master Proteins')
                            mods_in_master = True
                            
                        elif 'Modifications' in line:
                            modifications_row = line_columns.index('Modifications')
                            mods_in_master = False
                        
                        #IMPORTANT: some files do not have 'Positions in Master Proteins', we need to do a workaround.
                        
                        if 'Positions in Master Proteins' in line:
                            positions_in_master_proteins_row = line_columns.index('Positions in Master Proteins')
                            
                        else:
                            positions_in_master_proteins_row = False
                            
                        peptide_PSMs_row = line_columns.index('# PSMs')
                        
                        #get peptide groups:
                        if '# Protein Groups' in line:
                            protein_groups_row = line_columns.index('# Protein Groups')
                        else:
                            protein_groups_row = None
                            
                        #get master protein accessions:
                        if 'Master Protein Accessions' in line:
                            master_accessions_row = line_columns.index('Master Protein Accessions')
                        else:
                            master_accessions_row = None

                    elif line != '':
                        #this line has something on it.
                        #this is a row containing peptide information
                        
                        #make sure this protein is not a contaminant protein. We do not want to analyze those peptides..
                        
                        if contaminant_protein is False:

                            #get the peptide (do not include upstream and downstream AA information that is given in the MS data...
                            peptide = re.search(r'^(?:\[[A-Z\-]{0,}\]\.)?([A-Z]+)(?:\.\[[A-Z\-]{0,}\])?', line_columns[peptide_row]) #searches any peptide 'word' with 2 or more characters, i.e. leaves out upstream and downstream info

                            if peptide is None:
                                print(f"ERROR: the following file does not contain peptide information in the Peptides & Proteins tab: {self.filename} \n\nPlease fix and re-run the program. The program will now quit.")
                                sys.exit()
                            
                            PSM_value = int(line_columns[peptide_PSMs_row])
                            
                            #check if the peptide is not already in the dictionary for this protein
                            if peptide is not None and peptide.groups()[0] in self.proteins[last_protein]['peptides'].keys():

                                #this peptide already exists, add to the PSMs (only IF it passes our minimum PSM cutoff):
                                                                
                                if PSM_value > minimum_PSM:
                                    
                                    self.proteins[last_protein]['peptides'][peptide.groups()[0]]['PSMs'].append(PSM_value)
                                
                                    #get phosphorylation status.
                                    self.proteins[last_protein]['peptides'][peptide.groups()[0]]['mods'].append(self.mods(line_columns[modifications_row], last_protein, peptide.groups()[0], mods_in_master, minimum_confidence))

                            else:

                                #this peptide does not exist already. add it and its position in protein, along with it's PSM count (if greater than cutoff), as well as mods if it has any, and number of protein groups.
                                
                                #only keep peptides with greater than minimum PSM:
                                if PSM_value > minimum_PSM:
                                
                                    position_start = None
                                    position_end = None
                                    
                                    if protein_groups_row is not None:
                                        protein_groups = int(line_columns[protein_groups_row])
                                        master_protein_accessions = line_columns[master_accessions_row]
                                    else:
                                        protein_groups = 1
                                        master_protein_accessions = proteins.Protein(last_protein).identifier

                                    if positions_in_master_proteins_row is not False:

                                        position_start = re.search(r'(?:\[)(\d+)', line_columns[positions_in_master_proteins_row]) #peptide start position
                                        position_end = re.search(r'(\d+)(?:\])', line_columns[positions_in_master_proteins_row])

                                        position_start = int(position_start.groups()[0])
                                        position_end = int(position_end.groups()[0])

                                    else:

                                        #the file does not have position info for the peptide so we need to search it.
                                        positions = self.get_peptide_position(last_protein, peptide.groups()[0], protein_groups, master_protein_accessions)

                                        if positions is not None:

                                            position_start = positions[0]
                                            position_end = positions[1]                                                

                                    self.proteins[last_protein]['peptides'][peptide.groups()[0]] = { 

                                        'start': position_start,
                                        'end': position_end,
                                        'PSMs': [ int(PSM_value) ], #put PSMs into a list. Each peptide is *slightly* different (unique peptides).
                                        'mods': [ self.mods(line_columns[modifications_row], last_protein, peptide.groups()[0], mods_in_master, minimum_confidence) ], #put mods into a similar list. the position of each item corresponds to the same PSM for that peptide, if the peptide repeats in the analysis. this function will add a mod if it detects one for that peptide.
                                        'groups': protein_groups
                                    }
                            
                            
                    else:
                        
                        print(f" !! !! Error: The file \"{self.filename}\" is not in a format that the program expects.")
                        
            count += 1
            
        return self #returns the analyzed file as an object
        
        self.filehandle.close()
        
    def get_phosphosites(self) :
        
        #   this gets KNOWN phosphosites (as opposed to phosphorylation with unknown residue). returns a dictionary with protein name as key and phosphosites in values
        #   with confidence values and which unique peptide it was identified in (in order from first to last, according to the excel sheet)
        
        self.phosphoproteins = {}
        
        for protein in self.proteins.keys():
            
            self.phosphoproteins[protein] = {} 
        
        for protein, protein_info in self.proteins.items():
            
            for peptide, peptide_info in protein_info['peptides'].items():
                
                position_in_list = 0
                protein_groups = int(peptide_info['groups'])
                
                for unique_peptide_mod in peptide_info['mods']:
                    
                    if unique_peptide_mod is not None:
                        #there is a modification on this peptide.
                        
                        for mod_type, mod_info in unique_peptide_mod.items():
                            
                            if mod_type == 'phos':
                                
                                #get amino acid and position, and PSMs for that unique peptide for each phosphorylation, and # protein groups
                                
                                for key, phos_info in mod_info.items():
                                
                                    #the values may not exist if the phosphorylation is known but the AA is not.
                                    if phos_info['pos'] is not None:
                                    
                                        aa = phos_info['AA']
                                        pos = phos_info['pos']
                                        conf = phos_info['conf']
                                        
                                    else:
                                        
                                        aa = phos_info['AA']
                                        pos = ''
                                        conf = None
                                        
                                    val = aa + str(pos)
                                    
                                    if val in self.phosphoproteins[protein].keys():
                                        
                                        #this amino acid is already in the dictionary. add the peptide counter, it's in a different unique peptide!
                                        
                                        self.phosphoproteins[protein][val]['unique_peptide_counter'].append(position_in_list)
                                        self.phosphoproteins[protein][val]['conf'].append(conf)
                                        
                                        
                                    else: 
                                        
                                        self.phosphoproteins[protein][val] = {

                                            'unique_peptide_counter': [position_in_list],
                                            'conf': [conf],
                                            'groups': protein_groups

                                        }
                                                                    
                    
                    position_in_list += 1
                    
        return self.phosphoproteins
                    

    def phospholist(self, protein) :
        
        #returns a simple list of the position of known phosphosites in the protein, or none if there aren't any.
        
        protein = protein
        dictionary = self.get_phosphosites()
        
        l = None
        
        if protein in dictionary.keys() :
            
            l = list()
        
            for key in dictionary[protein].keys():

                #only include the known sites, with known position
                search = re.search(r'(?:[A-Z])([\d]+)', key)

                if search is not None:

                    l.append(int(search.groups()[0]))
            
            if len(l) == 0:
                
                return None
            
            else:    
                
                l.sort()
                return l
        
        else:
            
            return None
        
    def mod_peptide_info(self, protein, mod_type, mod_position):
        
        #returns peptide (or peptides) info for a given mod, or None if that mod is not found.
        #returns a dictionary with peptide as key, and a nested dict as values, containing PSMs and confidences for each unique peptide.
        
        protein = protein
        mod_type = mod_type # 'phos' or 'acetyl'
        mod_position = mod_position #int, the amino acid position
        
        return_val = None
        
        if protein in self.proteins.keys():
            
            #the protein exists in the MS file (i.e. was detected):
            
            for peptide, peptide_info in self.proteins[protein]['peptides'].items():

                PSMs = peptide_info['PSMs'] #a list of the PSMs of each unique spectral peptide from this peptide
                mods = peptide_info['mods'] #a list, with each position containing a dictionary of mod info, respectively referring to each unique spectral peptide

                #make sure there are mods for this peptide. If none, skip.
                

                position_in_list = 0

                for unique_peptide_mod in mods:
                    
                    if type(unique_peptide_mod) is dict:

                        #there's a mod for this peptide.
                        #check to see if the mod is the mod we are interested in getting info for...

                        if mod_type in unique_peptide_mod.keys():

                            for mod_number, mod_info in unique_peptide_mod[mod_type].items():

                                if mod_info['pos'] is not None:

                                    if mod_info['pos'] == mod_position:

                                        #remember, there can be more than one peptide sequence containing the same phos, for example a peptide that has two lysines can be cleaved by trypsin in two spots.
                                        if return_val is None: #there is not a dictionary for the peptide, create one.

                                            return_val = {

                                                peptide: {

                                                    'PSMs': [PSMs[position_in_list]],
                                                    'conf': [mod_info['conf']]

                                                }

                                            }
                                        elif type(return_val) is dict:
                                            
                                            if peptide in return_val.keys():
                                                
                                                return_val[peptide]['PSMs'].append(PSMs[position_in_list])
                                                return_val[peptide]['conf'].append(mod_info['conf'])
                                                
                                            else:
                                                
                                                return_val[peptide] = {

                                                    'PSMs': [PSMs[position_in_list]],
                                                    'conf': [mod_info['conf']] 

                                                }
                                            
                                        
                    position_in_list += 1

        return return_val
    
    def total_PSMs_for_residue(self, protein, position):
        
        #total the number of spectral matches that cover a particular residue at a position.
        
        #RETURN: integer
        
        protein = protein
        position = int(position)
        
        if protein in self.proteins.keys():
            
            #this protein is detected in the mass spec.
            #loop through each peptide and see if 'position' is between start and end of peptide.
            
            total_PSMs = 0
            
            for peptide, peptide_info in self.proteins[protein]['peptides'].items():
                
                start   = peptide_info['start']
                end     = peptide_info['end']
                
                if position >= start and position <= end:
                    
                    #this peptide covered the residue we are interested in. add up the PSMs.
                    
                    #PSMs are stored in a list, and separate values indicate PSMs for a particular unique peptide. get the list, extract, and add up.
                    for PSMs in peptide_info['PSMs']:
                    
                        total_PSMs += PSMs
                        
            return total_PSMs
            
        else:
            
            return 0
        
    def get_peptide_position(self, protein, peptide, protein_groups = None, master_protein_accessions = None):
        
        #returns the start and end position of a peptide within a protein. Returns a tuple: (start, end) or None if not found.
        
        if protein_groups is None:
            protein_groups = 1
        
        match = re.search(r'(' + peptide + ')', proteins.Protein(protein).sequence)
        
        if match is not None:
            
            span = match.span()
            #the starting residue is position span[0]+1 but the ending residue is correct.
            return (span[0]+1, span[1])
            
        else:
            #we could not find the peptide... is it in another protein in the same group?
            
            if protein_groups > 1:
                
                did_we_find_a_match = False
                master_protein_accessions = master_protein_accessions.split('; ')
                for accession in master_protein_accessions:
                    match = re.search(r'(' + peptide + ')', proteins.Protein(accession).sequence)
                    
                    if match is not None:
                        span = match.span()
                        return (span[0]+1, span[1])
                        did_we_find_a_match = True
                        break #stop the loop
                        
                if did_we_find_a_match is False:
                    return None
            
            else:
                return None
        