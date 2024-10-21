import sys, os, re, regex, requests
regex.DEFAULT_VERSION = regex.VERSION1
#regex is a third party regex which allows overlapping regular expressions, unlike Python's re module.

#create a master list (dictionary) that has all protein info in it from yeast SGD:
global master_list
master_list = {}

FASTA_file = 'orf_trans_all.fasta'  #FASTA file with all yeast proteins (ORFs, trans and dubious), downloaded from SGD.
file_path = os.path.dirname(os.path.realpath(__file__)) #the path to this file right here.

filehandle = open(file_path + '/' + FASTA_file, 'r')
lines = filehandle.readlines() #read line by line...

protein_count = 0
last_key = None

for line in lines:

    #read the FASTA file line-by-line
    #look for the identifier (Y ID)

    #if the line is the protein identifier information in the FASTA file:
    match = re.search(r">(\w+[-]*\w+)", line)

    if match is not None:

        #extract gene name to set the dictionary key
        gene_name = re.search(r"(?:>[A-Z0-9\-]+\s)([A-Z0-9\-]+)", line)
        gene_name = gene_name.group(1)

        #identifier name
        identifier_name = match.group(1)

        master_list[gene_name] = { 'seq': '', 'identifier': identifier_name }

        last_key = gene_name

        protein_count += 1

    #if the line is a protein sequence following the protein identifier:

    protein_seq = re.search(r"^([A-Z])+|(?:[\*])+$", line)

    if protein_seq:
        if protein_seq.group(0) != '*':
            #it has protein sequence. add it to the dictionary.
            if last_key is not None:
                master_list[last_key]['seq'] = master_list[last_key]['seq'] + protein_seq.group(0) #append the sequence to the existing sequence.

#now every gene is in a dictionary named master_list with gene name (or, if no gene name, identifier) as key!

class Protein:
    
    def __init__(self, name) :
        
        self.name = name.upper()#make uppercase in case user does not do so...
        self.name.strip()
        self.sequence = None
        self.length = None
        self.identifier = None
        
        #CONSENSUS SITES, in a dictonary, of each known kinase:
        #important: the phosphorylated residue should always be in a match group (i.e. in parentheses) and the rest of the regular expression should not be.
        
        self.consensus = {
        
            'Mps1': '[DENC][A-Z]([ST])[^DP]',
            'Mps1 (if -2 phos)': '[ST][A-Z]([ST])[^DP]',
            'Ipl1': '[RK][A-Z]([ST])[^P]',
            'Cdk1': '([ST])[P](?:[A-Z][KR])?', #sometimes the x-KR at the end is optional.
            'Cdc5': '[DENC][A-Z]([ST])(?:[^PDEN]|[FMYI])',
            'Cdc5 (if +1 phos)': '[DENC][A-Z]([ST])[ST]'
            
        }
        
        self.master_list = master_list
        
        #look for the protein in the master list.
        try:
            master_list[self.name]
        except:
            
            self.sequence = None #default these to None.
            self.length = None
            self.identifier = None
            
            # maybe the user supplied the accession/identifier
            # Fetch the gene name that I have in my list:
            
            found_the_gene = False
            
            for gene_name, gene_info in master_list.items():

                if gene_info['identifier'] == self.name:
                        
                    self.sequence = master_list[gene_name]['seq']
                    self.length = len(master_list[gene_name]['seq'])
                    self.identifier = master_list[gene_name]['identifier']
                    self.name = gene_name
                    
                    found_the_gene = True
                    
            if found_the_gene is False:
                
                #we failed to find the gene. search the SGD:
                
                r = requests.get(f"https://www.yeastgenome.org/backend/locus/{self.name}")
            
                if r.status_code == 200:
                
                    #let's set the self.name as the identifier/accession
                    response = r.json()
                    self.name = response['gene_name']
                    self.identifier = response['format_name']
                    self.sequence = master_list[self.name]['seq']
                    self.length = len(master_list[self.name]['seq'])

        else:
            self.sequence = master_list[self.name]['seq']
            self.length = len(master_list[self.name]['seq'])
            self.identifier = master_list[self.name]['identifier']
            
    def subsequence(self, start_pos, end_pos) :
        
        start_pos = start_pos - 1
        end_pos = end_pos - 1
        
        substr = self.sequence[start_pos:end_pos]
        
        return substr
    
    def find_sites(self, reg_expression) :
        
        #searches a particular protein for a regular expression pattern, and returns some values
        #RETURN self: multiple variables pertaininig to the search. These variables are lists. Empty lists are returned if no sites are found.
        #the length of the positions list is equal to the number of matches.
        #IMPORTANT: the reg_expression must contain a single match group (i.e. a residue in parentheses) - this pertains to the phosphoresidue. The rest of the sequence should not be included in the match group.

        reg_expression = reg_expression.upper()

        does_anything_match = re.search(reg_expression, self.sequence) # a quick check to make sure there is a match.

        self.positions = []
        self.subsequences = []
        self.phosphoresidues = []

        if does_anything_match is not None:

            matches = regex.finditer(reg_expression, self.sequence, overlapped=True) #uses third party regex function to allow overlapping finds. returns an interable object containing all matches       

            for single_match in matches:

                self.positions.append(single_match.span()) #put start and ending tuples in the list!
                self.subsequences.append(self.subsequence(single_match.span(0)[0]+1, single_match.span(0)[1]+1))

                #get the phosphoresidue:

                self.phosphoresidues.append(single_match.span(1)[1])

        return self
        
        re.purge()   
    
    def kinase (self, AA_site = None) :
        
        #important: returns a dictionary of predicted phosphorylation sites for the protein, UNLESS AA_site is supplied, then it returns a dictionary with 'kinase': 'matched consensus' for that site, or None if no known kinase is predicted.
        
        AA_site = int(AA_site) #(position within the protein)
            
        AA_site_dictionary = {}

        for enzyme, consensus_sequence in self.consensus.items():

            does_it_have_a_match = re.search(consensus_sequence, self.sequence)

            if does_it_have_a_match is not None:

                #a kinase consensus site exists. let's find all of them, including overlapping ones:
                
                matches = regex.finditer(consensus_sequence, self.sequence, overlapped=True)

                for match in matches: 
                    
                    site_position = int(match.span(1)[1])
                    sequence_that_matched = match.group(0)
                    
                    #check to see if this AA site has been added to the AA_site_dictionary yet:
                    
                    if site_position not in AA_site_dictionary.keys():
                        AA_site_dictionary[site_position] = { 'enzyme': [enzyme], 'consensus': [sequence_that_matched] }
                        
                    else:
                        #the AA site already exists, and is predicted to be more than one kinase, so add it.
                        AA_site_dictionary[site_position]['enzyme'].append(enzyme)
                        AA_site_dictionary[site_position]['consensus'].append(sequence_that_matched)
            
        #sort by residue order:
        
        AA_sites_list = list(AA_site_dictionary.keys())
        AA_sites_list.sort()
        AA_sites_ordered = {}
        
        for site in AA_sites_list:
            
            AA_sites_ordered[site] = AA_site_dictionary[site] #another dictionary with kinase and sequence info for that site."""
             
        if AA_site is not None:
            
            #return dictionary with 'kinase': and 'consensus'. Keep in mind more than one kinase could be predicted.
            
            if len(AA_site_dictionary) == 0:
                
                return None #no known enzyme is predicted to phosphorylate this site.
            
            else:
                
                if AA_site in AA_site_dictionary.keys():
                
                    return AA_site_dictionary[AA_site]
                
                else: 
                    
                    return None
            
        else:
            
            #the user did not supply a site, return all consensus sites for the protein (i.e. predicted_sites dictionary):
            
            return AA_site_dictionary
        
        re.purge()
        
    
        
    def surrounding_seq(self, AA_site):
        
        # generate the surrounding sequence (10 amino acids total) of a particular phosphosite
        # this function will add X's in place of amino acids if the site is near the N or C terminus and a 10-AA sequence cannot be generated.
        # that ensures the length is always 10.
        
        #the user can supply amino acid with position or just position.
        AA_site = str(AA_site) #regex gives error if not a string first..
        AA_site = re.search(r'(\d+)', AA_site)
        
        if AA_site is not None:
            
            AA_site = AA_site.group()
            AA_site = int(AA_site)
            
        AA_site_substr = AA_site - 1
        minus_5 = AA_site_substr - 5
        plus_4 = AA_site + 4
        minus_5_abs = None #to add X's to the sequence if the sequence runs off the beginning
        plus_4_abs = None #to add X's to the sequence if the sequence runs off the end.
        surrounding_seq = ''
        
        if minus_5 < 0:
            minus_5_abs = abs(minus_5)
            minus_5 = 0
        
        if plus_4 > len(self.sequence):
            plus_4_abs = abs(plus_4-len(self.sequence))
            plus_4 = len(self.sequence)
            
        #add X's to the beginning or end of subsequence in the case it runs off the end, so that it is always length of 10 amino acids.
        
        if minus_5_abs is not None:
            for i in range(0, minus_5_abs):
                surrounding_seq += 'X'
        
        surrounding_seq += self.sequence[minus_5:plus_4]
        
        if plus_4_abs is not None:
            for i in range(0, plus_4_abs):
                surrounding_seq += 'X'
                
        #make sure the site actually exists...
        if AA_site > len(self.sequence) or AA_site < 1:
            return None
        else:
            return surrounding_seq
        
    def get_all_ser_thr(self):
        
        #returns a list of sites that are serine or threonine:
        search = regex.search('([ST])', self.sequence)
        
        if search is not None:
            all_ST = list()
            searches = regex.finditer(r'([ST])', self.sequence)
            
            for match in searches:
                pos = match.span()[1]
                all_ST.append(pos)
            
            return all_ST
            
        else:
            
            return list(None)