import kinetochore as KT, json, time, requests, sys, proteins as p

print('-----------> Updating SGD data... please wait...')

# this file gets all reported phosphorylation (on SGD) for each kinetochore protein. 
# stores the information obtained from the SGD as a JSON file (published_phos.json) and will update every 6 months or so.

#first, check the last time the json file was updated. We want to update upon the following conditions:
#   --> A new protein is added or subtracted from the kinetochore protien dictionary (within the kinetochore.py file) = update only that one protein.
#   --> it has been more than 6 months since the last update occured (update all proteins)

class SGD :
    
    def __init__(self):
    
        #open the JSON file:
    
        json_file = open('published_phos.json', 'r')
        read = json_file.read()
        self.json_dict = json.loads(read) #the contents of the json file, converted to a python dictionary.
        json_file.close()
        
        unix = time.time()
    
        self.json_proteins = self.json_dict['proteins']

        # 1) check to see when the last update was.
        
        #time since last update (in seconds):
        time_since = unix - self.json_dict['last_update']
        six_months = 60*60*24*7*4*6 # roughly six months, in seconds

        if time_since > six_months:
            #we should update the JSON file completely.
            self.update_all()
        
        #Check the protein list in the kinetochore.py file for any proteins that may have been added by the user since the last update:
        KT_accessions_dict = KT.Kinetochore().accessions
        KT_accessions = KT_accessions_dict.items()
    
        for protein, accession in KT_accessions:

            #Add proteins to the json file if they do not exist in the json list but exist in the kinetochore protein list:
            if accession not in self.json_proteins.keys():
                #the protein does not exist in the JSON file. Add it.
    
                req = self.request(accession)
                self.update(accession, req)
            
                self.json_dict['proteins'][accession] = req
    
        self.json_proteins = self.json_dict['proteins'] #a dictionary containing accession: PTM data
    
    def request(self, accession):

        r = requests.get(f"https://www.yeastgenome.org/backend/locus/{accession}/posttranslational_details")
        if r.status_code != 200:
            print('Warning: Something went wrong searching for this gene on SGD. Check your internet status or the gene name and try again.')
            sys.exit()
        else:
            #the request went through. let's put it in a format that allows to append to the phosphorylation json file.
            response = json.loads(r.text)

            self.output = {}
    
            for each_site in response:
                # each_site is a dictionary that contains info per residue that has been discovered to have a PTM.
                # site_index is the residue number
                # only include phosphorylated residues:

                if each_site['type'] == 'phosphorylated residue':

                    #there could be multiple reports of phosphorylation on this residue.

                    site_index = int(each_site['site_index'])

                    if site_index not in self.output.keys():
                        #the residue is not already in the output dictionary. add it, and it's reference, into a list.

                        self.output[site_index] = {

                            'refs': [each_site['reference']['display_name']],
                            'ref_PMIDs': [int(each_site['reference']['pubmed_id'])]

                        }

                    else:
                        #the residue is already in the output dictionary. update the references to include this reference.

                        self.output[site_index]['refs'].append(each_site['reference']['display_name'])
                        self.output[site_index]['ref_PMIDs'].append(int(each_site['reference']['pubmed_id']))

            
        time.sleep(2) #give the SGD 2 seconds of peace and quiet. Helps not flood requests if doing iteratively.
    
        return self.output
    
    def update(self, accession, req_output):
        #update a single protein.
        
        #protein name:
        name = p.Protein(accession).name
        print(f' - - - - - > Updating SGD data to include {accession} / {name}...')
        
        #do not update the time, since we are just updating one.
    
        self.json_dict['proteins'][accession] = req_output
    
        json_file_write = open('published_phos.json', 'w')
        json_file_write.write(json.dumps(self.json_dict))
        json_file_write.close()
        
    def update_all(self):
        #updates all accessions within the published_phos.json file.
        #we want to create a temp published_phos.json file. Update that. And then delete the old file, and rename the new file as the original.
        
        #json_dict['last_update'] = int(time.time()) #update with current unix timestamp.
        new_json_dict = {
        
            'last_update': None, #fill this in later
            'proteins': {} #empty
            
        }
        
        name = p.Protein(accession).name
        
        number_of_proteins = len(self.json_proteins.items())
        count = 1
        
        for accession, phos_info in self.json_proteins.items():
            new_json_dict['proteins'][accession] = self.request(accession)
            
            print(f' - - - - - > Updating {accession} / {name} phospho-sites from SGD... {count}/{number_of_proteins}')
            count += 1
                
        new_json_dict['last_update'] = int(time.time()) #update with current unix timestamp.
        
        update_file = open('published_phos.json', 'w')
        update_file.write(json.dumps(new_json_dict))
        update_file.close()
        
    def get_phos(self, accession):
        
        #returns a list of reported phosphosites of an accession. simple!
        
        if accession in self.json_proteins.keys():
            #the protein already exists in published_phos.json
        
            list_of_sites = list(self.json_proteins[accession].keys())
        else:
            #the protein does not exist in the file.. let's add it to the file.
            request = self.request(accession)
            self.update(accession, request)
            
            list_of_sites = list(request.keys())
        
        return list_of_sites
    
#run the SGD class:
SGD()