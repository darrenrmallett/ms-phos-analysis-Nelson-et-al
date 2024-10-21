import os, sys, regex, re, mass_spec, proteins as p, xlsxwriter, time, kinetochore as KT, sgd_phosphorylation as sgd, json
from xlsxwriter.utility import xl_col_to_name, xl_rowcol_to_cell

current_output_path = os.path.dirname(os.path.realpath(__file__)) + '/output'

kinetochore = KT.Kinetochore()
unix = int(time.time()) #UNIX timestamp.

#compares files.
#remember, master_list is the master protein list of all yeast proteins.

def output_file_names(prefix, files, cond) :
    #takes prefix and adds filenames to it so it's easy for the user to know which files were used to generate the output file
        
    file_string = str(unix)
    cond = cond
    return prefix + '_' + file_string + '__' + cond + '.xlsx'

def compare_all_phosphorylation(files = dict) :
    
    files = files #files = a dictionary with MassSpec objects as values.
    
    #goal: generate a dictionary with protein name as key and a new dictionary as a value, with keys in this nested dictionary being the file, and the values as phosphosites.
    
    all_proteins = {} #a dictionary with all proteins detected between the MS files that are being compared.
    all_phosphosites = {} #a dictionary with all proteins as key and a nested dictionary with keys being all phosphorylation sites between all files
    
    for file_name, ms_data in files.items():
        
        for protein_name in ms_data.proteins.keys():
        
            all_proteins[protein_name] = {}
            all_phosphosites[protein_name] = {}
            
    #now the all_proteins and all_phosphosites dictionaries have all proteins detected across all MS files in it.
    
    for protein_name in all_proteins.keys():
        
        #loop through each protein, put in value of dictionary a new dictionary containing file_number as key, and list of phosphosites as values
        #also put each phosphosite as key in nested dictionary in all_phosphosites
        
        file_counter = 0
        for file_name, ms_data in files.items():
                            
            l = ms_data.phospholist(protein_name)
                
            all_proteins[protein_name][file_counter] = l
            
            file_counter += 1
            
            #put each phosphosite as key in nested dictionary in all_phosphosites. Put value as the max protein groups. Default to 0.
            
            if l is not None:
                
                for site in l:
                    
                    all_phosphosites[protein_name][site] = 0 #protein groups will be filled below.
                    protein_grps = 0
                    
                    for ms_data in files.values():
                        
                        phosphosites = ms_data.get_phosphosites()
                        
                        if protein_name in phosphosites.keys():
                            
                            site_with_aa = p.Protein(protein_name).sequence[site-1]
                            site_with_aa = site_with_aa + str(site)
                            
                            if site_with_aa in phosphosites[protein_name].keys():
                                protein_grps = max(protein_grps, phosphosites[protein_name][site_with_aa]['groups']) #take the max value. I dont know if protein groups can be 1 in one file and 2 in another. probably not but it's better this way.
                    all_phosphosites[protein_name][site] = protein_grps
                    
    
    #generate output file
    
    output_header = [
    
        'Accession',
        'Protein',
        'Site',
        'Protein groups',
        'Reported on SGD?',
        'Predicted kinase (regex)',
        'Kinase consensus (regex)',
        'Match',
        'Surrounding seq (+/- 5)'

    ]
    
    global_fields = len(output_header) #this is a list of just global fields, not including file-specific fields.
    
    all_files = list()
     
    for file in files.keys():
        
        file = file.split('.')
        file.pop()
        
        file_name_without_extension = ''
        
        for list_item in file:
            
            file_name_without_extension += list_item
        
        all_files.append(file_name_without_extension)
        
    file_specific_fields = [
     
        'detected?',
        'Total PSMs covering residue',
        'Total Phospho-PSMs',
        'Peptide(s)',
        '% Confidence (PSMs)'
        
    ]
    
    field_specific_bg_color = [
     
        '#0A2647',
        '#144272',
        '#205295',
        '#2C74B3',
        '#004254',
        '#006071',
        '#007B91',
        '#0D7C99',
        '#4C0033',
        '#790252',
        '#AF0171',
        '#E80F88',
        '#2D132C',
        '#801336',
        '#C72C41',
        '#EE4540'
        
    ]
        
    
    num_of_files = str(len(files.keys()))
    
    for file in files.keys():
        min_PSM = int(files[file].minimum_PSM)
        min_PSM = min_PSM+1
        min_conf = files[file].minimum_conf
        cond = f"min_PSM={min_PSM}__min_confidence={min_conf}"
        break;
    
    workbook = xlsxwriter.Workbook(current_output_path + '/' + output_file_names('compare_all_phos__' + num_of_files + '_cond_', files, cond), {'constant_memory': True}) #constant_memory = writes row-by-row to dump memory after writing
    worksheet = workbook.add_worksheet('Detected phos') #first worksheet.
    worksheet.freeze_panes(1, 0)
    
    header_format = {}
    
    how_many_columns_will_there_be = len(output_header) + (len(file_specific_fields*len(files.keys())))
    color_value = 0
    file_count = 1
    for i in range(0, how_many_columns_will_there_be):
        
        header_format[i] = workbook.add_format()
        header_format[i].set_align('center')
        header_format[i].set_align('vcenter')
        header_format[i].set_text_wrap()
        header_format[i].set_color('#FFFFFF')
        header_format[i].set_border(1)
        header_format[i].set_border_color('#CACACA')
        
        if i <= (len(output_header)-1):
        
            header_format[i].set_bg_color('#000000')
        
        else:
            
            header_format[i].set_bg_color(field_specific_bg_color[color_value])
            
            if color_value == len(field_specific_bg_color)-1 or file_count == len(files.keys()): #max color reached, or max file reached. reset back to the first color
                color_value = 0
            else:            
                color_value += 1
                
            if file_count == len(files.keys()):
                file_count = 1
            else:
                file_count += 1
    
    for field in file_specific_fields:
        
        for file in all_files:
            
            output_header.append(file + ' - ' + field)
            
        
    output_header_string = ''
    
    count = 1
    for column in output_header:
        
        output_header_string += column
        
        if count != len(output_header):
            output_header_string += '\t' #no tab characters after last column. only want \n
        else:
            output_header_string += '\n'
        
        count += 1
    
    #conditional formatting for the rest of the Excel data:
    
    # Light red fill with dark red text.
    red_format = workbook.add_format({      'bg_color':   '#FFC7CE',
                                            'font_color': '#9C0006'})

    # Light yellow fill with dark yellow text.
    yellow_format = workbook.add_format({   'bg_color':   '#FFEB9C',
                                            'font_color': '#9C6500'})

    # Green fill with dark green text.
    green_format = workbook.add_format({    'bg_color':   '#C6EFCE',
                                            'font_color': '#006100'})
    
    # Dark blue fill with white text.
    blue_format = workbook.add_format({    'bg_color':   '#2154A5',
                                            'font_color': '#FFFFFF',
                                            'bold': True})
    
    #for the trigger warning, if there is a warning, make it red.
    
    column = xl_col_to_name(8)

    worksheet.conditional_format(f'{column}1:{column}1048576', { #entire column = H1:H1048576 (the last possible row is 1048576).

        'type': 'cell',
        'criteria': 'equal to',
        'value': '"Warning"',
        'format': red_format
    })
    
    worksheet.conditional_format(f'E1:E1048576', { #entire column = H1:H1048576 (the last possible row is 1048576).

        'type': 'cell',
        'criteria': 'equal to',
        'value': '"Reported"',
        'format': blue_format
    })
    
    worksheet.conditional_format(f'D2:D1048576', { #entire column = H1:H1048576 (the last possible row is 1048576).

        'type': 'cell',
        'criteria': 'greater than',
        'value': '1',
        'format': yellow_format
    })
     

    #for the 'detected' columns, I want conditional formatting. The number of columns will change with # of files. Start at column H (column 7).
    
    col = global_fields #starting column for file-specific outputs.
    
    for f in all_files: 
        
        column = xl_col_to_name(col)

        worksheet.conditional_format(f'{column}1:{column}1048576', { #entire column = H1:H1048576 (the last possible row is 1048576).

            'type': 'cell',
            'criteria': 'equal to',
            'value': '"NO"',
            'format': red_format
        })

        worksheet.conditional_format(f'{column}1:{column}1048576', {

            'type': 'cell',
            'criteria': 'equal to',
            'value': '"YES"',
            'format': green_format
        })

        worksheet.conditional_format(f'{column}1:{column}1048576', {

            'type': 'cell',
            'criteria': 'equal to',
            'value': '"No PSMs"',
            'format': yellow_format
        })
        
        col += 1
    
    #write the header string to the Excel worksheet:
    
    row = 0
    col = 0
    
    for header_item in output_header:
        
        worksheet.write(row, col, header_item, header_format[col])
        
        col += 1 #move to next column every header item
        
    #set row for output strings:
    row = 1
    
    #now generate the output string, a strenuous task. Output line-by-line.:
    
    output_list = list()
    kinase_statistics = {} #for the second worksheet - statistics for each kinase.
    
    for protein, sites in all_phosphosites.items():
        
        all_sites = list(sites.keys())
        all_sites.sort()
        
        for site in all_sites:
            
            accession = p.Protein(protein).identifier
            protein = protein
            site_substr = site - 1
            amino_acid = p.Protein(protein).sequence[site_substr]
            AA_site = amino_acid + str(site)
            protein_groups = all_phosphosites[protein][site]
            
            if str(site) in sgd.SGD().get_phos(accession):
                reported_in_sgd = "Reported"
            else:
                reported_in_sgd = " "
            
            predicted_kinase = p.Protein(protein).kinase(site) #returns a dictionary with enzyme: match_sequence, or None if no known enzyme.
            
            predicted_kinase_string = ''
            kinase_match = ''
            kinase_consensus = ''
            
            if predicted_kinase is None:
                predicted_kinase_string = '' #leave blank, we don't know.
            else:
                #it could return more than one kinase.
                count = 1
                if protein not in kinase_statistics.keys():
                    kinase_statistics[protein] = {} #add the protein to the statistics dictionary, if it's not there already.
                
                kinase_statistics[protein][site] = [] #add the site to the statistics dictionary, with a list of all kinases that match this seq.
                
                for kinase in predicted_kinase['enzyme']:
                    
                    predicted_kinase_string += kinase
                    kinase_consensus += p.Protein(protein).consensus[kinase]
                    
                    kinase_statistics[protein][site].append(kinase) #now this dictionary is all set for the second worksheet (tab).
                    
                    if len(predicted_kinase['enzyme']) >=2 and count != len(predicted_kinase['enzyme']):
                        predicted_kinase_string += '; '
                        kinase_consensus += '; '
                    
                for match_seq in predicted_kinase['consensus']:
                    
                    kinase_match += match_seq
                    
                    if len(predicted_kinase['consensus']) >=2 and count != len(predicted_kinase['consensus']):
                        kinase_match += '; '
                    
                    count +=1
                    
            site_minus_1 = site-1
            five_minus = (site-1)-5
            five_plus = (site-1)+6
            
            if five_minus < 0:
                five_minus = 0 #we don't want to go below zero, the starting methionine.
            
            if five_plus > (p.Protein(protein).length-1):
                five_plus = p.Protein(protein).length-1
            
            surrounding_seq = p.Protein(protein).sequence[five_minus:site] + '*' + p.Protein(protein).sequence[site:five_plus]
            surrounding_seq = '[' + str(five_minus+1) + ']-' + surrounding_seq + '-[' + str(five_plus) + ']'            
            
            #now let's put together the output_list so far...
            output_list = [
             
                accession,
                protein,
                AA_site,
                protein_groups,
                reported_in_sgd,
                predicted_kinase_string,
                kinase_consensus,
                kinase_match,
                surrounding_seq
                
            ]
            
            
            non_file_spec_output_list_count = len(output_list)
            
            #Now the file-specific information for each site...
            
            file_specific_output = {}
            
            for file_name, ms_data in files.items():
                
                file_specific_output[file_name] = {}
            
                # ------> Can we detect the phosphorylation in this file? 
                
                mod_peptide_info = ms_data.mod_peptide_info(protein, 'phos', site)
                
                if mod_peptide_info is None:
                    
                    file_specific_output[file_name]['detected'] = "NO"
                
                else:
                    
                    file_specific_output[file_name]['detected'] = "YES"
                    
                # ------> Get peptides for this phosphorylation, and confidences, and PSMs. Also get all PSMs that detect the phosphorylation. Also get total unique phosphopeptides.
                
                peptides = ms_data.mod_peptide_info(protein, 'phos', site)
                
                file_specific_output[file_name]['peptides'] = ''
                file_specific_output[file_name]['confidences_PSMs'] = '' #condense all the confidences in one string, write as: Conf1 (PSMs1); Conf2 (PSMs2)
                file_specific_output[file_name]['phospho_PSMs'] = 0
                file_specific_output[file_name]['total_unique_phosphopeptides'] = 0
                
                if peptides is not None:
                
                    n = 1

                    for peptide, values in peptides.items():

                        file_specific_output[file_name]['peptides'] += peptide + '; '

                        key = 0

                        for conf in values['conf']:
                            #values['conf'] is a list of confidence values. create the string from this, and the corresponding PSM for that confidence.

                            file_specific_output[file_name]['confidences_PSMs'] += str(conf) + ' (' + str(values['PSMs'][key]) + '); '
                            file_specific_output[file_name]['phospho_PSMs'] += values['PSMs'][key]

                            key += 1
                            
                            #total unique phosphopeptides should be equal to key at this point!

                            file_specific_output[file_name]['total_unique_phosphopeptides'] = key
            
                # ------> Get total PSMs that covered the residue (NOT just phosphorylation):
                
                file_specific_output[file_name]['total_PSMs_for_residue'] = ms_data.total_PSMs_for_residue(protein, site)
                
                if  file_specific_output[file_name]['total_PSMs_for_residue'] == 0:
                    
                    file_specific_output[file_name]['detected'] = "No PSMs"

                
                # Now write to the Excel file:
                
            #now add all those file_specific_outputs to the output_list:
            
            file_specific_parameters = [

                'detected',
                'total_PSMs_for_residue',
                'phospho_PSMs',
                'peptides',
                'confidences_PSMs'

            ]
            
            for parameter in file_specific_parameters:
                
                for file_name in file_specific_output.keys():
                    
                    output_list.append( file_specific_output[file_name][parameter] )
            
            #figure out how many columns there will be:
            
            #col_count = (len(output_list) + (len(file_specific_parameters)*len(files.keys())))-1
            col_count = len(output_list)
            
            formatting = {}
            
            file_spec_col_count = 1
            file_count = 0
            alternate_bg = 0
            
            for i in range(0, col_count):
            
                formatting[i] = workbook.add_format()
                formatting[i].set_border_color('#C0BEBE')
                formatting[i].set_border(1) 
                
                #column-specific formats:
                
                center_cols = [0,1,2,3,4,5,6,7,8]
                if i in center_cols:
                    
                    formatting[i].set_align('center')
                    
                courier_font = [6,8]
                if i in courier_font:
                    
                    formatting[i].set_font_name('Courier')
                    
                #get file-specific columns:
                
                if i+1 == non_file_spec_output_list_count:
                    
                    formatting[i].set_right(2)
                    formatting[i].set_right_color('#959595')
                
                elif i+1 > non_file_spec_output_list_count:
                    
                    #this is a file-specific column
                    file_count += 1
                    
                    #they are all centered.
                    formatting[i].set_align('center')
                    
                    if alternate_bg == 0:

                        formatting[i].set_bg_color('#F4F9FF')

                    else:

                        formatting[i].set_bg_color('#EAF4FF')
                        
                    if file_spec_col_count == 4 or file_spec_col_count == 5:

                        #the last two file-specific outputs I want align left:
                        formatting[i].set_align('left')
                    
                    if file_count == len(files.keys()):
                        
                        file_count = 0 #reset.
                        file_spec_col_count += 1
                        
                        #we want a thicker border to the left.
                        
                        formatting[i].set_right(2)
                        formatting[i].set_right_color('#959595')
                        
                        if alternate_bg == 0:
                            alternate_bg = 1
                        elif alternate_bg == 1:
                            alternate_bg = 0
                        
                else:
                    formatting[i].set_border(1) 
                    

            #set column widths (of global columns i.e. non-file-specific columns):
            worksheet.set_column(0, 0, 14) #column A = 14
            worksheet.set_column(1, 1, 12) #column B
            worksheet.set_column(2, 2, 11) #column C
            worksheet.set_column(3, 3, 7) #column C
            worksheet.set_column(4, 4, 11) #column D
            worksheet.set_column(5, 5, 20) #column E
            worksheet.set_column(6, 6, 30)
            worksheet.set_column(7, 7, 15)
            worksheet.set_column(8, 16, 35)
            
            #for the file-specific columns, make sure they all get the same formatting.
            
            c = global_fields #the starting column of the file-specific output columns.
            
            for param in file_specific_parameters:
                
                
                for f in all_files:
                    
                    if param == 'detected' or param == 'total_PSMs_for_residue' or param == 'phospho_PSMs' or param == 'total_unique_phosphopeptides':
                        
                        #formatting[c] = workbook.add_format({'align': 'center'})
                        worksheet.set_column(c, c, 12)
                        
                    if param == 'peptides' or param == 'confidences_PSMs':
                    
                       #formatting[c] = ''
                        worksheet.set_column(c, c, 35)
                        
                    
                    c += 1
            
            col = 0
            
            output_string = ''
                
            for output_item in output_list:

                if col in formatting.keys():
                
                    output_string += str(output_item) + '; '
                    worksheet.write(row, col, output_item, formatting[col])
                    
                else:
                    
                    worksheet.write(row, col, output_item)

                if col == len(output_list)-1:
                    
                    row += 1 #go to next row!
                    
                col += 1

    # # # second worksheet - all phos summary # # #
    
    worksheet2 = workbook.add_worksheet('Kinase stats') #kinase statistics worksheet.
    
    # To read the JSON file containing all published phosphorylation sites already known:    
    json_file = open('published_phos.json', 'r')
    read = json_file.read()
    json_dict = json.loads(read) #the contents of the json file, converted to a python dictionary.
    json_file.close()
    
    #pre-set vars:
    total_phosphosites = 0 #count all phosphosites.
    total_phosphosites_with_kinase_assigned = 0 #count all sites that have been assigned at least one kinase.
    
    for protein, site_dict in all_phosphosites.items():
        total_phosphosites += len(site_dict.keys())
    
    for protein, site_dict in kinase_statistics.items():
        total_phosphosites_with_kinase_assigned += len(site_dict.keys())
    
    #create a dictionary with kinases as key and dict as value with parameters..
    site_stats = {}
    site_stats_incl_priming = {} #combined if sites are primed or not.
    for kinase in p.Protein('NDC80').consensus.keys(): #randomly chose Ndc80 just to get the consensus dictionary since a protein has to be called.
        site_stats[kinase] = {
            'unique': 0,
            'num_reported_sgd': 0,
            'shared': 0
        }
        
        #to generate priming site_stats:
        search = regex.search(r'\(if', kinase)
        if search is None:
            site_stats_incl_priming[kinase] = {
                'unique': 0,
                'num_reported_sgd': 0,
                'shared': 0
            }
        
    for protein, site_info in kinase_statistics.items():
        
        accession = p.Protein(protein).identifier
        
        for site, kinase_list in site_info.items():
            #for the sites that are separate from the primed sites
            for kinase in site_stats.keys():
                #loop through each kinase in the site_stats dict and determine whether this site is in this list.
                
                if kinase in kinase_list:
                    #yes, it's in the list... is it the only kinase to match this site?
                    if len(kinase_list) == 1:
                        #yes it is. add it to the unique count.
                        site_stats[kinase]['unique'] += 1
                    else:
                        #it is shared with another kinase.
                        site_stats[kinase]['shared'] += 1
                        
                    #how many of these sites are reported in the SGD?
                    if str(site) in sgd.SGD().get_phos(accession):
                        site_stats[kinase]['num_reported_sgd'] += 1
            
            #for the sites that are combined (primed and not primed):
            new_list = {} #put unique kinases in this dict as keys after the priming has been removed from kinase name (to generate a unique list)

            for kinase in site_stats_incl_priming.keys():
                #loop through each kinase in the site_stats dict and determine whether this site is in the list
                
                #let's remove duplicates if there exists the main kinase plus its primed site.
                for kinase_name in kinase_list:
                    match = regex.search(r'^(\w+)', kinase_name) #just pull out the kinase. Will match primed or not.
                    new_list[match[1]] = 0 #put in dictionary as key.
                    
                if kinase in new_list.keys():
                    #yes, it's in the list. is it the only kinase in this list?
                    if len(new_list.keys()) == 1:
                        #yes it is. add it to the unique count.
                        site_stats_incl_priming[kinase]['unique'] += 1
                    else:
                        #it is shared with another kinase.
                        site_stats_incl_priming[kinase]['shared'] += 1
                
                    #how many of these sites are reported in the SGD?
                    if str(site) in sgd.SGD().get_phos(accession):
                        site_stats_incl_priming[kinase]['num_reported_sgd'] += 1
                        
    worksheet2.set_column(0, 0, 25) #column A = 15
    
    header_format = workbook.add_format()
    header_format.set_align('center')
    header_format.set_align('vcenter')
    header_format.set_text_wrap()
    
    #write the HEADER (row, col):
    worksheet2.write(0, 0, 'Kinase (priming separate)', header_format)
    worksheet2.write(0, 1, '# unique', header_format)
    worksheet2.write(0, 2, '# shared', header_format)
    worksheet2.write(0, 3, 'total', header_format)
    worksheet2.write(0, 4, '# report on SGD', header_format)
    
    row = 1
    
    for kinase, stats in site_stats.items():
        
        total = stats['unique'] + stats['shared']
        
        worksheet2.write(row, 0, kinase, '')
        worksheet2.write(row, 1, stats['unique'], '')
        worksheet2.write(row, 2, stats['shared'], '')
        worksheet2.write(row, 3, total, '')
        worksheet2.write(row, 4, stats['num_reported_sgd'], '')
        
        row += 1   

    #write the next header
    worksheet2.write(row, 0, 'Kinase (priming included)', header_format)
    worksheet2.write(row, 1, '# unique', header_format)
    worksheet2.write(row, 2, '# shared', header_format)
    worksheet2.write(row, 3, 'total', header_format)
    worksheet2.write(row, 4, '# report on SGD', header_format)
    
    row += 2
    
    for kinase, stats in site_stats_incl_priming.items():
        
        total = stats['unique'] + stats['shared']
        
        worksheet2.write(row, 0, kinase, '')
        worksheet2.write(row, 1, stats['unique'], '')
        worksheet2.write(row, 2, stats['shared'], '')
        worksheet2.write(row, 3, total, '')
        worksheet2.write(row, 4, stats['num_reported_sgd'], '')
        
        row += 1
        
    worksheet2.write(row+3, 0, 'Total phosphosites detected', '')
    worksheet2.write(row+3, 1, total_phosphosites, '')
    
    worksheet2.write(row+4, 0, 'Total sites with kinase assigned', '')    
    worksheet2.write(row+4, 1, total_phosphosites_with_kinase_assigned, '')
            
    
    # # # third worksheet - contains predicted phosphorylation sites for the whole kinetochore # # # 
    
    worksheet3 = workbook.add_worksheet('All pred. kin. sites') #all kinetochore proteins
    worksheet3.freeze_panes(1, 0)
    
    #set column widths:
    worksheet3.set_column(0, 0, 15) #column A = 15
    worksheet3.set_column(1, 1, 13) #column B
    worksheet3.set_column(2, 2, 10) #column C
    worksheet3.set_column(3, 3, 12) #column D
    worksheet3.set_column(4, 4, 18)                                   
    worksheet3.set_column(5, 5, 13)
    worksheet3.set_column(6, 6, 28)
    worksheet3.set_column(7, 7, 33)
    worksheet3.set_column(8, 8, 12)
    worksheet3.set_column(9, 9, 100)
    
    #conditional formatting:
    worksheet3.conditional_format(f'E1:E1048576', { #entire column = H1:H1048576 (the last possible row is 1048576).

        'type': 'cell',
        'criteria': 'equal to',
        'value': '"NO"',
        'format': red_format
    })
    
    worksheet3.conditional_format(f'E1:E1048576', { #entire column = H1:H1048576 (the last possible row is 1048576).

        'type': 'cell',
        'criteria': 'equal to',
        'value': '"YES"',
        'format': green_format
    })
    
    worksheet3.conditional_format(f'I1:I1048576', { #entire column = H1:H1048576 (the last possible row is 1048576).

        'type': 'cell',
        'criteria': 'equal to',
        'value': '"Reported"',
        'format': blue_format
    })
    
    #write the header string to the Excel worksheet:
    
    row = 0
    col = 0
    
    output_header = [
     
        'Complex/type',
        'Accession',
        'Protein',
        'Potential phosphosite',
        'Detected?',
        'Predicted kinase',
        'Kinase consensus (regex)',
        'Sequence +/- 5',
        'Literature report? (SGD)',
        'Citation'
        
    ]
    
    header_format = workbook.add_format()
    header_format.set_align('center')
    header_format.set_align('vcenter')
    header_format.set_text_wrap()

    for header_item in output_header:
        
        worksheet3.write(row, col, header_item, header_format)
        
        col += 1 #move to next column every header item
    
    row = 1
    col = 0
    
    current_protein = '' #the name of the protein we are currently looping through. Important for the collapsible row feature in Excel.
    new_protein = False
    l_blue = True #formatting for the collapsable rows.
    
    # To read the JSON file containing all published phosphorylation sites already known:    
    json_file = open('published_phos.json', 'r')
    read = json_file.read()
    json_dict = json.loads(read) #the contents of the json file, converted to a python dictionary.
    json_file.close()
    
    for desc, complx_dict in kinetochore.proteins.items():
        for complx, proteins in complx_dict.items():
        
            for protein in proteins:                
                all_predicted_phos_sites = {} #a dictionary with site (number) as key, and another dictionary containing info as values
                all_predicted_phos_sites_sorted_dict = {} #the same dictionary as before, except sorted by residue order.
                
                accession = kinetochore.accessions[protein] #use the accession, since there are discrepancies with our common KT protein names and SGD common name
                protein = p.Protein(accession).name #now convert the protein identified from the accession back to the common name.
                
                if current_protein != protein:
                    new_protein = True #we are going to insert a collapsible row in this iteration of the loop.
                else:
                    new_protein = False
                
                current_protein = protein
                
                #loop through each kinase consensus (regex):
                for kinase, consensus in p.Protein(protein).consensus.items():


                    find_sites = p.Protein(protein).find_sites(consensus)
                    sites = find_sites.phosphoresidues

                    if len(sites) >= 1:
                        #there are matched consensus sequences!
                        
                        for site in sites:
                            
                            #put this site in the all_predicted_phos_sites dictionary
                            #remember a site can be predicted to be phos by more than one kinase!
                            
                            if site not in all_predicted_phos_sites.keys():
                                
                                all_predicted_phos_sites[site] = {

                                    'aa': p.Protein(protein).sequence[site-1],
                                    'kinase': [kinase],
                                    'consensus': [consensus]

                                }
                                
                            else:
                                
                                #the site is already in the dictionary. must match more than one kinase consensus...
                                
                                all_predicted_phos_sites[site]['kinase'].append(kinase)
                                all_predicted_phos_sites[site]['consensus'].append(consensus)
                                
                # now we have looped through all the kinases.
                # the phosphosites are going to be out of order, since we searched by kinase sequence.
                # generate a list of phosphosites, then sort them

                all_predicted_phos_sites_list = list(all_predicted_phos_sites.keys())
                all_predicted_phos_sites_list.sort()
                

                for site in all_predicted_phos_sites_list:

                    all_predicted_phos_sites_sorted_dict[site] = all_predicted_phos_sites[site] #now they are in order.
                    
                formatting = {

                    0: workbook.add_format({'align': 'center'}),
                    1: workbook.add_format({'align': 'center'}),
                    2: workbook.add_format({'align': 'center'}),
                    3: workbook.add_format({'align': 'center'}),
                    4: workbook.add_format({'align': 'center'}),
                    5: '',
                    6: workbook.add_format({'font_name': 'Courier'}),
                    7: workbook.add_format({'align': 'center', 'font_name': 'Courier'}),
                    8: workbook.add_format({'align': 'center'}),
                    9: ''
                }
                
                #create the collapsible row:
                
                worksheet3.outline_settings(True, False, True, False)
                
                if new_protein is True:
                    
                    #worksheet3.set_row(row, None, None, {'level': 1, 'collapsed': True})
                    
                    sum_of_predicted_sites = len(all_predicted_phos_sites_list)
                    
                    #etract sum of sites reported:
                    #some sites are repeatedly reported. therefore put them as dictionary keys to get the unique number:
                    unique_reported_sites = {}
                    for resideaux in all_predicted_phos_sites_list:
                        if str(resideaux) in json_dict['proteins'][accession].keys():
                            unique_reported_sites[resideaux] = 0                                           
                                           
                    sum_of_sites_reported = len(unique_reported_sites.keys())
                    
                    collapse_output_item = [
                        complx,
                        accession,
                        protein,
                        sum_of_predicted_sites,
                        '',
                        '',
                        '',
                        '',
                        sum_of_sites_reported,
                        '',
                        '',
                        '',
                        '',
                        '',
                        '',
                        '',
                        '' #added extra columns to make the blue a little longer.
                    ]
                    
                    cc = 0 #column
                    for item in collapse_output_item:
                        
                        #header rows formatting
                        
                        light_blue = workbook.add_format({'bg_color': '#E1EFFF', 'align': 'center', 'border': 1, 'border_color': '#A5C9F5'})
                        darker_blue = workbook.add_format({'bg_color': '#D5E8FF', 'align': 'center', 'border': 1, 'border_color': '#A5C9F5'})
                        
                        if l_blue is True:
                            worksheet3.write(row, cc, item, light_blue)
                        else:
                            worksheet3.write(row, cc, item, darker_blue)
                        
                        cc += 1
                        
                    if l_blue is True:
                        l_blue = False
                    else:
                        l_blue = True
                    
                    row += 1 #increase the row.
                
                #now loop through the sorted list:
                
                for site, site_info in all_predicted_phos_sites_sorted_dict.items():
                
                    # was the site detected to be phosphorylated in any of the MS files?                            
                    detected = 'NO' #default

                    if protein in all_phosphosites.keys():
                        for phosphosite in all_phosphosites[protein].keys():

                            if phosphosite == site:
                                detected = 'YES'

                    aa_site = site_info['aa'] + str(site)
                    
                    predicted_kinases = ''
                    ii = 1
                    for kinases in site_info['kinase']:
                        predicted_kinases += kinases
                        
                        if ii < len(site_info['kinase']):
                            predicted_kinases += '; '
                    
                    kinase_consensuses = ''
                    ii = 1
                    for consensus in site_info['consensus']:
                        kinase_consensuses += consensus
                        
                        if ii < len(site_info['consensus']):
                            kinase_consensuses += '; '
                        
                    #surrounding sequence:                    
                    site_minus_1 = site-1
                    five_minus = (site-1)-5
                    five_plus = (site-1)+6

                    if five_minus < 0:
                        five_minus = 0 #we don't want to go below zero, the starting methionine.

                    if five_plus > (p.Protein(protein).length-1):
                        five_plus = p.Protein(protein).length-1

                    surrounding_seq = p.Protein(protein).sequence[five_minus:site] + '*' + p.Protein(protein).sequence[site:five_plus]
                    surrounding_seq = '[' + str(five_minus+1) + ']-' + surrounding_seq + '-[' + str(five_plus) + ']'
                    
                    #literature report? :
                    
                    literature_report = ''
                    citation = ''
                    PMID = ''
                    
                    for acc, reported_phos_info in json_dict['proteins'].items():
                        
                        if accession == acc:
                            
                            #check whether this site has been reported before:
                            if str(site) in reported_phos_info.keys():
                                
                                literature_report = 'Reported'
                                
                                ii = 0                                
                                for ref in reported_phos_info[str(site)]['refs']:
                                    PMID = reported_phos_info[str(site)]['ref_PMIDs'][ii]
                                    PMID = str(PMID)
                                    
                                    citation += ref + ' (PMID: ' + PMID + ')'
                                    
                                    if ii < len(reported_phos_info[str(site)]['refs'])-1:
                                        citation += '; '
                                    ii+=1
                    
                    output_items = [

                        '',
                        '',
                        protein,
                        aa_site,
                        detected,
                        predicted_kinases,
                        kinase_consensuses,
                        surrounding_seq,
                        literature_report,
                        citation

                    ]
                    
                    
                    for output_item in output_items:                        
                    
                        worksheet3.set_row(row, None, None, {'level': 1, 'hidden': True})
                        worksheet3.write(row, col, output_item, formatting[col])
                        
                        col += 1

                    row += 1
                    col = 0 #reset back to 0 for each site.
                    
                    worksheet3.write(row, col, 'Note that this tab only contains sites that have been predicted to be phosphorylated by a kinase provided. This is not a full list of all reported phosphorylation.')
                    
    # # # fourth worksheet (the regex info):
    
    worksheet4 = workbook.add_worksheet('Regex used') #all kinetochore proteins
    
    formatting = {
    
        0: workbook.add_format({'align': 'center'}),
        1: workbook.add_format({'align': 'left', 'font_name': 'Courier'})
        
    }
    
    #set column widths:
    worksheet4.set_column(0, 0, 20) #column A = 15
    worksheet4.set_column(1, 1, 60) #column B
    
    worksheet4.write(0, 0, 'Kinase')
    worksheet4.write(0, 1, 'Regular expression used')
    
    row = 1
    
    for kinase, reg_exp in p.Protein('NDC80').consensus.items(): #it just needs a protein name in order to work... I just put Ndc80.
     
        worksheet4.write(row, 0, kinase, formatting[0])
        worksheet4.write(row, 1, reg_exp, formatting[1])
        
        row += 1
        
    
    workbook.close()
    
    path = current_output_path + '/' + output_file_names('compare_all_phos__' + num_of_files + '_cond_', files, cond)
    
    if os.path.isfile(path):
        
        return True
    
    else:
        
        return False