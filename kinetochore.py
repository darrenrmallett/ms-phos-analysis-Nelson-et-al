#kinetochore.py - Kinetochore() class - contains all proteins in the kinetochore.

import sys, os, proteins as p, re, requests, xlsxwriter, time, json

current_output_path = os.path.dirname(os.path.realpath(__file__)) + '/output'
current_filepath = os.path.dirname(os.path.realpath(__file__))
unix = int(time.time()) #UNIX timestamp.

class Kinetochore:

    def __init__(self):
        
        #generate data structures that store kinetochore protein info:
        
        self.proteins = {
        
            'Inner Kinetochore': {
                
                'CBF3 complex': ['Ndc10', 'Cep3', 'Ctf19', 'Skp1'],
                'Nucleosome': ['Cse4', 'HHF1', 'HHF2', 'HTA1', 'HTA2', 'HTB1', 'HTB2'],
                'CCAN' : ['Mcm16', 'Ctf3', 'Mcm22', 'Iml3', 'Chl4', 'Nkp1', 'Nkp2', 'Mif2', 'Ctf19', 'Mcm21', 'Cnn1', 'Wip1'],
                'OA': ['Okp1', 'Ame1']
                
            },
            
            'Outer kinetochore': {

                'Mis12 complex': ['Mtw1', 'Nnf1', 'Nsl1', 'Dsn1'],
                'Ndc80 complex': ['Ndc80', 'Nuf2', 'Spc24', 'Spc25'],
                'Spc105 complex': ['Spc105', 'Kre28'],
                'Dam1 complex': ['Dam1', 'Dad1', 'Dad2', 'Dad3', 'Dad4', 'Ask1', 'Duo1', 'Hsk3', 'Spc19', 'Spc34'],
                'SAC': ['Bub3', 'Bub1', 'Mad1', 'Mad2', 'Mps1', 'Glc7']
            },
            
            'MAPs': {
            
                'MAPs': ['Stu2', 'Stu1', 'Slk19'],
                'Motors': ['Kar3', 'Cin8']
                
            },
            
            'Other': {
            
                'CPC': ['Ipl1', 'Bir1', 'Nbl1', 'Sli15'],
                'regulators/other': ['Psh1', 'Smt3', 'Cdc5', 'Sgo1', 'Cdc14', 'Cdk1', 'Cdc55', 'Ubr2', 'Fin1']
                
            }
            
        }
        
        
        self.accessions = {  'Ame1': 'YBR211C',
                             'Ask1': 'YKL052C',
                             'Bir1': 'YJR089W',
                             'Bub1': 'YGR188C',
                             'Bub3': 'YOR026W',
                             'Cdc14': 'YFR028C',
                             'Cdc5': 'YMR001C',
                             'Cdc55': 'YGL190C',
                             'Cdk1': 'YBR160W',
                             'Cep3': 'YMR168C',
                             'Chl4': 'YDR254W',
                             'Cin8': 'YEL061C',
                             'Cnn1': 'YFR046C',
                             'Cse4': 'YKL049C',
                             'Ctf19': 'YPL018W',
                             'Ctf3': 'YLR381W',
                             'Dad1': 'YDR016C',
                             'Dad2': 'YKR083C',
                             'Dad3': 'YBR233W-A',
                             'Dad4': 'YDR320C-A',
                             'Dam1': 'YGR113W',
                             'Dsn1': 'YIR010W',
                             'Duo1': 'YGL061C',
                             'Glc7': 'YER133W',
                             'HHF1': 'YBR009C',
                             'HHF2': 'YNL030W',
                             'HTA1': 'YDR225W',
                             'HTA2': 'YBL003C',
                             'HTB1': 'YDR224C',
                             'HTB2': 'YBL002W',
                             'Hsk3': 'YKL138C-A',
                             'Iml3': 'YBR107C',
                             'Ipl1': 'YPL209C',
                             'Kar3': 'YPR141C',
                             'Kre28': 'YDR532C',
                             'Mad1': 'YGL086W',
                             'Mad2': 'YJL030W',
                             'Mcm16': 'YPR046W',
                             'Mcm21': 'YDR318W',
                             'Mcm22': 'YJR135C',
                             'Mif2': 'YKL089W',
                             'Mps1': 'YDL028C',
                             'Mtw1': 'YAL034W-A',
                             'Nbl1': 'YHR199C-A',
                             'Ndc10': 'YGR140W',
                             'Ndc80': 'YIL144W',
                             'Nkp1': 'YDR383C',
                             'Nkp2': 'YLR315W',
                             'Nnf1': 'YJR112W',
                             'Nsl1': 'YPL233W',
                             'Nuf2': 'YOL069W',
                             'Okp1': 'YGR179C',
                             'Psh1': 'YOL054W',
                             'Sgo1': 'YOR073W',
                             'Skp1': 'YDR328C',
                             'Sli15': 'YBR156C',
                             'Slk19': 'YOR195W',
                             'Smt3': 'YDR510W',
                             'Spc105': 'YGL093W',
                             'Spc19': 'YDR201W',
                             'Spc24': 'YMR117C',
                             'Spc25': 'YER018C',
                             'Spc34': 'YKR037C',
                             'Stu1': 'YBL034C',
                             'Stu2': 'YLR045C',
                             'Wip1': 'YDR374W-A',
                             'Ubr2': 'YLR024C',
                             'Fin1': 'YDR130C'
                          }