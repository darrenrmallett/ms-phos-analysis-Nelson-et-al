import sys, os, pandas, proteins, mass_spec

current_filepath = os.path.dirname(os.path.realpath(__file__)) + '/msfiles'

class FileReader:
    
    def __init__(self, filename, filepath = current_filepath):
        
        self.Error = None
        self.filename = filename
        self.filepath = filepath
        self.current_filepath = current_filepath
        
        self.accepted_file_exts = ["xlsx", "xls", "xlsb", "xlsm"]
        
        #get the file extension for the file:
            
        self.file_ext = self.filename.split(".")
        self.file_ext = self.file_ext[-1]

        self.filename_without_extension = filename.split(self.file_ext)[0]
        self.file_ext = self.file_ext.lower()

        self.convert_file_to_csv(filename)

    def convert_file_to_csv(self, file):
        
        #filehandle where the new CSV file will be stored:        
        #check to make sure the file exists first:
    
        if file in os.listdir(self.current_filepath):
        
            self.filename = file
            self.csv_filename = self.filename_without_extension + 'csv'
            self.filehandle = open(self.filepath + '/' + self.csv_filename, 'w')
            
            #if the file is a MS Excel file, then convert it to CSV file:        
            if self.file_ext in self.accepted_file_exts:

                #use Pandas library to read Excel file and convert to CSV.
                data_frame = pandas.read_excel(self.filepath + '/' + self.filename)
                data_frame.to_csv(self.filehandle, sep='\t', index=False, index_label=False, lineterminator='\n')
                
                self.Error = None
            
            else: 
                
                self.Error = f"invalid file type for the following file: {file}. Please use MS Excel file (.xlsx or .xls)"

            self.filehandle.close()
            
        else:
            
            self.Error = f"cannot find the following file: {file}. Are you sure this file exists? Please try again."