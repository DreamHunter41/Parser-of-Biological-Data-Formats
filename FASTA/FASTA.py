# By FreshGirls
class FastaReader():
    def __init__(self,file):
        """Initializing the class, importing the FASTA file."""
        self.file = open(file,'r') #Import the file that we want to read
        self.line = "" #Prepare to start reading lines
    def read_header(self):
        self.header = self.file.readline()
        if self.header[0] == '>':
            return self.header
        else:
            print('Please check the format of the FASTA file.')
            return None
    def read_sequence(self):
        self.sequence = ""
        while self.file.readline():
            line_on_hold = self.file.readline().strip()
            if line_on_hold[0] != ";":
                self.sequence += line_on_hold
        return self.sequence
'''improvement方向:
1. 对header进行拆解,提取header中的信息并作为features矩阵呈现(比如将accession number单拉出来 CDS = 0 or 1, mRNA = 0 or 1);
    【但是因为FASTA的header多样性实在太高,比较难拆,需要较多判断】
2. 对sequence进行初步分析,概括如sequence的长度、类型等,作为properties呈现
'''