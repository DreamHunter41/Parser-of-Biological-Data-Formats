# By FreshGirls
class FastqReader():
    def __init__(self,file):
        self.file = open(file,'r')
        self.id = ""
    
    def split(self):
        '''Split the file into several records and save them as a list of FastqRecord.'''
        # initiating the reading process
        self.line = self.file.readline()
        if self.line[0] == "@":
            header = self.line
        else:
            print("Please check the format of the file.")
        self.record=[]
        record_1 = FastqRecord()
        record_1.header=header
        record_1.sequence=self.file.readline()
        self.file.readline() #for "+"
        record_1.quality_c = self.file.readline()
        self.record.append(record_1)
        nextline = self.file.readline()
        while nextline:
            line = nextline
            if line[0]=="@":
                record_temp = FastqRecord()
                record_temp.header = line
                record_temp.sequence = self.file.readline()
                self.file.readline()
                record_temp.quality_c = self.file.readline()
                self.record.append(record_temp)
            nextline = self.file.readline()
        return self.record
        # May need to use the function len(self.record) to see how many records are there in the file.

class FastqRecord():
    '''A fastq file may include the record of several different runs.
    We therefore need to split the target file into several single records.
    '''
    def __init__(self):
        self.header = ""
        self.sequence = ""
        self.quality_c = ""
        self.quality_P = []
    
    def parse_header(self):
        '''Return the specific information stored in the header, in the form of a dictionary.'''
        self.header_items = []
        item = ""
        classical = 0 # To judge whether a header is in the classical format
        equal = 0
        for i in range(1,len(self.header)):
            if self.header[i] != ':' and self.header[i] != " ":
                item += self.header[i]
                classical = 0
            elif self.header[i] == ":":
                classical = 1
                self.header_items.append(item)
                item = ""
            elif self.header[i] == " ":
                classical=0
                self.header_items.append(item)
                if "=" in item:
                    equal = 1
                item = ""

        # for some fastaq files, the classical header format is adopted
        # however in some other fastq files, users use equations to represent the information     
        if classical == 1:
            self.header_c=['id','ri','fci','fcl','fcltn','x','y','pair','filtered','cb','Barcode']
            self.header_parsed = dict()
            for i in range(len(self.header_items)):
                key = self.header_c[i]
                value = self.header_items[i]
                self.header_parsed[key]=[value]
        elif equal == 1:
            self.header_parsed = dict()
            #substract all the features
            for i in range(1,len(self.header_items)):
                key = ""
                value = ""
                ctrl = 0
                for j in range(len(self.header_items[i])):
                    chrct = self.header_items[i][j]
                    if chrct != "=" and ctrl == 0:
                        key += chrct
                    elif chrct != "=" and ctrl == 1:
                        value +=chrct
                    else:
                        ctrl = 1
                self.header_parsed[key]=[value]
            self.header_parsed["id"]=[self.header_items[0]]
        else:
            self.header_parsed["id"]=[self.header_items[0]]
        self.id = self.header_parsed["id"] #for a quick look-up
        return self.header_parsed

    def trans_quality(self,quality_chr):
        '''
        A function for calculating the quality from single characters.
        Only for +33 PHRED score (as Illumina uses)! Might need to be changed in different situations.
        '''
        quality_asc = ord(quality_chr)
        # .
        quality_PHRED = quality_asc-33
        return quality_PHRED
    
    def calculate_quality(self):
        '''Calculate and store the quality of each base.'''
        for i in range(len(self.quality_c)):
            rst = self.trans_quality(self.quality_c[i])
            self.quality_P.append(rst)
        return self.quality_P
    
    def avrg_quality(self):
        '''The overall quality of the run.'''
        ttl = 0
        length = len(self.quality_P)
        for i in range(length):
            ttl += self.quality_P[i]
        self.avrg = ttl/length
        return self.avrg
    
    def identifier(self,location):
        '''Generates base-quality pairs for further analyses.'''
        index = int(location)-1
        base = self.sequence[index]
        qlty = self.quality_P[index]
        return (base,qlty)


# Test
test = FastqReader('sample.fastq')
spt = test.split()
record_num = len(spt)
demo_header = spt[0].parse_header()
print(demo_header['runid'])
spt[0].calculate_quality()
demo_avrg_quality = spt[0].avrg_quality()
print(demo_avrg_quality)
demo_identifier = spt[0].identifier(4)
print(demo_identifier)
