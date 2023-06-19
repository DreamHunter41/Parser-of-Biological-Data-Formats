# By FreshGirls
# 制表符tab \t
class VCFHeader():
    usual_header = ["INFO","FILTER","FORMAT"]
    def __init__(self,line):
        self.line = line
        # judge whether the line is a header
        if line[0] != "#" or line[1] !=line[0]:
            print("Input is not a header, please check the format.")
        # judge whether the header is in usual format, meta-information or others
    def ana_headers(self):
        item = 0
        self.item_name = ""
        self.item_cont1 = ""
        item_cont_ad = 0
        line = self.line
        for i in range(2,len(line)):
            if item == 0 and line[i]!= "=":
                self.item_name += line[i]
            elif line[i] == "=":
                self.name = self.item_name
                item = 1
            elif item == 1 and line[i] != "=" and line[i] !="<" and item_cont_ad == 0:
                self.item_cont1 += line[i]
            elif item == 1 and line[i] == "<" and item_cont_ad == 0:
                self.item_cont2 = MetaInfo()
                item_cont_ad = 1
            elif item == 1 and line[i] != "=" and line[i] != "<" and item_cont_ad == 1:
                self.item_cont2.save(line[i])
            elif item == 1 and line[i] == "=" and line[i] !="<" and item_cont_ad == 1:
                self.item_cont2.mode = 1
                self.item_cont2.name.append(self.item_cont2.current_name)
            elif item == 1 and line[i] == "," and item_cont_ad ==1:
                self.item_cont2.mode = 0
                self.item_cont2.cont.append(self.item_cont2.current_cont)
        if item_cont_ad == 0:
            return (self.item_name,self.item_cont1)
        else:
            return (self.item_name,self.item_cont2)
        
    def classificate(self,metainfo):
        '''The intergration of the most common comments.'''
        self.info_lst = []
        self.filter_lst= []
        self.format_lst = []
        for i in range(len(metainfo.name)):
            if metainfo.name[i] == "INFO":
                self.info_lst.append(metainfo.cont[i])
            elif metainfo.name[i] == "FILTER":
                self.filter_lst.append(metainfo.cont[i])
            elif metainfo.name[i] == "FORMAT":
                self.format_lst.append(metainfo.cont[i])
        return (self.info_lst,self.filter_lst,self.format_lst)

class MetaInfo(VCFHeader):
    def __init__(self):
        self.name = []
        self.current_name = ""
        self.cont = []
        self.current_cont= ""
        self.mode = 0 # 0 when saving name and 1 when saving the content
    def save(self,char):
        if self.mode == 0:
            self.current_name += char
        if self.mode == 1:
            self.current_cont += char
        
class VCF():
    def __init__(self,file):
        self.file = open(file,"r")
    def parse_header(self):
        # create an empty dictionary for headers
        self.headers = dict()
        line = self.file.readline()
        while line[0] == line[1] and line[0] == "#":
            header = VCFHeader(line).ana_headers()
            self.headers[header[0]]=header[1]
            line = self.file.readline() 
        return self.headers
    def parse_main(self):
        '''Suppose you have finished parsing the headers. 
        You are now arriving at the first line of main content which starts with a line led by "#".
        '''
        # parse the header for main content
        line = self.file.readline()
        self.main_items = []
        if line[0] == "#":
            temp = ""
            for i in range(1,len(line)):
                if line[i] != "\t":
                    temp += line[i]
                else:
                    self.main_items.append(temp)
                    temp = ""
        # parse the main content
        line = self.file.readline()
        temp2 = ""
        self.main_content=[]
        self.all_main_content=[]
        while line:
            for i in range(0,len(line)):
                if line[i] != "\t":
                    temp2 += line[i]
                else:
                    self.main_content.append(temp2)
                    temp2 = ""
            self.all_main_content.append(self.main_content)
            line = self.file.readline()
        return self.all_main_content

# test
test = VCF("sample.vcf")
test.parse_header()
test.parse_main()