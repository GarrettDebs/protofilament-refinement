from info_file import InfoFile
import pandas as pd
import numpy as np

class StarFile(InfoFile):
    def __init__(self, file):
        self.readStar(file)
        
    def readStar(self, file):
        f=open(file)
        lines=f.readlines()
        f.close()
        header, start=self.readHeader(lines)
        body=self.readBody(start, lines)
        ###Convert to a DataFrame
        self.df=pd.DataFrame(body, columns=header)
        
    def readHeader(self, lines):
        self.top=[]
        header=[]
        cont=True
        for i in range(len(lines)):
            ###Make sure to ignore the non-column names, but retain them
            ###for printing leader
            if cont and not lines[i].startswith('_rln'):
                self.top.append(lines[i])
                continue
            ##Keep the actual column names and remove the prefix
            elif lines[i].startswith('_rln'):
                header.append(lines[i].split()[0][4::])
                cont=False
            else:
                break
        
        ###Return the header and the line that the body starts at    
        return header, i
    
    def readBody(self, start, lines):
        body=[]
        
        #while lines[-1].split()==[]:
        #    print 'ha'
        #    lines.pop(-1)
        
        for i in range(start, len(lines)):
            temp=lines[i].split()
            ###Make sure that we aren't adding blank lines to the body
            if temp:
                body.append(temp)
            
        return np.array(body)
    
    def getHeader(self):
        header=self.top
        i=1
        for col in self.df.columns:
            header.append('_rln%s #%g\n'%(col, i))
            i+=1
            
        return header
    
    def writeStar(self, output, data=None):
        if data is None:
            data=self.df.to_numpy().tolist()
            
        n=open(output,'w')
        header=self.getHeader()
        n.write(''.join(header))
        
        new=[]
        if data[0][-1] is not '\n':
            for i in range(len(data)):
                data[i].append('\n')
                new.append(' '.join(data[i]))
        else:
            for i in range(len(data)):
                new.append(' '.join(data[i]))
            
        n.write(''.join(new))
        n.close()