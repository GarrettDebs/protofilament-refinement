from info_file import InfoFile
import pandas as pd
import numpy as np

class StarFile(InfoFile):
    def __init__(self, file=None):
        try:
            self.readStar(file)
        except Exception:
            pass
        
        self.name=file
        
    def readStar(self, file):
        ### Read the enitre file
        f=open(file)
        lines=f.readlines()
        f.close()
        
        ### Read the header and then the body
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
        
        for i in range(start, len(lines)):
            temp=lines[i].split()
            ###Make sure that we aren't adding blank lines to the body
            if temp:
                body.append(temp)
            
        return np.array(body)
    
    def getHeader(self):
        header=[]
        
        ### Add all the weird stuff to the top of the header
        header.extend(self.top)
        i=1
        
        ### Add the _rln format to the column names along with column number
        for col in self.df.columns:
            header.append('_rln%s #%g\n'%(col, i))
            i+=1
            
        return header
    
    def writeStar(self, output, data=None):
        ### If no data is provided, use the data from the DataFrame
        if data is None:
            data=self.df.to_numpy().tolist()
            
        ### Check if input is a DataFrame, otherwise assume it is a list
        try:
            data=data.to_numpy().tolist()
        except Exception:
            pass
            
        n=open(output,'w')
        
        ### First write out the header
        header=self.getHeader()
        n.write(''.join(header))
        
        new=[]
        ### Next write out the body
        if data[0][-1] is not '\n':
            for i in range(len(data)):
                data[i].append('\n')
                new.append(' '.join(data[i]))
        else:
            for i in range(len(data)):
                new.append(' '.join(data[i]))
            
        n.write(''.join(new))
        n.close()