import Tkinter as tk
import tkFileDialog
import logging
#from Tkinter.tkk import Notebook

logging.basicConfig(level=logging.DEBUG)

class gdGui():
    def __init__(self, title, **kwargs):
        ###Title of the gui
        self.top=tk.Tk()
        self.top.title(title)
        
        ##Add variables
        self.addVars(**kwargs)
        
        ##Determine layout
        self.packVars()
        
        ##Add button to return parameters
        self.addButton()
        
        self.top.mainloop()
                
    def addVars(self, **kwargs):
        self.bool_vars={}
        self.bool_box={}
        
        self.entry_label={}
        self.entry_vars={}        
        
        for key, item in kwargs.items():
            if type(item) is bool:
                ###Initialize the checkbox
                self.bool_vars[key]=tk.BooleanVar()
                self.bool_vars[key].set(item)
                
                self.bool_box[key]=tk.Checkbutton(self.top, text=key, \
                                             variable=self.bool_vars[key])
                
            elif item is 'XXXFILENAMEXXX':
                self.file_var=tk.Entry(self.top)
                self.file_button=tk.Button(self.top, text='Browse Files', \
                                           command=self.addFile)
            
            else:
                pretty_var=' '.join(key.split('_')).title()
                self.entry_label[key]=tk.Label(self.top, text=pretty_var)
                self.entry_vars[key]=tk.Entry(self.top)
                self.entry_vars[key].insert(0, item)
                
    def addFile(self):
        self.top.filename = tkFileDialog.askopenfilename\
        (initialdir = ".",title = "Select file",filetypes = \
         (("star files","*.star"),("all files","*.*")))
        self.file_var.insert(0, self.top.filename)
    
    def packVars(self):
        num_bools=len(self.bool_vars)
        num_entry=len(self.entry_vars)
        
        total_vars=num_bools+num_entry
        
        logging.info('For now only dealing with two columns')
        self.i=0
        
        try:
            self.file_var.grid(row=0, column=0)
            self.file_button.grid(row=0, column=1)
            self.i+=1
        except:
            pass                
        
        for key in self.bool_vars:
            self.bool_box[key].grid(row=self.i, column=0)
            self.i+=1
            
        for key in self.entry_vars:
            self.entry_label[key].grid(row=self.i, column=0, padx=5, pady=5)
            self.entry_vars[key].grid(row=self.i, column=1, padx=5, pady=5)
            self.i+=1
            
    def addButton(self):
        button=tk.Button(self.top, text='Run', command=self.returnValues).grid(row=self.i, \
                                                                column=1)
        
    def returnValues(self):
        self.vals={}
        
        ###Check if there is a file name variable
        try:
            self.vals['filename']=self.top.filename
        except:
            pass
        
        for key in self.bool_vars:
            self.vals[key]=self.bool_vars[key].get()
        
        for key in self.entry_vars:
            self.vals[key]=self.entry_vars[key].get()
            
        self.top.destroy()
            
    def sendValues(self):
        return self.vals  
             
