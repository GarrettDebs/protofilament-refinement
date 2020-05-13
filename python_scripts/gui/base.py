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
        
        self.file_var={}
        self.file_button={}
        
        self.entry_label={}
        self.entry_vars={}        
        
        for key, item in kwargs.items():
            if type(item) is bool:
                ###Initialize the checkbox
                self.bool_vars[key]=tk.BooleanVar()
                self.bool_vars[key].set(item)
                
                self.bool_box[key]=tk.Checkbutton(self.top, text=key, \
                                             variable=self.bool_vars[key])
                
            elif item.startswith('XXX'):
                self.addFile(key, item)
            
            else:
                pretty_var=' '.join(key.split('_')).title()
                self.entry_label[key]=tk.Label(self.top, text=pretty_var)
                self.entry_vars[key]=tk.Entry(self.top)
                self.entry_vars[key].insert(0, item)
    
    def addMrcFile(self, key):            
        self.top.filename = tkFileDialog.askopenfilename\
            (initialdir = ".",title = "Select file",filetypes = \
             (("MRC files","*.mrc"),("all files","*.*")))
        
        self.file_var[key].delete(0,'end')
        self.file_var[key].insert(0, self.top.filename)
    
    def addStarFile(self, key):            
        self.top.filename = tkFileDialog.askopenfilename\
            (initialdir = ".",title = "Select file",filetypes = \
            (("star files","*.star"),("all files","*.*")))
            
        self.file_var[key].delete(0,'end')
        self.file_var[key].insert(0, self.top.filename)
        
    def addAllFile(self, key):            
        self.top.filename = tkFileDialog.askopenfilename\
            (initialdir = ".",title = "Select file")
            
        self.file_var[key].delete(0,'end')
        self.file_var[key].insert(0, self.top.filename)
        
    def addFile(self, key, item):
        self.file_var[key]=tk.Entry(self.top)
        pretty_var=' '.join(key.split('_')).title()
        self.file_var[key].insert(0, pretty_var)
        
        if 'MRC' in item:
            self.file_button[key]=tk.Button(self.top, text='Browse Files', \
                                        command=lambda var=key: self.addMrcFile(var))

        elif 'STAR' in item:
            self.file_button[key]=tk.Button(self.top, text='Browse Files', \
                                        command=lambda var=key: self.addStarFile(var))
            
        else:
            self.file_button[key]=tk.Button(self.top, text='Browse Files', \
                                        command=lambda var=key: self.addAllFile(var))
        
                    
    def packVars(self):
        num_bools=len(self.bool_vars)
        num_entry=len(self.entry_vars)
        
        total_vars=num_bools+num_entry
        
        logging.info('For now only dealing with two columns')
        self.i=0
        
        for key in self.file_var:
            self.file_var[key].grid(row=self.i, column=0)
            self.file_button[key].grid(row=self.i, column=1)
            self.i+=1             
        
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
        
        for key in self.file_var:
            self.vals[key]=self.file_var[key].get()
        
        for key in self.bool_vars:
            self.vals[key]=self.bool_vars[key].get()
        
        for key in self.entry_vars:
            self.vals[key]=self.entry_vars[key].get()
            
        self.top.destroy()
            
    def sendValues(self):
        return self.vals  
             
