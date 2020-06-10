from pkg_resources import resource_stream
from PIL import ImageTk, Image
import Tkinter as tk
import tkFileDialog
import os

class gdGui():
    def __init__(self, title, **kwargs):
        ###Title of the gui
        self.top=tk.Tk()
        self.top.title(title)
        
        ##Add variables
        self.addVars(**kwargs)
        
        ##Determine layout
        self.addImage()
        self.packVars()
        
        ##Add button to return parameters
        self.addButton()
        
        self.top.mainloop()
                
    def addVars(self, **kwargs):
        self.bool_vars={}
        self.bool_box={}
        self.bool_label={}
        
        self.file_var={}
        self.file_button={}
        
        self.entry_label={}
        self.entry_vars={}        
        
        for key, item in kwargs.items():
            if type(item) is bool:
                ###Initialize the variable. Set to string variable for display
                self.bool_vars[key]=tk.StringVar()
                self.bool_vars[key].set(str(item))
                
                pretty_var=' '.join(key.split('_')).title()
                self.bool_label[key]=tk.Label(self.top, text=pretty_var)
                self.bool_box[key]=tk.OptionMenu(self.top, self.bool_vars[key],
                                                 *['True', 'False'])
                
            elif type(item) is str and item.startswith('XXX'):
                self.addFile(key, item)
            
            else:
                pretty_var=' '.join(key.split('_')).title()
                self.entry_label[key]=tk.Label(self.top, text=pretty_var)
                self.entry_vars[key]=tk.Entry(self.top)
                self.entry_vars[key].insert(0, item)
    
    def prettyFile(self):
        try:
            self.cwd
        except Exception:
            self.cwd=os.getcwd()
            
        return os.path.relpath(self.top.filename, self.cwd)
    
    def addMrcFile(self, key):            
        self.top.filename = tkFileDialog.askopenfilename\
            (initialdir = ".",title = "Select file",filetypes = \
             (("MRC files","*.mrc"),("all files","*.*")))
        
        filename=self.prettyFile()
        self.file_var[key].delete(0,'end')
        self.file_var[key].insert(0, filename)
    
    def addStarFile(self, key):            
        self.top.filename = tkFileDialog.askopenfilename\
            (initialdir = ".",title = "Select file",filetypes = \
            (("star files","*.star"),("all files","*.*")))
            
        filename=self.prettyFile()
        self.file_var[key].delete(0,'end')
        self.file_var[key].insert(0, filename)
        
    def addAllFile(self, key):            
        self.top.filename = tkFileDialog.askopenfilename\
            (initialdir = ".",title = "Select file")
            
        filename=self.prettyFile()
        self.file_var[key].delete(0,'end')
        self.file_var[key].insert(0, filename)
        
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
        
                    
    def addImage(self):
        imname=resource_stream('pf_refinement','data/mt50.jpg')
        self.img=ImageTk.PhotoImage(Image.open(imname))
        self.imlabel=tk.Label(self.top, image=self.img)
        self.imlabel.grid(row=0, column=0, columnspan=1, rowspan=30, padx=5)
        
    def packVars(self):
        num_bools=len(self.bool_vars)
        num_entry=len(self.entry_vars)
        
        total_vars=num_bools+num_entry
        
        self.i=0
        
        for key in self.file_var:
            self.file_var[key].grid(row=self.i, column=1)
            self.file_button[key].grid(row=self.i, column=2)
            self.i+=1             
        
        for key in self.bool_vars:
            self.bool_label[key].grid(row=self.i, column=1)
            self.bool_box[key].grid(row=self.i, column=2)
            self.i+=1
            
        for key in self.entry_vars:
            self.entry_label[key].grid(row=self.i, column=1, padx=5, pady=5)
            self.entry_vars[key].grid(row=self.i, column=2, padx=5, pady=5)
            self.i+=1
            
    def addButton(self):
        button=tk.Button(self.top, text='Run', command=self.returnValues)
        button.grid(row=self.i, column=2)
        
    def returnValues(self):
        self.vals={}
        
        for key in self.file_var:
            self.vals[key]=self.file_var[key].get()
        
        for key in self.bool_vars:
            self.vals[key]=self.bool_vars[key].get()=='True'
        
        for key in self.entry_vars:
            self.vals[key]=self.entry_vars[key].get()
            
        self.top.destroy()
            
    def sendValues(self):
        try:
            return self.vals
        except Exception:
            raise AttributeError('No variables were not defined in the GUI')  
             
