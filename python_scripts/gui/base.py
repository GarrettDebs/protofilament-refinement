import Tkinter as tk
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
        self.entry_label={}
        self.entry_vars={}        
        
        for key, item in kwargs.items():
            if type(item) is bool:
                #bool_vars[key]=tk.BooleanVar(default=item)
                logging.warning('Have not programmed in Boolean Variables '\
                             'for the gui yet!')
                exit()
                
            else:
                pretty_var=' '.join(key.split('_')).title()
                self.entry_label[key]=tk.Label(self.top, text=pretty_var)
                self.entry_vars[key]=tk.Entry(self.top)
                self.entry_vars[key].insert(0, item)
                
    def packVars(self):
        num_bools=len(self.bool_vars)
        num_entry=len(self.entry_vars)
        
        total_vars=num_bools+num_entry
        
        logging.info('For now only dealing with two columns')
        self.i=0
        for key in self.bool_vars:
            logging.warning('No boolean variables in packVars yet!')
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
        for key in self.bool_vars:
            logging.warning('No booleans in returnValues yet!')
        
        for key in self.entry_vars:
            self.vals[key]=self.entry_vars[key].get()
            
        self.top.destroy()
            
    def sendValues(self):
        return self.vals  
             
