###For Mac users
from sys import platform as sys_pf

if sys_pf == 'darwin':
    print('MacOS detected.')
    import matplotlib
    matplotlib.use("TkAgg")

from scipy.optimize import least_squares
from pf_refinement import StarOp
import matplotlib.pyplot as plt
from scipy import stats
from gui import gdGui
import pandas as pd
import numpy as np
import os

class PlotDeltaPhi(StarOp):
        def __init__(self, star_file=None, gui=False):
            if gui:
                init_vals=self.initVals()
                temp=gdGui('Plot Distortions', **init_vals)
                self.vals=temp.sendValues()
                super(PlotDeltaPhi, self).__init__(self.vals['input_star'])
            else:
                super(PlotDeltaPhi, self).__init__(star_file)
                
            self.delPhiXs()
            self.setParams()
            self.dirname='plot_dphis'
            
        def initVals(self):
            vals={
                'input_star': 'XXXSTARFILEXXX',
                'plot_indiv_mt': False
                }
            return vals
            
        def setParams(self):
            ### Plotting Parameters, feel free to adjust
            self.bin_size=0.7
            self.xlim=[10,50]
            self.x_spacing=10
            self.bin_use=np.linspace(self.xlim[0], self.xlim[-1], 
                                  (self.xlim[-1]-self.xlim[0])/self.bin_size+1)
            self.im_size=(5,4)
            self.lw=2
            self.dphie=360.0/self.num_pfs
            
            self.min_peak_separation=2
            
        def __call__(self):
            self.plotAll()
            self.plotSeam()
            self.plotNonseam()
            if self.vals['plot_indiv_mt']:
                self.plotIndiv()

        def plotStuff(self, data, output_root):
            plt.figure(figsize=self.im_size)
            ax=plt.axes()
            plt.ylabel('Frequency')
            plt.xlabel(r'$\Delta\phi$ (Degrees)')
            
            plt.xlim(self.xlim)
            xlabels=[]
            xrange=range(self.xlim[0], self.xlim[1]+1, self.x_spacing)
            xrange=map(str, xrange)
            xrange.append('%0.1f'%self.dphie)
            xrange.sort()
            for i in xrange:
                xlabels.append(r'$%s\degree$'%i)
            
            ax.set_xticks(map(float,xrange))
            ax.set_xticklabels(xlabels)
            
            plt.axvline(self.dphie, 0, 0.1, linestyle='-', 
                        color='black', linewidth=1.5)
            plt.hist(data, bins=self.bin_use, density=True, color='dimgrey')
            
            plt.savefig(output_root+'.jpg')
            plt.savefig(output_root+'.pdf')
            
            self.best_fit(data)
            plt.savefig(output_root+'_bestfit.jpg')
            plt.savefig(output_root+'_bestfit.pdf')
            plt.close()          
            
        def best_fit(self, data):
            x=np.linspace(self.xlim[0], self.xlim[-1], 500)
            kde=stats.gaussian_kde(data)
            y=kde.evaluate(x)
            
            y_double, y1, y2, mu1, mu2 =self.fit_double_gauss(x,y)

            
            if np.abs(mu2-mu1)>self.min_peak_separation:
                ys=np.array([y_double, y1, y2]).transpose()
                plt.plot(x, ys)
                
                title=r'Peak1=$%0.1f\degree$, Peak2=$%0.1f\degree$'\
                %(mu1, mu2)
                plt.title(title)

            else:
                y_single, mu1_single =self.fit_gauss(x,y)
                plt.plot(x, y_single)
                
                title=r'Peak1=$%0.1f\degree$'%(mu1_single)
                plt.title(title)                
                
        def fit_double_gauss(self, x, y):
            def double_gaussian(x, params):
                (c1, mu1, sigma1, c2, mu2, sigma2) = params
                res=c1*np.exp(-(x-mu1)**2.0/(2.0*sigma1**2.0)) + \
                c2*np.exp(-(x-mu2)**2.0/(2.0*sigma2**2.0))
                return res
    
            def double_gaussian_fit(params):
                fit=double_gaussian(x,params)
                return (fit-y)

            def gaussian(x, params):
                (c1, mu1, sigma1) = params
                res=c1*np.exp(-(x-mu1)**2.0/(2.0*sigma1**2.0))
                return res
            
            bounds=([-100, -100, -100, -100, -100, -100], 
                    [100, 100, 100, 100, 100, 100])
            fit=least_squares(double_gaussian_fit,
                              [1.0, self.dphie-1.0, 1.0, 
                               0.5, self.dphie+3, 2.0], 
                              bounds=bounds)
            
            mu1=fit['x'][1]
            mu2=fit['x'][-2]
            
            y_fit=double_gaussian(x, fit['x'])
            y1_fit=gaussian(x, fit['x'][:-3])
            y2_fit=gaussian(x, fit['x'][3:])
            
            return y_fit, y1_fit, y2_fit, mu1, mu2
            
        def fit_gauss(self, x, y):
            def gaussian(x, params):
                (c1, mu1, sigma1) = params
                res=c1*np.exp(-(x-mu1)**2.0/(2.0*sigma1**2.0))
                return res
    
            def gaussian_fit(params):
                fit=gaussian(x,params)
                return (fit-y)      

            bounds=([-100, -100, -100], [100, 100, 100])
            fit=least_squares(gaussian_fit,
                              [1.0, self.dphie, 1.0], 
                              bounds=bounds)
            
            mu1=fit['x'][1]
            
            y_fit=gaussian(x, fit['x'])
            
            return y_fit, mu1


        def plotAll(self):
            if not os.path.isdir(self.dirname):
                os.mkdir(self.dirname)
                
            output_root=self.dirname+'/all_dphis'
            self.plotStuff(self.dphi_flat, output_root)
            
        def plotSeam(self):
            data=self.dphi_flat[self.num_pfs-1::self.num_pfs].copy()
            output_root=self.dirname+'/seam_dphis'
            self.plotStuff(data, output_root)
            
        def plotNonseam(self):
            data=self.dphi_flat.copy()
            data=np.delete(data,np.s_[self.num_pfs-1::self.num_pfs])
            output_root=self.dirname+'/nonseam_dphis'
            self.plotStuff(data, output_root)
            
        def plotIndiv(self):
            dirname='indiv_mt_plots'
            if not os.path.isdir(dirname):
                os.mkdir(dirname)
            
            start=5
            line=np.ones(self.num_pfs)*start
            spacing=start*2
            y=range(self.num_pfs)
            dphie=360.0/self.num_pfs
            
            ylim=[y[0], y[-1]]
            
            
                
            for key in self.dphixs.keys():
                plt.figure(figsize=(8,5))
                plt.ylim(ylim)
                name=key[0].split('/')[-1].split('.')[0]+'_'+key[-1]
                filename=dirname+'/'+name+'.jpg'
                
                data=self.dphixs[key]-dphie+start
                
                xlim=[data[0].min()-1, (data[-1]+1+(len(data)-1)*spacing).max()]
                plt.xlim(xlim)
                
                plt.xlabel(r'$\Delta\Phi (\degree)$')
                plt.ylabel('Protofilament #')
                
                for i in range(len(data)):
                        x_data=data[i,:]+spacing*i
                        x_ref=line+spacing*i
                        plt.plot(x_data, y, '-k')
                        plt.plot(x_ref, y, '--b')
                        
                
                plt.savefig(filename)
                plt.close()
            