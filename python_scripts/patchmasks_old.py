from basics import Basics
import os
import cupy as np
from util3d import get_center_of_pdb
from helix.basics import helical_wedge_mask
#from old_basics import helical_wedge_mask
from EMImage import EMImage
from util3d import spherical_cosmask
#from rotate import rotshift3D_spline
from cupyx.scipy.ndimage.interpolation import affine_transform

class RemoveWedges(Basics):
    def __init__(self,info_file='info.txt',shift=False):
        self.readInfo(info_file)
        self.subunitsToUse(shift)

    def maskVolume(self, volume, soft_mt_mask, ref_pdb=None, num_subunits=3, \
                   width=3):
        ###Use the mt mask from relion refinement
        vol=EMImage(volume)
        mtmask=EMImage(soft_mt_mask)
        vol.mult(mtmask)
        del mtmask
        vol_dim=np.asarray(vol.data.shape)
        self.createWedge(vol_dim=vol_dim,num_subunits=1, width=1, ref_pdb=ref_pdb,pf_num=0)
        self.soft_wedge.write_mrc('pf0_wedge.mrc')
        pf_num=1
        self.createWedge(vol_dim=vol_dim,num_subunits=num_subunits, width=width, \
                         ref_pdb=ref_pdb,pf_num=pf_num)
        z_shifts=(self.subunits_to_use[:,2]-self.subunits_to_use[pf_num,2])/self.pixel_size
        rot=self.twist_per_subunit*np.linspace(0,self.num_pfs-1,self.num_pfs)- \
            self.twist_per_subunit*pf_num
        sizex=self.soft_wedge.data.shape[0]
        padded=self.soft_wedge.clip(int(np.ceil(sizex*1.25)))
        for i in range(self.num_pfs):
            ##Data structure is set up weird. need to transpose data in order to 
            ##get desired phi rotation
            if i==0:
                mask=(self.rotshift3D_spline(padded.data.transpose(), \
                                             np.array([rot[i],0,0]))).transpose()
                #print type(z_shifts),z_shifts, type(z_shifts[i]),z_shifts[i]
                mask=np.roll(mask,int(round(z_shifts[i])),axis=0)
                temp=EMImage(mask)
                del mask
                temp=temp.clip(sizex)
            elif i==1:
                temp=self.soft_wedge
            elif i==self.num_pfs-1:
                pf_num=self.num_pfs-2
                self.soft_wedge=EMImage('pf%g/pf%g_wedge.mrc'%(i-1,i-1))
                #self.createWedge(vol_dim=vol_dim,num_subunits=num_subunits, width=width, \
                #                 ref_pdb=ref_pdb, pf_num=i-1)
                z_shifts=(self.subunits_to_use[:,2]-self.subunits_to_use[pf_num,2])/self.pixel_size
                rot=self.twist_per_subunit*np.linspace(0,self.num_pfs-1,self.num_pfs)- \
                    self.twist_per_subunit*pf_num
                padded=self.soft_wedge.clip(int(np.ceil(sizex*1.25)))
                mask=(self.rotshift3D_spline(padded.data.transpose(), \
                                             np.array([rot[i],0,0]))).transpose()
                #print type(z_shifts),z_shifts, type(z_shifts[i]),z_shifts[i]
                mask=np.roll(mask,int(round(z_shifts[i])),axis=0)
                temp=EMImage(mask)
                del mask
                temp=temp.clip(sizex)
            else:
                self.createWedge(vol_dim=vol_dim,num_subunits=num_subunits, width=width, \
                                 ref_pdb=ref_pdb, pf_num=i)
                temp=self.soft_wedge
            if not os.path.isdir('pf%g'%i):
                os.mkdir('pf%g'%i)
            temp.write_mrc('pf%g/pf%g_wedge.mrc'%(i,i))
            (temp.add(-1)).mult(-1)
            vol.mult(temp)
            vol.write_mrc('pf%g/pf_masked.mrc'%(i))
            del vol, temp
            vol=EMImage(volume)
            
    def pfModel(self, volume, ref_pdb, soft_mt_mask):
        vol=EMImage(volume)
        mtmask=EMImage(soft_mt_mask)
        vol.mult(mtmask)
        del mtmask
        subunits_to_use=self.subunitsToUse()
        vol_dim=np.asarray(vol.data.shape)
        num_repeats=int(1+np.ceil(vol_dim[0]*self.pixel_size/self.dimer_repeat_dist))
        self.createWedge(vol_dim,num_repeats,1,ref_pdb=ref_pdb,pf_num=0)
        vol.mult(self.soft_wedge)
        vol.write_mrc('protofilament.mrc')

    def patchPfModel(self,volume,ref_pdb,soft_mt_mask):
        vol=EMImage(volume)
        mtmask=EMImage(soft_mt_mask)
        vol.mult(mtmask)
        del mtmask
        subunits_to_use=self.subunitsToUse()
        vol_dim=np.asarray(vol.data.shape)
        num_subunits=int(1+np.ceil(vol_dim[0]*self.pixel_size/self.dimer_repeat_dist))
        com=get_center_of_pdb(ref_pdb)
        wedge=np.zeros((vol_dim[0],vol_dim[1],vol_dim[2]))
        subunit_num=self.subunits_to_use[0,0]*self.num_pfs
        num=[]
        for i in range(num_subunits):
            num.append(((-1)**i)*round(float(i)/2)*self.num_pfs+subunit_num)
            
        num.sort()
        print num
        for i in range(3):
            num.pop(num_subunits/2-1)
            
        print num
        for i in num:
            print i
            wedge[:]=wedge+helical_wedge_mask(vol_dim,self.pixel_size,self.rise_per_subunit,\
                                                 self.twist_per_subunit,num_pfs=self.num_pfs,\
                                                 num_starts=self.num_starts,ref_com=com,subunit_num=i)
            
        mask=wedge
        edge_resolution=20
        edge_width = self.pixel_size * np.ceil(edge_resolution/(2*self.pixel_size))
        cosmask_filter = np.fft.fftshift(spherical_cosmask(vol_dim, 0, edge_width / self.pixel_size))
        cosmask_filter = np.fft.fftn(cosmask_filter) / np.sum(cosmask_filter)
        soft_m = np.real(np.fft.ifftn( cosmask_filter * np.fft.fftn(mask)))
        
        soft_m[soft_m<0]=0
        ##Need to transpose to be on the same axes as before
        soft_wedge=EMImage(soft_m.transpose())
        vol.mult(soft_wedge)
        vol.write_mrc('patch_proto.mrc')

        
    def rotshift3D_spline(self, v, eulers, shifts=np.array([0,0,0]), mode='wrap'):
    # With a nod to:
    #  http://stackoverflow.com/questions/20161175/how-can-i-use-scipy-ndimage-interpolation-affine-transform-to-rotate-an-image-ab
        print shifts
        print eulers
        rot_origin = 0.5*np.array(v.shape)
        rot_rad = -np.deg2rad(eulers)
        phi = np.array([[np.cos(rot_rad[0]), np.sin(rot_rad[0]), 0],
                               [-np.sin(rot_rad[0]),np.cos(rot_rad[0]), 0],
                               [0                , 0              , 1]])
        theta=np.array([[np.cos(rot_rad[1]), 0,-np.sin(rot_rad[1])],[0,1,0],\
                            [np.sin(rot_rad[1]), 0, np.cos(rot_rad[1])]])
        psi=np.array([[np.cos(rot_rad[2]), np.sin(rot_rad[2]), 0],\
                          [-np.sin(rot_rad[2]),np.cos(rot_rad[2]), 0],\
                          [0,0,1]])
        rot_matrix=np.dot(np.dot(phi,theta),psi)
        offset = -(rot_origin-rot_origin.dot(rot_matrix)).dot(np.linalg.inv(rot_matrix))
        offset = offset - shifts
        transformed_v = affine_transform(v,rot_matrix,offset=offset,mode=mode)
    #    print '1done'
    #    transformed_v = scipy.ndimage.interpolation.shift(transformed_v,shifts,mode=mode)
       
    
        return transformed_v

    def createWedge(self, vol_dim, num_subunits, width, ref_pdb, pf_num=6, edge_resolution=20):
        ###Create a wedge mask centered at pf6 to be applied to reference volume
        com=get_center_of_pdb(ref_pdb)
        wedge=np.zeros((vol_dim[0],vol_dim[1],vol_dim[2]))
        subunit_num=self.subunits_to_use[pf_num,0]*self.num_pfs+pf_num
        for i in range(num_subunits):
            num=((-1)**i)*round(float(i)/2)*self.num_pfs+subunit_num
            print num
            wedge[:]=wedge+helical_wedge_mask(vol_dim,self.pixel_size,self.rise_per_subunit,\
                                                 self.twist_per_subunit,num_pfs=self.num_pfs,\
                                                 num_starts=self.num_starts,ref_com=com,subunit_num=num)
            for j in range(1,width):
                num2=num+((-1)**j)*round(j/2.0)
                wedge[:]=wedge+helical_wedge_mask(vol_dim,self.pixel_size,self.rise_per_subunit,\
                                                     self.twist_per_subunit,num_pfs=self.num_pfs,\
                                                     num_starts=self.num_starts,ref_com=com,subunit_num=num2)
                print num2
        mask=wedge
        edge_width = self.pixel_size * np.ceil(edge_resolution/(2*self.pixel_size))
        cosmask_filter = np.fft.fftshift(spherical_cosmask(vol_dim, 0, edge_width / self.pixel_size))
        cosmask_filter = np.fft.fftn(cosmask_filter) / np.sum(cosmask_filter)
        soft_m = np.real(np.fft.ifftn( cosmask_filter * np.fft.fftn(mask)))
        
        soft_m[soft_m<0]=0
        ##Need to transpose to be on the same axes as before
        self.soft_wedge=EMImage(soft_m.transpose())

        
    def maskSubunit(volume,subunit_pdb,num_subunits=1,num_pfs=13,num_starts=1.5, rise_per_subunit=9.4,\
                    width=1,voxel_size=1.247,edge_resolution=20,subunit_num=0,\
                    helical_twist=None,save=False,pf_num=None,new_name=None):


        mask=wedge
        edge_width = voxel_size * np.ceil(edge_resolution/(2*voxel_size))
        cosmask_filter = np.fft.fftshift(cf.spherical_cosmask(vol_dim, 0, edge_width / voxel_size))
        cosmask_filter = np.fft.fftn(cosmask_filter) / np.sum(cosmask_filter)
        soft_m = np.real(np.fft.ifftn( cosmask_filter * np.fft.fftn(mask)))
        
        inv=1-soft_m
        masked=inv*vol
        subunit=soft_m
        subunit[subunit<0]=0
        if new_name is None:
            new_name=volume
            
        if save:
            second_name='wedge_'+new_name
            new_name='masked_'+new_name
        #    third_name='removed_subunit.mrc'
            gc.saveImage(masked,new_name)
            gc.saveImage(subunit,second_name)
        #    gc.saveImage(inv_masked,third_name)
        else:
            return soft_m, masked
        

                
        