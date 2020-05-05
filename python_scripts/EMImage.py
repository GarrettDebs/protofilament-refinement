''' EMImage class.
    class of image represent an electron microscopy image, 1d, 2d or 3d

'''
import sys
import os
import logging
import warnings
import mrcfile
import numpy as np
import scipy.ndimage
import scipy.ndimage.fourier as fourier
import scipy.misc

# image type constant
FORMAT_MRC = 1
FORMAT_SPIDER = 2
FORMAT_SPIDER_SINGLE = 3
FORMAT_MRCS = 4
FORMAT_JPEG = 5

logger = logging.getLogger(__name__)
# logger.setLevel(logging.DEBUG)
# logger.addHandler(logging.StreamHandler(sys.stdout))

class _EMImage(object):
    '''The EM Image class.
        If don't want to modify original instance, use .copy().
        The data is represented with numpy array, and
        can be accessed by img.data. (numpy.ndarray)
        To get the dimensions, use img.get_xsize(), img.get_ysize(), img.get_zsize()
        x is column, y is row, z is slice.
        EMImage np array is default to np.float, which is 32 bit float number.
        This is made compatable to mrc file. If 64 bits double is needed,
        the save routine should specify the type
    '''
    def __init__(self, filenameornx, ny=None, nz=None, 
                 fmt=None, isFFT=False, mmap=False, 
                 mmap_file=None, mode='r', dtype=np.float32):
        '''
            self.data: The numpy.ndarray of the image. With numpy index convention.
            Order of index is slice, row, col.
            usage:
            img = EMImage("filename.mrc") # create a EMImage using filename
            # create a mmap image, which only read file when needed. 
            # Default mode is read only
            img=  EMImage("filename.mrc", mmap=True) 
            
            # create a new mmap image, with dimensions. 
            # This is used when the memory is not large enough for the volume
            img = EMImage(128, 128, 128, mmap=True, 
                        mmap_file="/tmp/a.mrc", mode='w+')
            
            img = EMImage(128, 128, 128) # create a empty 3D volume
            
            img = EMImage(128, 128) # create a empty 2D image
            
            # create a image with its fft data
            img = EMImage(fftdata, isFFT = True)
            
            # create a image with numpy array. borrowed reference of nparray

            img = EMImage(nparray)        
            
            Access numpy array data:
            fftdata = img.fftdata
            rdata = img.rdata
            data = img.data  # depending on the isFFT flag. 
            
            Set numpy data:
            img.fftdata = otherfftdata
            img.rdata = otherrdata
            img.data = otherdata # depending on the isFFT flag 
        
        Args:
            filenameornx: A `str` or number or ndarray.
                if is str, then read from file.
                if is int, then this is the x dimension, or number of columns
            ny: A `int`. The y dimension of image (row numbers)
            nz: A `int`. The z dimension. (slice numbers)
        '''
        self.headers = {}
        self._isfft = isFFT
        self._mmap = mmap # memmapped file, used for large file access
        self._mmap_file = mmap_file
        self._mmap_mode = mode
        self.dtype = dtype
        # The _voxel size is used for calculate the resolution.
        self._voxel_size = 1
        # data. 
        self._fftdata = None
        self._rdata = None
       
        if isinstance(filenameornx, str):
            self.filename = filenameornx
            self.read_image(self.filename, fmt=fmt)
        elif self._mmap and self._mmap_file is None:
            raise Exception("Memmap file must have filename assoicated. pass mmap_file as parameter")
        elif type(filenameornx) is int or type(filenameornx) is float:
            # this is initialization of empty image
            self._init_data(filenameornx, ny, nz)
        elif isinstance(filenameornx, list) or isinstance(filenameornx, tuple):
            # The dimention is given as list of tuple
            self._init_data(*filenameornx)
        elif isinstance(filenameornx, np.ndarray): 
            # memmap array is subclass of ndarray
            if isFFT:
                self.fftdata = filenameornx
            else:
                self.rdata = filenameornx
        else:
            raise Exception("Invalid constructor")

    def _init_data(self, nx, ny=None, nz=None):
        ''' Init data with zeros. All arguments are truncated to int silently.

        Args:
            nx: A `int`. The x dimension of the image
            ny: A `int`. The y dimension of the image. default is 1
            nz: A `int`. The z dimension of the image. Default is 1
        '''
        nx = int(nx)
        if nx <= 0:
            raise ValueError("Wrong value for x dimension")
        dims = [nx]
        if ny is not None:
            ny = int(ny)
            if ny <= 0:
                raise ValueError("Wrong value for y dimension")
            dims = [ny] + dims
            if nz is not None:
                nz = int(nz)
                if nz <= 0:
                    raise ValueError("Wrong value for z dimension")
                dims = [nz] + dims
        if not self._mmap:
            self.data = np.zeros(dims, self.dtype)
        else:
            mrc = mrcfile.mmap(self._mmap_file, mode=self._mmap_mode)
            mrc._open_memmap(np.float32, tuple(dims))
            mrc.update_header_from_data()
            self.data = mrc.data

    @property
    def voxelSize(self):
        return self._voxel_size

    @voxelSize.setter
    def voxelSize(self, value):
        ''' Set voxel Size
        '''
        if value <= 0.:
            return
        self._voxel_size = value

    @property
    def ndim(self):
        ''' Return the number of dimensions
        @return: The number of dimensions
        '''
        if self._rdata is not None: return self._rdata.ndim
        return self._fftdata.ndim

    @property
    def center(self):
        ''' Get center of image. size / 2
        '''
        return [x / 2 for x in self.get_size()]

    def read_image(self, filename, fmt=None):
        ''' Read an image from file. Support mrc format through mrcfile
        '''
        if fmt is None:
            fmt = self.guess_file_format(filename)
        if fmt == FORMAT_MRC or fmt.lower() == "mrc":
            self.read_mrc(filename)
        return self

    @staticmethod
    def guess_file_format(filename):
        ''' Guess file type from file name
        '''
        if filename.lower().endswith(".mrc") or filename.lower().endswith(".mrcs"):
            return FORMAT_MRC
        elif filename.lower().endswith(".jpeg") or filename.lower().endswith(".jpg"):
            return FORMAT_JPEG
        return None

    def read_mrc(self, filename):
        ''' Read mrc file through
        '''
        # permissive is to go throught invalid mrc file like generated by
        # motioncor2.... MotionCor2 does not write "MAP" into header...
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            mrc = mrcfile.mmap(filename, self._mmap_mode, permissive = True)
            if self._mmap:
#               logger.debug("Using mmaped file: %s" % filename)
                self._mmap_file = filename
                self.data = np.memmap(self._mmap_file,
                        dtype = mrc.data.dtype,
                        mode = self._mmap_mode,
                        offset = mrc.data.offset,
                        shape = mrc.data.shape
                        )
            else:
                self.data = mrc.data.copy()
        return self

    def get_xsize(self):
        ''' get x size, number of columns
        '''
        if self._rdata is not None:
            return self._rdata.shape[-1]
        return self._fftdata.shape[-1]

    def get_ysize(self):
        ''' get y size. number of rows
        '''
        if self.ndim < 2: return 1
        if self._rdata is not None: return self._rdata.shape[-2]
        return self._fftdata.shape[-2]

    def get_zsize(self):
        ''' get z size. number of slices
        '''
        if self.ndim < 3: return 1
        if self._rdata is not None: return self._rdata.shape[-3]
        return self._fftdata.shape[-3]

    def get_size(self):
        '''
            return the size as tuple.
            (nx, ny, nz) if 3 dim, (nx, ny) if 2 dim, (nx,) if one dim
        '''
        if self._fftdata is not None:
            return self._fftdata.shape[::-1]
        else:
            return self._rdata.shape[::-1]

    def write_image(self, filename, fmt=None):
        '''
            Write image to disk.
        '''
        if fmt is None:
            fmt = self.guess_file_format(filename)
        if (fmt == FORMAT_MRC) or (isinstance(fmt, str) and fmt.lower() == "mrc"):
            self.write_mrc(filename)
        elif ((fmt == FORMAT_JPEG)
              or (isinstance(fmt, str) and (fmt.lower() == "jpg"
                                            or fmt.lower() == "jpeg"))):
            self.write_jpeg(filename)
        else:
            raise Exception("EMImage does not support format: %s" % fmt)
        return self
    def write_jpeg(self, filename):
        ''' Write as jpeg format to file name
        '''
        if self.ndim > 2 and (self.ndim == 3 and self.get_zsize() != 1):
            raise Exception("EMImage does not support 3D for jpeg")
        if self.ndim == 3:
            data = self.data.reshape([self.get_ysize(), self.get_xsize()])
        else:
            data = self.data
        scipy.misc.toimage(data).save(filename)  # @UndefinedVariable
        return self

    def write_mrc(self, filename):
        ''' Write as mrc format to file name
        '''
        #header = self.get_mrc_header()
        # mmaped file already written to file. no need to write again.
        if self._mmap: 
            return
        mrc = mrcfile.new(filename, overwrite=True)
        if self.data.dtype != np.float32:
            array_data = self.data.astype(np.float32, order='C')
        else:
            array_data = self.data.copy(order='C')
        if array_data.ndim == 1:
            array_data = array_data.reshape([array_data.shape[0], 1])
        mrc.set_data(array_data)
        mrc.voxel_size = self.voxelSize
        mrc.header.nxstart, mrc.header.nystart, mrc.header.nzstart = (-self.xsize / 2, -self.ysize/2, -self.zsize/2)
        mrc.close()

    def append_to_mrc(self, filename):
        ''' Append data to mrc file
        '''
        # if file does not exists yet, write new mrc file
        if not os.path.exists(filename): 
            self.write_mrc(filename)
            return
        # use memmapped mrcfile to append to mrc
        mrc = mrcfile.mmap(filename, mode='r+')
        shape = mrc.data.shape
        shape = list(shape)
        data = self.data
        
        # Get shape informations
        # promote dimension if 1d append to 1d or 2d append to 2d
        if ((self.ndim == 1 and len(shape) == 1)
                or (self.ndim == 2 and len(shape) == 2)):
            shape = [1] + shape
            data = data.reshape(shape)
            shape[0] += 1
        # if 3d append to 2d or 2d append to 1d:
        elif ((self.ndim == 2 and len(shape)==1)
                or (self.ndim == 3 and len(shape) == 2)):
            shape = [self.data.shape[0] + 1] + shape
        # if 1d append to 2d or 2d append to 3d:
        elif ((self.ndim == 2 and len(shape) == 3) 
                or (self.ndim == 1 and len(shape) == 2)):
            shape[0] += 1
            data = data.reshape([1] + list(self.data.shape))
        elif ((self.ndim == 3) and len(shape) == 3):
            shape[0] = shape[0] + self.data.shape[0]
        else:
            raise Exception("Does not support %dD append to %dD" % (self.ndim, len(shape)))
        # append data
        mrc._open_memmap(dtype = mrc.data.dtype, shape = tuple(shape))
        mrc.data[shape[0] - data.shape[0]:, :, :] = data
        mrc.update_header_from_data()
        return self


    # addition
    def add(self, other):
        ''' Add other image / constant to this image. Inplace
        '''
        if isinstance(other, EMImage):
            assert self.same_shape(other)
            other = other.data

        np.add(self.data, other, self.data, casting='safe')

        return self

    def __add__(self, other):
        tmp = self.copy()
        tmp.add(other)
        return tmp
    def __radd__(self, other):
        return self.__add__(other)
    def __iadd__(self, other):
        self.add(other)
        return self

    # subtraction
    def sub(self, other):
        ''' Subtract other image / constant from this image. Inplace
        '''
        if isinstance(other, EMImage):
            assert self.same_shape(other)
            other = other.data
        self.data -= other
        return self

    def __sub__(self, other):
        return self.copy().sub(other)

    def __rsub__(self, other):
        tmp = self.copy()
        tmp.mult(-1)
        tmp.add(other)
        return tmp

    def __isub__(self, other):
        self.sub(other)
        return self
    # multiplication
    def mult(self, other):
        ''' Multiply other image / constant to this image. Inplace
        '''
        if isinstance(other, EMImage):
            assert self.same_shape(other)
            other = other.data
        self.data *= other

    def __imul__(self, other):
        self.mult(other)
        return self
    def __mul__(self, other):
        tmp = self.copy()
        tmp.mult(other)
        return tmp
    def __rmul__(self, other):
        return self.__mul__(other)

    # division
    def div(self, other):
        if isinstance(other, EMImage):
            assert self.same_shape(other)
            other = other.data
        self.data /= other
    def __div__(self, other):
        tmp = self.copy()
        tmp.div(other)
        return tmp
    def __idiv__(self, other):
        self.div(other)
        return self

    @property
    def dim(self):
        ''' return dimensions in (x,y,z,...)
        '''
        return self.get_size()

    # data processing
    def insert_image(self, other, position):
        '''
            insert other image to this image, put the center other image in position of this image
            other image should have same number dimensions of this image
        '''
        self.paste(other, center1=position)

    def insert_sum(self, other, position, scale = 1.):
        ''' sum other at position,
        '''
        tmp = other.copy()
        tmp.paste(self, center2=position)
        tmp += other * scale
        self.paste(tmp, center1=position)

    def get_data(self, x_range=None, y_range=None, z_range=None):
        ''' Get image data as numpy array.
        @param x_range: The x range, tuple [start, end)
        @param y_range: The y range, tuple [start, end)
        @param z_range: The z range, tuple [start, end)
        '''
        x_range = x_range or [0, self.get_xsize()]
        y_range = y_range or [0, self.get_ysize()]
        z_range = z_range or [0, self.get_zsize()]
        if self.dim == 3:
            return self.data[z_range[0]:z_range[1], y_range[0]:y_range[1], x_range[0]:x_range[1]]
        elif self.dim == 2:
            return self.data[y_range[0]:y_range[1], x_range[0]:x_range[1]]
        elif self.dim == 1:
            return self.data[x_range[0]:x_range[1]]
        else:
            return None
    def shift(self, shx=0, shy=0, shz=0):
        ''' Shift the image
        '''
        return self.fourier_shift(shx, shy, shz)

    def shift_i(self, ix=0, iy = 0, iz = 0):
        ''' Shift image with integer pixels
        '''
        shft = []
        if self.ndim == 1:
            shft, axis = [ix], [0]
        elif self.ndim == 2:
            shft, axis = [iy, ix], [0, 1]
        elif self.ndim == 3:
            shft, axis = [iz, iy, ix], [0, 1, 2]
        self.rdata = np.roll(self.rdata, shft, axis)
        return self

    def fourier_shift(self, shx=0, shy=0, shz=0):
        ''' shift image in fourier space. Inplace
        '''
        shft = []
        if self.ndim == 1:
            shft = [shx]
        elif self.ndim == 2:
            shft = [shy, shx]
        elif self.ndim == 3:
            shft = [shz, shy, shx]
        else:
            raise Exception("Dimention %d not implemented for fourier_shift" 
                            % self.ndim)
        if self._isfft:
            self.data = fourier.fourier_shift(self.data, shft)
        else:
            self.data = EMImage(fourier.fourier_shift(self.fftdata, shft), 
                                isFFT=True).rdata
        return self

    def get_index_center(self):
        ''' get centered indexes
        '''
        nx = self.get_xsize()
        ny = self.get_ysize()
        nz = self.get_zsize()
        idxx = np.arange(-(nx / 2), (nx+1)/2)
        idxy = np.arange(-(ny / 2), (ny+1)/2)
        idxz = np.arange(-(nz / 2), (nz+1)/2)
        if self.ndim == 1:
            return idxx
        if self.ndim == 2:
            return np.meshgrid(idxx, idxy)
        if self.ndim == 3:
            return np.meshgrid(idxx, idxy, idxz)

    def resample(self, new_size=None, scale=None):
        '''
            resample image to new size
            if new_size is given, scale is ignored.
            if new_size is not given, scale is used.
            if new_size or scale is scalar, aspect_ratio is kept
            otherwise len(new_size) or len(scale) must be same as self.ndim
            @param new_size: The new size of the image, int
            @param scale: The scale factor, result is ceil(ori_size * scale)
        '''
        if new_size is None:
            if scale is None:
                logger.warn("Image is not resampled due to both par is None")
                return
            else:
                try:
                    new_size = list(self.size)
                    new_size = [new_size[i]*scale[i] for i in range(self.ndim)]
                except TypeError:
                    new_size = [new_size[i]*scale for i in range(self.ndim)]
        try:
            _ = new_size[self.ndim-1]
        except TypeError:
            scale = float(new_size) / self.xsize
            new_size = [self.size[i] * scale for i in range(self.ndim)]

        except IndexError:
            raise Exception("The scale is provided for each dimension, but not enough")
        new_size = map(int, new_size)
        scale = float(new_size[0]) / self.xsize
        if abs(scale - 1) < 0.001:
            return
        # pad/clip fft
#        d = self.data
#        sv = np.sum(d[0,:]) + np.sum(d[-1,:]) + np.sum(d[1:-1,0]) + np.sum(d[1:-1,-1])
#        if self.ndim == 2:
#            sv /= sum(self.size) * 2 - 4
#        elif self.ndim == 3:
#            sv /= sum([_*_ for _ in self.size]) * 4 - 16
        from cryofilia.filter import FourierFilter as ff
        fr = 2./scale 
        edge = fr / 10
        rl, rh = fr + edge, fr - edge
        bl, bh = 1. / rl, 1. / rh
        tmp = self.copy()
        ff.butterworth(bl, bh, image=tmp)
        if self.xsize > new_size[0]:
            self.data = tmp.decimate(self.xsize / new_size[0], average=True).data
            self.voxelSize /= scale

    def decimate(self, factor, average = False):
        '''
            Downsample the image with factor, factor must be integer
        '''
        if not isinstance(factor, int):
            logger.warning("factor must be integer, try to convert to integer")
            factor = int(factor)
        if factor <= 1:
            return
        if self.ndim == 1:
            arr = self.data[::factor].copy()
            if average:
                for i in range(1, factor):
                    tmp = self.data[i::factor]
                    arr[:tmp.shape[0]] += tmp
                arr /= factor
            return EMImage(arr)
        elif self.ndim == 2:
            arr = self.data[::factor,::factor].copy()
            if average:
                for i in range(0, factor):
                    for j in range(0, factor):
                        if i == 0 and j == 0: continue
                        tmp = self.data[i::factor, j::factor]
                        arr[:tmp.shape[0], :tmp.shape[1]] += self.data[i::factor,j::factor]
                arr /= factor * factor
            return EMImage(arr)
        elif self.ndim == 3:
            arr = self.data[::factor,::factor,::factor].copy()
            if average:
                for i in range(0, factor):
                    for j in range(0, factor):
                        for k in range(0, factor):
                            if i == 0 and j == 0 and k == 0: continue
                            tmp = self.data[i::factor, j::factor]
                            arr[:tmp.shape[0], :tmp.shape[1], :tmp.shape[2]] += self.data[i::factor,j::factor, k::factor]
                arr /= factor * factor * factor
            return EMImage(arr)
        
    def inverse(self):
        '''  inverve the density
        '''
        self.data *= -1

    def peakfind2d(self, mask=None):
        ''' Find peak from the image
        '''
        from cryofilia import util2d
        offset = np.min(self.data[:])
        if mask is not None:
            data = self.data * mask - offset
        else:
            data = self.data - offset
        x_fit, y_fit, ccc = util2d.peakfind2d(data)
        return x_fit, y_fit, ccc + offset

    def copy(self):
        '''
        make a copy of this image
        '''
        from copy import copy as _copy
        data = self.data.copy()
        copied = EMImage(data, mmap=self._mmap, mode=self._mmap_mode)
        copied.headers = _copy(self.headers)
        return copied

    def paste(self, other, center1=None, center2=None):
        ''' paste other image to this image
        '''
        if self.ndim != other.ndim:
            raise Exception("Dimension not match: %s vs %s" % (self.ndim, other.ndim))
        mydim = self.dim
        otherdim = other.dim

        if center1 is None:
            center1 = [_x / 2 for _x in mydim]
        if center2 is None:
            center2 = [_x / 2 for _x in otherdim]
        if isinstance(center1, float):
            center1 = int(center1)
        if isinstance(center2, float):
            center2 = int(center2)
        if isinstance(center1, float):
            center1 = [center1] * self.ndim
        if isinstance(center2, float):
            center2 = [center2] * self.ndim
        ovps = [self.get_overlap(mydim[i], center1[i], otherdim[i], center2[i])
                for i in range(self.ndim)]
        rgs1, rgs2 = zip(*ovps)
#        if self.data.dtype.type is np.complex128:
#            asd
        if self.ndim == 1:
            self.data[rgs1[0][0]:rgs1[0][1]] = other.data[rgs2[0][0]:rgs2[0][1]]
        elif self.ndim == 2:
            self.data[rgs1[1][0]:rgs1[1][1],
                      rgs1[0][0]:rgs1[0][1]] = other.data[rgs2[1][0]:rgs2[1][1],
                                                          rgs2[0][0]:rgs2[0][1]]
        elif self.ndim == 3:
            self.data[rgs1[2][0]:rgs1[2][1],
                      rgs1[1][0]:rgs1[1][1],
                      rgs1[0][0]:rgs1[0][1]] = other.data[rgs2[2][0]:rgs2[2][1],
                                                          rgs2[1][0]:rgs2[1][1],
                                                          rgs2[0][0]:rgs2[0][1]]

    @property
    def rfftdata(self):
        if self.ndim == 1:
            return self.fftdata[:self.xsize / 2 + 1]
        if self.ndim == 2:
            return self.fftdata[:, :self.xsize / 2 + 1]
        if self.ndim == 3:
            return self.fftdata[:, :, :self.xsize / 2 + 1]
        raise RuntimeError("dimension: %d not supported" % self.ndim)

    @property
    def fftdata(self):
        ''' FFT data getter
        '''
        if self._fftdata is not None:
            return self._fftdata
        # calculate fft
        self._fftdata = self._fft(self._rdata)
        return self._fftdata

    @fftdata.setter
    def fftdata(self, value):
        """ set fft data will invalid the real data, if not the same
        """
        if self._fftdata is not value:
            self._fftdata = value
            self._rdata = None

    @property
    def rdata(self):
        ''' Property of real data
        '''
        if self._rdata is not None:
            return self._rdata
        # calculate the rdata
        self._rdata = self._ifft(self._fftdata).real
        return self._rdata

    @rdata.setter
    def rdata(self, value):
        ''' set real data. if real data changed, reset fft data
        @param value: set real data value
        '''
        if self._rdata is not value:
            self._rdata = value
            self._fftdata = None


    def rotate2d(self, angle, shifts=np.array([0, 0]), mode='wrap', order=3):
        '''
            Rotate 2D image.
            psi is the rotation angle, in degree, CCC is positive
        '''
        rot_origin = np.array([self.get_ysize() / 2, self.get_xsize() / 2])
        rot_rad = -angle * np.pi / 180.
        rot_matrix = np.array([[np.cos(rot_rad), np.sin(rot_rad)],
                               [-np.sin(rot_rad), np.cos(rot_rad)]])
        offset = -(rot_origin - rot_origin.dot(rot_matrix)).dot(np.linalg.inv(rot_matrix))
        offset = offset - shifts
        if self.ndim == 2:
            transformed_img = scipy.ndimage.interpolation.affine_transform(
                self.data, rot_matrix, offset=offset, mode=mode, order=order)
            self.data = transformed_img

        elif self.ndim == 3:
            for islice in range(self.get_zsize()):
                transformed_img = scipy.ndimage.interpolation.affine_transform(
                    self.data[islice, :, :], rot_matrix, offset=offset, mode=mode)
                self.data[islice, :, :] = transformed_img[:, :]

    def ramp(self):
        ''' remove the ramp from image. Only work in real image
        '''
        if self.isFFT:
            raise Exception("Could not do ramp on fft image")
        from cryofilia.util2d import ramp as ra
        tmp = ra(self.data)
        self.data = tmp

    def normalize(self, stdv=1):
        ''' normalize the image
        '''
        ori_std = np.std(self.data)
        ori_mean = np.average(self.data)
        self.data -= ori_mean
        self.data /= ori_std / stdv

    def same_shape(self, other):
        ''' Check if the image has same shape of other image
        '''
        return self.size == other.size

    def clip(self, sizex, center=None, padValue=0):
        '''
            used to clip a volume to smaller size, or pad to bigger volume
        '''
        if isinstance(sizex, int):
            sizex = [sizex] * self.ndim
        if len(sizex) != self.ndim:
            raise Exception("new size must have same number of dimensions as ori: %s vs %s"
                            % (len(sizex), self.ndim))
        img = EMImage(sizex)
        if padValue != 0:
            img+=padValue
        img.paste(self, center2=center)
        return img

    def crop(self, sizex, center=None):
        ''' Crop the image with new size, same as clip
        '''
        return self.clip(sizex, center)

    def box_particles(self, coords, box_size=None, prerotate=None):
        '''
            box particles centerd at coords, and return the stack
            @param coords: The list of (x,y) coordinates, center of each ptcls
            @param box_size: The output box size.
            @return: The boxed ptcls
        '''
        nboxes = len(coords)
        ptcl_dim = len(coords[0])
        result = []
        for box_idx in range(nboxes):
            try:
                if ptcl_dim == 1:
                    ptcl = EMImage(box_size)
                elif ptcl_dim == 2:
                    ptcl = EMImage(box_size, box_size)
                elif ptcl_dim == 3:
                    ptcl = EMImage(box_size, box_size, box_size)
                coords[box_idx] = [_ - 1 for _ in coords[box_idx]]
                if prerotate and ptcl_dim > 1:
                    self.box_particle1(ptcl, coords[box_idx], prerotate=prerotate[box_idx])
                else:
                    self.box_particle1(ptcl, coords[box_idx])
                result.append(ptcl)
            except ValueError as ve:
                import traceback as tb
                tb.print_exc()
                logger.error("Could not box at coords: %s, with box_size %s, and prerotate: %s, from micrograph with size: %s" % (coords[box_idx], box_size, prerotate, self.size))
                raise ve
        return result

    @classmethod
    def get_overlap(cls, shape1, center1, shape2, center2, boxsize=-1):
        ''' Get overlap between shape1 and shape2
            center is defined as bs / 2
            @param shape1: the (start, end) for first range
            @param shape2: The (start, end) for second range
            @param center1: The center in first range
            @param center2: the center in second range
            @param bs: The expected size. -1 if to get all overlap
            @return: The range in shape1, the range in shape2.
                    [[left1, right1], [left2, right2]]
                    Left index is included, right index is excluded.
        '''
        if isinstance(shape1, int):
            shape1 = [0, shape1]
        if isinstance(shape2, int):
            shape2 = [0, shape2]
        new_shape1 = [_x - center1 for _x in shape1]
        new_shape2 = [_x - center2 for _x in shape2]
        left = max(new_shape1[0], new_shape2[0])
        if boxsize > 0:
            left = max(left, -(boxsize / 2))
        right = min(new_shape1[1], new_shape2[1])
        if boxsize > 0:
            right = min(right, (boxsize + 1) / 2)
        return ([_x + center1 for _x in [left, right]],
                [_x + center2 for _x in [left, right]])

    def box_particle1(self, ptcl, center, prerotate=None):
        ''' box particle to ptcl, centered at center.
            if center is fractional, the boxed ptcl would be fourier shifted.
            To avoid the shift, use whole number instead
        '''
        center_i = [int(c) for c in center]
        center_f = [_c - c for (_c, c) in zip(center, center_i)]
        ptcl.paste(self, center2=center_i)
        #if any([abs(_x) > 0.01 for _x in center_f]):
        #    ptcl.fourier_shift(*center_f)
        if prerotate is not None:
            ptcl.rotate2d(prerotate)

    @staticmethod
    def box_particle(ptcl, micrograph, centerx, centery, prerotate=0.):
        '''
            box particle from micrograph mg, with center cx, cy.
            to numpy 2d array ptcl.
            @param ptcl: The preallocated numpy array,
            @param micrograph: The numpy ndarray of the micrograph
            @param centerx: The center of ptcl in mg. Column
            @param centery: The center of ptcl in mg. Row
        '''
        ptcl_img = EMImage(ptcl)
        EMImage(micrograph).box_particle1(ptcl_img, [centerx, centery], prerotate=prerotate)

    def box(self, center, size):
        ''' Get boxed particle from image
        '''
        ptcl = EMImage(*size)
        self.box_particle1(ptcl, center)
        return ptcl

    def corrmap(self, other, fftshift=False):
        ''' Get cross correlation map with other image.
        Must be same dimension
        '''
        cmap = self._ifft(self.fftdata * other.fftdata.conj()).real
        if fftshift:
            cmap = np.fft.fftshift(cmap)
        return EMImage(cmap)

    def correlation(self, other, mask=None):
        ''' Get cross correlation coefficient
        '''
        if mask is None:
            return np.corrcoef(self.rdata.flatten(), other.rdata.flatten())[0, 1]
        else:
            return np.corrcoef(self.rdata[mask.data > 0.5].flatten(),
                               other.rdata[mask.data > 0.5].flatten())[0, 1]
    
    def gaussian_filter(self, resolution = None):
        import numpy as np
        N = self.xsize 
        r = self.voxel_size * N / resolution
        mask = TestImage.gaussian(self.size, r, 4)
        self.fftdata = self.fftdata * np.fft.ifftshift(mask.data)
        return self
    @property
    def xsize(self):
        ''' property xsize
        '''
        return self.get_xsize()

    @property
    def ysize(self):
        ''' property ysize
        '''
        return self.get_ysize()

    @property
    def zsize(self):
        ''' property zsize
        '''
        return self.get_zsize()

    @property
    def size(self):
        ''' property size
        '''
        return self.get_size()

    @property
    def data(self):
        ''' property data
        '''
        if self._isfft:
            return self.fftdata
        else:
            return self.rdata

    @data.setter
    def data(self, value):
        if self._isfft:
            self.fftdata = value
        else:
            self.rdata = value
    @property
    def _ifft(self):
        ''' ifft function for image
        '''
        if self.ndim == 1:
            return np.fft.ifft
        elif self.ndim == 2:
            return np.fft.ifft2
        elif self.ndim == 3:
            return np.fft.ifftn

    @property
    def _fft(self):
        ''' fft function for image fourier transform.
        Currently use numpy fft function
        '''
        if self.ndim == 1:
            return np.fft.fft
        elif self.ndim == 2:
            return np.fft.fft2
        elif self.ndim == 3:
            return np.fft.fftn

    @property
    def isFFT(self):
        ''' property isFFT
        '''
        return self._isfft
    @isFFT.setter
    def isFFT(self, value):
        ''' FFT flag setter. Convenient way to convert fft image to real image
        '''
        self._isfft = True if value else False
#     def insert_scaled_sum(self, *args, **kwargs):
#         from EMAN2 import EMNumPy
#         EMNumPy.numpy2em(self.data).insert_scaled_sum(*args, **kwargs)
    @staticmethod
    def get_image_size(image_file):
        ''' Get image size.
        @param image_file: Input image file
        '''
        return EMImage(image_file).get_size()
    def to_stack(self):
        """ Convert a 3D image to a list of 2D image
        """
        return [EMImage(self.data[i,:]) for i in range(self.zsize)]

    @staticmethod
    def stack(ptcls):
        ''' Stack the particles into volume
        Suppose all image in ptcls has same dimension. 
        Only works on real data domain, not the fft domain
        :return: A EMImage instance, has one dimension higher than ptcls
        '''
        if len(ptcls) == 0: return None
        nptcls = len(ptcls)
        stk_shape = list(ptcls[0].size) + [nptcls]
        result = EMImage(stk_shape)
        for ptcl_idx in range(nptcls):
            result.data[ptcl_idx, :] = ptcls[ptcl_idx].data
        return result

class TestImage(object):
    '''
        generate synthetic images using pattern
    '''
    @staticmethod
    def cylinder(box_size, radius=None, radius_2=None, axis='z'):
        ''' Generate a cylinder volume
        '''
        img = EMImage(box_size, box_size, box_size)
        if radius is None:
            radius = box_size / 2
        if radius_2 is None:
            radius_2 = 0
        bs2 = box_size / 2
        idx_x, idx_y = np.mgrid[-bs2:bs2, -bs2:bs2]
        idx_r = idx_x*idx_x + idx_y*idx_y
        radius2 = radius * radius
        r22 = radius_2 * radius_2
        idx_c = 1
        if axis == 'z':
            img.data[:, np.logical_and(idx_r <= radius2, idx_r >= r22)] = idx_c
        elif axis == 'x' or axis == 'y':
            img.data[np.logical_and(idx_r <= radius2, idx_r >= r22), :] = idx_c
            if axis == 'y':
                img.data = img.data.transpose((0, 2, 1))
        return img

    @staticmethod
    def sphere(box_size, radius=None):
        ''' Generate a sphere volume
        '''
        if radius is None:
            radius = box_size / 2
        bs2 = box_size / 2
        idx_z, idx_y, idx_x = np.mgrid[-bs2:bs2, -bs2:bs2, -bs2:bs2]
        idx_r = idx_x*idx_x + idx_y*idx_y + idx_z*idx_z
        radius2 = radius * radius
        idx_r[idx_r < radius2] = 1
        idx_r[idx_r >= radius2] = 0
        return EMImage(idx_r)
    
    @staticmethod
    def circle(box_size, radius=None, inner_radius = None):
        ''' Generate a circle in 2D
        '''
        if radius is None:
            radius = box_size / 2
        if inner_radius is None:
            inner_radius = 0
        bs = box_size
        idx_y, idx_x = np.mgrid[-(bs/2):(bs + 1)/2, 
                -(bs/2):(bs+1)/2]
        idx_r = idx_x*idx_x + idx_y*idx_y
        radius2 = radius * radius
        iradius2 = inner_radius * inner_radius
        r = np.zeros([box_size, box_size])
        r[idx_r < radius2] = 1 
        r[idx_r < iradius2] = 0
        return EMImage(r)
    
    @staticmethod
    def gaussian(box_size, radius = None, width = 0):
        try:
            _ = box_size[0]
        except TypeError:
            box_size = [box_size]
        if radius is None: 
            radius = box_size[0] / 2 - width
        r = TestImage._radius(box_size, dtype = np.float64)
        arr = r.copy()
        if radius < 0:
            arr = -radius - arr
            if width > -radius: width = -radius
            arr.data[arr.data<0] = 0
        else:
            arr = arr - radius
            arr.data[arr.data<0] = 0
        weight = 1/float(width)
        arr = np.exp(-0.5 * arr.data * arr.data * weight)
        if radius > 0:
            arr[r.data>(radius + width)] = 0
        else:
            arr[r.data < (-radius - width)] = 0
        return EMImage(arr)
    
    
    @staticmethod
    def index(box_size, dtype = np.int32, center = None):
        return TestImage._index(box_size, dtype, center)
    
    @staticmethod
    def _index(box_size, dtype = np.int32, center = None):
        
        if isinstance(box_size, int):
            box_size = [box_size]
        if center is None: center = [b / 2 for b in box_size]
        if len(box_size) == 1: 
            arr = np.mgrid[0:box_size[0]] - center[0]
        elif len(box_size) == 2:
            arr = np.mgrid[0:box_size[0], 0:box_size[1]]
            arr[0] -= center[0]
            arr[1] -= center[1]
        elif len(box_size) == 3:
            arr = np.mgrid[0:box_size[0], 0:box_size[1], 0:box_size[2]]
            arr[0] -= center[0]
            arr[1] -= center[1]      
            arr[2] -= center[2] 
        return arr
    
    @staticmethod
    def radius(box_size, dtype = np.float32):
        return TestImage._radius(box_size, dtype)
    
    @staticmethod
    def _radius(box_size, dtype = np.int32):
        ''' Return a radius array with center at box_size/2
        @param box_size: The box size array. 
        '''
        arr = TestImage._index(box_size, dtype)
        arr = np.sqrt(np.sum([x*x for x in arr], axis = 0))
        if dtype == np.float64:
            return EMImage(arr)
        else:
            return EMImage(arr.astype(np.float32))

EMImage = _EMImage

if __name__ == "__main__":
    pass
