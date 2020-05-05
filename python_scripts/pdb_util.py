'''
    PDB file manipulation
    usage:

    p = PDB()
    p.read_file(pdbfilename)
    p1 = p.symmetrize([[1,0,0], [0, 1, 0], [1, 0, 1]])
    cent_x, cent_y, cent_z = p1.get_center()
    p1.to_center()
    vol = p1.to_volume(pixel_size = 1.3, box_size = 80) # EMImage
    vol.write_image('vol.mrc')
    p1.write("out.pdb")  # write to file
'''
from collections import OrderedDict
import cupy as np

# The common used element in pdb
ELEMENT_DEFS = {'H' : (1.0, 1.00794),
              'C': (6.0, 12.0107),
              'N': (7.0, 14.00674),
              'O': (8.0, 15.9994),
              'P': (15.0, 30.973761),
              'S': (16.0, 32.066),
              'W': (18.0, 1.00794*2.0+15.9994),
              'AU': (79.0, 196.96655),
              'U': (0, 0)}

ATOM_ALIAS = {"NZ" : "N",
              "CB" : "C",
              "CD" : "C",
              "CE" : "C",
              "OG1" : "O",
              "OG2" : "O"
    }
class PDB(object):
    """ PDB
    """
    def __init__(self):
        '''
            one pdb contains multiple models
            one model is dictionary of chains
            one chain is list of atoms. atom.chainID should be same as chain
        '''
        self.models = []

    def read_file(self, filename):
        """ Read a PDB file

        Args:
            filename: A `str`. The PDB file to read.
                The filename usually ends with .pdb.
        Return:
            self. The PDB instance
        """
        self._init()
        with open(filename) as pdb_f:
            while True:
                line = pdb_f.readline()
                if len(line) == 0: # End of file
                    break
                if line.startswith("ATOM"):
                    if not self.models: # This means no MODEL line
                        self.models.append(OrderedDict())
                    atom = Atom.parse(line)
                    # append the atom to chain
                    if atom.chainID not in self.models[-1]:
                        self.models[-1][atom.chainID] = []
                    self.models[-1][atom.chainID].append(atom)
                elif line.startswith("MODEL"):
                    self.models.append(OrderedDict())
        # in case there is only MODEL line with no atom
        if self.models and len(self.models[-1]) == 0:
            self.models.pop()

    def write(self, filename):
        """ Write pdb to file

        Args:
            filename: A `str`. The output file name
        """
        with open(filename, 'w') as pdb_f:
            for atom in self.next_atom():
                pdb_f.write("%s\n" % atom)

    def _init(self):
        """ Initliaze the pdb instance. eg. Clear
        """
        self.models = []

    def copy(self):
        '''return a deep copy of pdb

        Returns:
            A `PDB` instance. Copied pdb from this pdb.
        '''
        p_copy = PDB()
        p_copy.models = []
        for model in self.models:
            new_model = OrderedDict()
            for chain_id in model:
                new_model[chain_id] = []
                for atom in model[chain_id]:
                    new_model[chain_id].append(atom.copy())
            p_copy.models.append(new_model)
        return p_copy

    def symmetrize(self, symm_a, inplace=False):
        '''Symmetrize the pdb according to the transformation matrix.

            helical symmetry is:
            cosa, sina, 0, 0
            -sina, cosa, 0, 0
            0, 0, 1, dz
        Args:
            symm_a: A `ndarray`. The transformation matrix. Should be 3x4 matrix.
                first 3x3 is rotation, last column is translation
            inplace: A `bool`. Inplace symmetrize if true.
        Return:
            The symmetrized pdb.
        '''
        # Check the symmetry matrix:
        if not isinstance(symm_a, np.ndarray):
            symm_a = np.array(symm_a)
        if symm_a.shape == (3, 4):
            symm_a = symm_a.T
        if symm_a.shape != (4, 3):
            raise Exception("The input matrix should be 3x4 transform matrix")
        if inplace:
            rpdb = self
        else:
            rpdb = self.copy()
        for model_ in rpdb.models:
            coords = []
            # get the coordinates
            coords = [[atom.x, atom.y, atom.z, 1] for atom in rpdb.next_atom_in_model(model_)]
            arr = np.array(coords)
            arr_t = np.dot(arr, symm_a)

            for idx, atom in enumerate(rpdb.next_atom_in_model(model_)):
                atom.x = arr_t[idx][0]
                atom.y = arr_t[idx][1]
                atom.z = arr_t[idx][2]

        return rpdb
#
#     def to_volume1(self,
#                    voxel_size=1.,
#                    resolution=None,
#                    box_size=64):
#         '''Depreated.
#             convert pdb to volume.
#             If want the volume be centered, center the pdb first use pdb.to_center()
#
#         Args:
#             voxel_size: A `float`. The voxel size for the volume
#             resolution: A `float`. The desired resolution to filter
#             box_size: A `int`. Even number. The output volume size.
#         '''
#         if resolution is None or resolution < voxel_size * 2:
#             resolution = voxel_size * 2
#         if isinstance(box_size, int): box_size = [box_size] * 3
#         if len(box_size) != 3:
#             raise Exception("Wrong format of boxSize")
#         from cryofilia.EMImage import EMImage, TestImage
#         outmap = EMImage(*box_size)
#         xt, yt, zt = map(lambda x: x / 2, box_size)
#         gaus = TestImage.sphere(64, 12)
# #        gaus.process_inplace("mask.gaussian",{"outer_radius":12.0})
# 
#         for atom in self.next_atom():
#             x, y, z = atom.x, atom.y, atom.z
#             x, y, z = x / voxel_size, y / voxel_size, z / voxel_size
#             ix, iy, iz = int(x), int(y), int(z)
#             ax, ay, az = x - ix, y - iy, z - iz
#             idxx, idxy, idxz = np.mgrid[-1:2, -1:2, -1:2]
#             at = (idxx-(2*idxx-1)*ax) * (idxy-2*(idxy-1)*ay) * (idxz-(2*idxz-1)*az)
#             at *= atom.weight
#
#             outmap.insert_sum(gaus, (ix+xt, iy+yt, iz+zt))
#
#         outmap.voxelSize = voxel_size
#         return outmap

    def to_volume(self, voxel_size=1., box_size=64):
        ''' Convert pdb to volume. using linear interpolation
            (0,0,0) in pdb is (box_size/2, box_size/2, box_size/2) in volume

        Args:
            voxel_size: A `float`. The voxel size of output volume
            box_size: A `int`. Even number. The size of output volume
        '''
        from cryofilia.EMImage import EMImage
        cent_x, cent_y, cent_z = box_size / 2, box_size / 2, box_size / 2
        x, y, z = [], [], []
        weights = []
        vol = np.zeros([box_size, box_size, box_size])
        for atom in self.next_atom():
            x.append(atom.x/voxel_size)
            y.append(atom.y/voxel_size)
            z.append(atom.z/voxel_size)
            weights.append(atom.weight)
        x, y, z = np.array(x), np.array(y), np.array(z)
        weights = np.array(weights)
        x0 = np.floor(x).astype(np.int32)
        y0 = np.floor(y).astype(np.int32)
        z0 = np.floor(z).astype(np.int32)
        idxes = np.logical_and.reduce([z0 >= -cent_x,
                                      y0 >= -cent_x,
                                      x0 >= -cent_x,
                                      z0 < cent_x - 1,
                                      y0 < cent_x  - 1,
                                      x0 < cent_x - 1])
        x, y, z = x[idxes], y[idxes], z[idxes]
        x0, y0, z0 = x0[idxes], y0[idxes], z0[idxes]
        weights = weights[idxes]
        f_x, f_y, f_z = x-x0, y-y0, z-z0
        mf_x, mf_y, mf_z = 1. - f_x, 1. - f_y, 1. - f_z

        d000 = mf_z * mf_y
        d001 = d000 * f_x
        d000 *= mf_x
        d010 = mf_z * f_y
        d011 = d010 * f_x
        d010 *= mf_x
        d100 = f_z * mf_y
        d101 = d100 * f_x
        d100 *= mf_x
        d110 = f_z * f_y
        d111 = d110 * f_x
        d110 *= mf_x
        x0 += cent_x
        y0 += cent_y
        z0 += cent_z
        x_i1, y_i1, z_i1 = x0+1, y0+1, z0+1

        vol[z0, y0, x0] += d000 * weights
        vol[z0, y0, x_i1] += d001 * weights
        vol[z0, y_i1, x0] += d010 * weights
        vol[z0, y_i1, x_i1] += d011 * weights
        vol[z_i1, y0, x0] += d100 * weights
        vol[z_i1, y0, x_i1] += d101 * weights
        vol[z_i1, y_i1, x0] += d110 * weights
        vol[z_i1, y_i1, x_i1] += d111 * weights

        img = EMImage(vol)
        img.voxelSize = voxel_size
        return img

    def transform(self, shx, shy, shz, inplace=True):
        """ Transform the pdb.
        """
        if inplace:
            result_pdb = self
        else:
            result_pdb = self.copy()
        for atom in result_pdb.next_atom():
            atom.x = atom.x + shx
            atom.y = atom.y + shy
            atom.z = atom.z + shz
        return result_pdb

    def rotate(self, phi, theta, psi, inplace=True):
        ''' Rotate pdb according to cryoem convention
            phi, theta, psi in degree

        Args:
            phi: A `float`. in degree. First rotation around z, CC positive
            theta: A `float`. in degree. Second rotation around Y, CC positive
            psi: A `float`. in degree. Third rotation around z, CC positive
            inplace: A `bool`. Modify inplace if true.

        Return:
            The rotated pdb. If inplace is true, self is returned.
        '''
        from math import cos, sin, pi
        phi = phi * pi / 180
        theta = theta * pi / 180
        psi = psi * pi / 180

        cos_x = cos(psi)
        sin_x = sin(psi)
        cos_t = cos(theta)
        sin_t = sin(theta)
        cos_p = cos(phi)
        sin_p = sin(phi)
        sym = np.array([[cos_x*cos_t*cos_p - sin_x*sin_p,
                         cos_x*cos_t*sin_p+sin_x*cos_p,
                         -cos_x*sin_t, 0],
                        [-sin_x*cos_t*cos_p - cos_x*sin_p,
                         -sin_x*cos_t*sin_p + cos_x*cos_p,
                         sin_x*sin_t, 0],
                        [-sin_t*cos_p, sin_t*sin_p, cos_t, 0]])
        return self.symmetrize(sym, inplace=inplace)

    def rename_chain(self, old_chain, new_chain):
        """ Rename pdb chain name from old_chain to new_chain

        Args:
            old_chain: A `str`. The old chain name
            new_chain: A `str`. The new chain name
        """
        if old_chain == new_chain:
            return
        for _model in self.models:
            new_mol = []
            if new_chain in _model:
                new_mol = _model[new_chain]
            if old_chain in _model:
                for atom in _model[old_chain]:
                    atom.chainID = new_chain
                    new_mol.append(atom)
                del _model[old_chain]
            if new_mol:
                _model[new_chain] = new_mol

    def extract_chains(self, chains):
        """ Extract chains to new pdb

        Args:
            chains: A `list` or `tuple`. The list of chains to extract
        Returns:
            A `PDB` instance. The extracted pdb
        """
        result_pdb = PDB()
        for _model in self.models:
            new_model = OrderedDict()
            for chain in chains:
                if chain in _model:
                    new_model[chain] = _model[chain][:]
            if new_model:
                result_pdb.models.append(new_model)
        return result_pdb

    def to_center(self, gravity=False):
        """ Center the pdb according to coordinate or gravity

        Args:
            gravity: A `bool`. Use gravity instead of center of coordinates.
        """
        cent_x, cent_y, cent_z = self.calc_center(gravity=gravity)
        self.transform(-cent_x, -cent_y, -cent_z)

    def next_atom_in_model(self, model):
        """ Generate to get next atom in model

        Args:
            model: A `OrderedDict`. The model to get the atoms from
        Returns:
            Generator to loop all atoms in model
        """
        for chain_id in model:
            model_ = model[chain_id]
            for atom in model_:
                yield atom

    def next_atom(self):
        """ Generator to get the next atom

        Returns:
            Generator to loop all atoms in PDB
        """
        for model_ in self.models:
            for chain_id in model_:
                model = model_[chain_id]
                for atom in model:
                    yield atom

    def calc_center(self, gravity=False):
        """ Calculate the center of pdb.

        Args:
            gravity: A `bool`. Calculate center of gravity instead of center of coordinate
        Returns:
            (cent_x, cent_y, cent_z).
        """
        x_sum, y_sum, z_sum = 0, 0, 0
        n_atoms = 0
        for atom in self.next_atom():
            x, y, z = atom.x, atom.y, atom.z
            if gravity:
                weight = atom.weight
                n_atoms += atom.weight
            else:
                weight = 1
                n_atoms += 1
            x_sum, y_sum, z_sum = x_sum + x*weight, y_sum + y*weight, z_sum + z*weight
        return x_sum / n_atoms, y_sum / n_atoms, z_sum / n_atoms

    @staticmethod
    def get_center(pdb, gravity=False):
        ''' Static method to get the center of pdb
        '''
        if isinstance(pdb, PDB):
            return pdb.calc_center(gravity)
        if isinstance(pdb, str):
            temp_p = PDB()
            temp_p.read_file(pdb)
            return temp_p.calc_center(gravity)


    @property
    def center(self):
        """ Property center.
        """
        return self.calc_center()

    @property
    def cog(self):
        """ Property center of gravity
        """
        return self.calc_center(gravity=True)

    def __add__(self, other):
        result = self.copy()
        result.models.extend(other.models)
        return result
    def __radd__(self, other):
        if other == 0:
            return self.copy()
        return other.__iadd__(self)
    def __iadd__(self, other):
        self.models.extend(other.models)
        return self

    def __str__(self):
        return "PDB with %d models" % (len(self.models))



class ElementNotFound(Exception):
    """ Exception for Elements
    """
    pass

class Element(object):
    """ Element definition
    """
    def __init__(self,
                 atom_name,
                 atom_index=None,
                 atomic_mass=None):
        """
        Args:
            atom_name: A `str`. The name of the element
            atom_index: A `str`. The index of the element in period table.
            atomic_mass: A `str`. The atomic mass of the element
        Raises:
            ElementNotFound. If atom_name not in ELEMENT_DEFS and atom_index not provided.
        """
        self._name = atom_name
        self.atom_index = atom_index
        self.weight = atomic_mass
        if self.atom_index is None:
            if atom_name not in ELEMENT_DEFS:
                raise ElementNotFound("Atom index not provided")
            self.atom_index = ELEMENT_DEFS[atom_name][0] or 0
        if self.weight is None:
            if atom_name not in ELEMENT_DEFS:
                raise ElementNotFound("Atom weight not provided")
            self.weight = ELEMENT_DEFS[atom_name][1] or 0
        if self.weight < 0:
            raise ElementNotFound("Atomic mass should larger than 0")
        if self.atom_index < 0:
            raise ElementNotFound("Atom index should be larger than 0")
    @property
    def name(self):
        """ Property name
        """
        return self._name

    @name.setter
    def name(self, atom_name):
        """ Property name setter
        """
        self._name = atom_name
        if self._name in ELEMENT_DEFS:
            self.atom_index = ELEMENT_DEFS[atom_name][0]
            self.weight = ELEMENT_DEFS[atom_name][1]

class Atom(object):
    '''http://www.wwpdb.org/documentation/file-format-content/format33/sect9.html#ATOM

    The name used here is to be consistent with reference. So the name is CamelCase
    '''
    field_width = (6, 5, -1, 4, 1, 3, -1, 1, 4, 1, -3, 8, 8, 8, 6, 6, -10, 2, 2)
    field_name = ('sig', 'serial', 'name', 'altLoc', 'resName', 'chainID',
                 'resSeq', 'iCode', 'x', 'y', 'z', 'occupancy', 'tempFactor',
                 'element', 'charge')
    field_type = (str, int, str, str, str, str,
                 int, str, float, float, float, float, float,
                 str, str)
    input_format = ['%s%s' % (abs(_), 'x' if _ < 0 else 's') for _ in field_width]
#     code = {str:'s', int:'d', float:'f'}
#     struct_fmt = "".join(['%s%s' % (abs(fw), 'x' if fw < 0 else 's')
#                           for fw in field_width])
#     output_format = ""
#     j = 0
#     for i in range(len(field_width)):
#         current_width = field_width[i]
#         if current_width > 0:
#             output_format += "%"
#             if field_name[j] in ['sig']:
#                 output_format += "-"
#             output_format += str(current_width)
#             if field_name[j] in ['x', 'y', 'z']:
#                 output_format += '.3'
#             elif field_name[j] in ['tempFactor', 'occupancy']:
#                 output_format += '.2'
#
#             output_format += code[field_type[j]]
#             j += 1
#         else:
#             output_format += " " * (-current_width)
    struct_fmt = "6s5s1x4s1s3s1x1c4s1s3x8s8s8s6s6s10x2s2s"
    output_format = "%-6s%5d %4s%1s%3s %1s%4d%1s   %8.3f%8.3f%8.3f%6.2f%6.2f          %2s%2s"

    def __init__(self,
                 atom_name="",
                 x=0,
                 y=0,
                 z=0,
                 atom_index=None,
                 atom_weight=None,
                 index=-1):
        """ Create a new Atom instance

        Args:
            atom_name: The name of the atom. Should be an element name.
            x: A `float`. The x coordinate of the atom
            y: A `float`. The y coordinate of the atom
            z: A `float`. The z coordinate of the atom
            atom_index: A `int`. The index of the atom in period table
            atom_weight: A `float`. The atomic mass of the element
            index: A `int`. The index of the atom in pdb
        """
        for idx in range(len(Atom.field_name)):
            setattr(self, Atom.field_name[idx], Atom.field_type[idx]())
        if atom_name:
            self._element = Element(atom_name, atom_index, atom_weight)
        else:
            self._element = None

        self.x = x
        self.y = y
        self.z = z
        self.occupancy = 1.
        self.serial = 0
        self.index = index
        self._name = atom_name
        self.sig = "ATOM"

    def copy(self):
        """ Copy to another Atom instance

        Returns:
            A `Atom`. Copied from this instance.
        """
        atom = Atom()
        for field in Atom.field_name:
            setattr(atom, field, getattr(self, field))
        atom.index = self.index
        return atom

    @staticmethod
    def parse(line):
        """Parse a line in pdb and return the `Atom` instance

        Args:
            line: A `str`. The pdb line to parse
        """
        import struct
        atom = Atom()
        line = line.strip("\n")
        struct_size = struct.calcsize(Atom.struct_fmt)
        # add extra space for non standard format
        if len(line) < struct_size:
            line += " " * (struct_size - len(line))
        result = struct.unpack(Atom.struct_fmt, line[:struct_size])
        for fname, ftype, value in zip(Atom.field_name, Atom.field_type, result):
            setattr(atom, fname, ftype(value))
        return atom

    def __str__(self):
        result = Atom.output_format % tuple([getattr(self, x) for x in Atom.field_name])
        return result

    def __repr__(self):
        return str(self)

    @property
    def element(self):
        """ property element name
        """
        return self._element.name

    @element.setter
    def element(self, value):
        """ setter for element name
        """
        if isinstance(value, Element):
            self._element = value
        elif isinstance(value, str):
            value = value.strip()
            if value in ELEMENT_DEFS:
                self._element = Element(value)
            elif value in ATOM_ALIAS:
                self._element = Element(ATOM_ALIAS[value])
            elif value != "":
                raise ElementNotFound("Atom %s not found in definition" % value)

    @property
    def name(self):
        """ property name
        """
        return self._name

    @name.setter
    def name(self, value):
        """ property name setter
        """
        self._name = value.strip()

    @property
    def index(self):
        """ property index
        """
        return self.serial
    @index.setter
    def index(self, value):
        """ property index setter
        """
        self.serial = value

    @property
    def weight(self):
        """ property weight
        """
        return self._element.weight
