from pf_refinement import InfoFile, StarFile, EMImage
from gui import gdGui
import os

class MaskedClassification(InfoFile):
    def __init__(self):
        super(MaskedClassification, self).__init__()
        
    def __call__(self):
        self.getVals()
        self.makeMask()
        self.relionSubtract()
        
    def getVals(self):
        vals={
            'input_star': 'XXXSTARFILEXXX',
            'protofilament_volume': 'XXXMRCFILEXXX',
            'protofilament_mask': 'XXXMRCFILEXXX',
            'ROI_mask': 'XXXMRCFILEXXX'}
        
        temp=gdGui('Focused Classification', **vals)
        self.vals=temp.sendValues()
        
    def makeMask(self):
        proto=EMImage(self.vals['protofilament_volume'])
        proto_mask=EMImage(self.vals['protofilament_mask'])
        roi=EMImage(self.vals['ROI_mask'])
        
        roi.mult(-1)
        roi.add(1)
        roi.mult(proto_mask)
        roi.mult(proto)
        
        roi.write_mrc('ROI_subtraction_volume.mrc')
        
    def relionSubtract(self):
        command='relion_project --i ROI_subtraction_volume.mrc '\
        '--o roi_for_classification --ctf --angpix %g --ang %s '\
        '--subtract_exp'%(self.pixel_size, self.vals['input_star'])
        check=os.system(command)
        if check!=0:
            raise RuntimeError('RELION did not run properly. Try running the '
                            'following command to troubleshoot \n\n'
                            '%s \n'%command)