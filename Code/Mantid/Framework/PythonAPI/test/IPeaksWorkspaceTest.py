import unittest

from MantidFramework import mtd
mtd.initialise()
from mantidsimple import *

class IPeaksWorkspaceTest(unittest.TestCase):
    """
    Test the python interface to PeaksWorkspace's
    """
    
    def setUp(self):
        LoadEventNexus(Filename='CNCS_7860_event.nxs', OutputWorkspace='cncs')
        CreatePeaksWorkspace(InstrumentWorkspace='cncs', OutputWorkspace='peaks')
        pass
    
    def test_interface(self):
        """ Rudimentary test to get peak and get/set some values """ 
        pws = mtd['peaks']
        self.assertEqual(pws.getNumberPeaks(), 1)
        p = pws.getPeak(0)
        
        # Try a few IPeak get/setters. Not everything.
        p.setH(234)
        self.assertEqual(p.getH(), 234)
        p.setHKL(5,6,7)
        self.assertEqual(p.getH(), 5)
        self.assertEqual(p.getK(), 6)
        self.assertEqual(p.getL(), 7)
        
        hkl = p.getHKL()
        assert hkl == V3D(5,6,7)
        
        p.setIntensity(456)
        p.setSigmaIntensity(789)
        self.assertEqual(p.getIntensity(), 456)
        self.assertEqual(p.getSigmaIntensity(), 789)
        
        # Finally try to remove a peak
        pws.removePeak(0)
        self.assertEqual(pws.getNumberPeaks(), 0)
        
        # Create a new peak at some Q in the lab frame
        qlab = V3D(1,2,3)
        p = pws.createPeak(qlab, 1.54)
        self.assertAlmostEquals( p.getQLabFrame().getX(), 1.0, 3)
        
        # Now try to add the peak back
        pws.addPeak(p)
        self.assertEqual(pws.getNumberPeaks(), 1)
        
        # Check that it is what we added to it
        p = pws.getPeak(0)
        self.assertAlmostEquals( p.getQLabFrame().getX(), 1.0, 3)
                
        
if __name__ == '__main__':
    unittest.main()

    
