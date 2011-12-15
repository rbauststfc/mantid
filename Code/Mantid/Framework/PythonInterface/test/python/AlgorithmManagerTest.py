import unittest

from mantid.api import AlgorithmManager
from mantid.api import (IAlgorithm, Algorithm, AlgorithmProxy, 
                        registerAlgorithm)
import sys

class IsAnAlgorithm(Algorithm):
    def init_(self):
        pass
    
class NotAnAlgorithm(object):
    pass

class AlgorithmManagerTest(unittest.TestCase):
    
    def test_create_default_version(self):
        try:
            alg = AlgorithmManager.Instance().create("ConvertUnits")
        except RuntimeError:
            self.fail(str(exc))
        # Tests
        self.assertNotEqual(alg, None)
        self.assertEquals(alg.name(), "ConvertUnits")
        self.assertEquals(alg.version(), 1)
        self.assertEquals(alg.category(), "Transforms\\Units")
        
    def test_create_unknown_alg_throws(self):
        self.assertRaises(RuntimeError, AlgorithmManager.Instance().create,"DoesNotExist")
        
    def test_created_alg_isinstance_of_IAlgorithm(self):
        alg = AlgorithmManager.Instance().create("ConvertUnits")
        self.assertTrue(isinstance(alg, IAlgorithm))
        
    def test_managed_cppalg_isinstance_of_AlgorithmProxy(self):
        alg = AlgorithmManager.Instance().create("ConvertUnits")
        self.assertTrue(isinstance(alg, AlgorithmProxy))

    def test_unmanaged_cppalg_isinstance_of_Algorithm(self):
        alg = AlgorithmManager.Instance().createUnmanaged("ConvertUnits")
        self.assertTrue(isinstance(alg, Algorithm))
        
    def test_pyalg_isinstance_of_Algorithm(self):
        alg = IsAnAlgorithm()
        self.assertTrue(isinstance(alg, Algorithm))
        self.assertTrue(isinstance(alg, IAlgorithm))
        
    def test_algorithm_registration_with_valid_object_succeeds(self):
        try:
            registerAlgorithm(IsAnAlgorithm)
            noerror = True
        except Exception:
            noerror = False
        self.assertTrue(noerror)

    def test_algorithm_registration_with_invalid_object_throws(self):
        try:
            registerAlgorithm(NotAnAlgorithm)
            error = False
        except ValueError:
            error = True
        self.assertTrue(error)

if __name__ == '__main__':
    unittest.main()        
