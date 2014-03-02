import unittest
import os
import sys
import glob
import unittest

def getSuite(auto=True):
    pythonFiles = glob.glob(os.path.join(os.path.dirname(__file__), "*.py"))
    sys.path.insert(0, os.path.dirname(__file__))
    testSuite = unittest.TestSuite()
    for fname in pythonFiles:
        if os.path.basename(fname) in ["__init__.py", "testAll.py"]:
            continue
        modName = os.path.splitext(os.path.basename(fname))[0]
        try:
            module = __import__(modName)
        except ImportError:
            print("Failed to import %s" % fname)
            continue
        if hasattr(module, "getSuite"):
            testSuite.addTest(module.getSuite(auto))
    return testSuite

def main(auto=True):
    unittest.TextTestRunner(verbosity=2).run(getSuite(auto=auto))

if __name__ == '__main__':
    if len(sys.argv) > 1:
        auto = False
    else:
        auto = True
    main(auto)
