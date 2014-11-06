#!/usr/bin/env python
import glob
import os
import sys
from distutils.core import setup, Extension
from distutils.sysconfig import get_python_lib
try:
    from Cython.Distutils import build_ext
except:
    build_ext = None
import numpy
# deal with required data

#for the time being there is no doc directory
FISX_DATA_DIR = os.getenv("FISX_DATA_DIR")
FISX_DOC_DIR = os.getenv("FISX_DOC_DIR")

if FISX_DATA_DIR is None:
    FISX_DATA_DIR = os.path.join('fisx', 'fisx_data')
if FISX_DOC_DIR is None:
    FISX_DOC_DIR = os.path.join('fisx', 'fisx_data')

def getFisxVersion():
    cppDir = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "src")
    content = open(os.path.join(cppDir, "fisx_version.h"), "r").readlines()
    for line in content:
        if "FISX_VERSION_STR" in line:
            version = line.split("FISX_VERSION_STR")[-1].replace("\n","")
            version = version.replace(" ","")
            return version[1:-1]

__version__ = getFisxVersion()

print("fisx X-Ray Fluorescence Toolkit %s\n" % __version__)

from distutils.command.build_py import build_py
class smart_build_py(build_py):
    def run (self):
        toReturn = build_py.run(self)
        global FISX_DATA_DIR
        global FISX_DOC_DIR
        global INSTALL_DIR
        defaultDataPath = os.path.join('fisx', 'fisx_data')
        defaultDocPath = os.path.join('fisx', 'fisx_data')
        if (FISX_DATA_DIR == defaultDataPath) or\
           (FISX_DOC_DIR == defaultDocPath):
            #default, just make sure the complete path is there
            install_cmd = self.get_finalized_command('install')
            INSTALL_DIR = getattr(install_cmd, 'install_lib')

        #packager should have given the complete path
        #in other cases
        if FISX_DATA_DIR == defaultDataPath:
            FISX_DATA_DIR = os.path.join(INSTALL_DIR, FISX_DATA_DIR)
        if FISX_DOC_DIR == defaultDocPath:
            #default, just make sure the complete path is there
            FISX_DOC_DIR = os.path.join(INSTALL_DIR, FISX_DOC_DIR)
        target = os.path.join(self.build_lib, "fisx", "DataDir.py")
        fid = open(target,'r')
        content = fid.readlines()
        fid.close()
        fid = open(target,'w')
        for line in content:
            lineToBeWritten = line
            if lineToBeWritten.startswith("DATA_DIR"):
                lineToBeWritten = "DATA_DIR = r'%s'\n" % FISX_DATA_DIR
            if line.startswith("FISX_DOC_DIR"):
                lineToBeWritten = "FISX_DOC_DIR = r'%s'\n" % FISX_DOC_DIR
            fid.write(lineToBeWritten)
        fid.close()
        return toReturn

from distutils.command.install_data import install_data
class smart_install_data(install_data):
    def run(self):
        global INSTALL_DIR
        global FISX_DATA_DIR
        global FISX_DOC_DIR
        #need to change self.install_dir to the library dir
        install_cmd = self.get_finalized_command('install')
        self.install_dir = getattr(install_cmd, 'install_lib')
        INSTALL_DIR = self.install_dir
        print("fisx to be installed in %s" %  self.install_dir)
        return install_data.run(self)

topLevel = os.path.dirname(os.getcwd())
fileList = glob.glob(os.path.join(topLevel, "fisx_data", "*.dat"))
fileList.append(os.path.join(topLevel, "LICENSE"))
fileList.append(os.path.join(topLevel, "README.md"))
data_files = [(FISX_DATA_DIR, fileList)]

# actual build stuff
FORCE = False
if build_ext:
    #make sure everything is totally clean?
    cython_dir = os.path.join(os.getcwd(), "cython")
    fileList  = glob.glob(os.path.join(cython_dir, "*.cpp"))
    fileList += glob.glob(os.path.join(cython_dir, "*.h"))
    if FORCE:
        for fname in fileList:
            if os.path.exists(fname):
                os.remove(fname)
    #this does not work:
    #src = glob.glob(os.path.join(cython_dir, "*pyx"))
    multiple_pyx = os.path.join(cython_dir, "_fisx.pyx")
    if os.path.exists(multiple_pyx):
        try:
            os.remove(multiple_pyx)
        except:
            print("WARNING: Could not delete file. Assuming up-to-date.")
    if not os.path.exists(multiple_pyx):
        pyx = glob.glob(os.path.join(cython_dir, "*pyx"))
        f = open(multiple_pyx, 'wb')
        for fname in pyx:
            inFile = open(fname, 'rb')
            lines = inFile.readlines()
            inFile.close()
            for line in lines:
                f.write(line)
        f.close()
    src = [multiple_pyx]
else:
    src = glob.glob(os.path.join(cython_dir,'*.cpp'))
src += glob.glob(os.path.join(os.path.dirname(os.getcwd()),'src', 'fisx_*.cpp'))

include_dirs = [numpy.get_include(),
                os.path.join(os.path.dirname(os.getcwd()), "src")]

if sys.platform == 'win32':
    extra_compile_args = ['/EHsc']
    extra_link_args = []
else:
    extra_compile_args = []
    extra_link_args = []

def buildExtension():
    module = Extension(name="fisx._fisx",
                    sources=src,
                    include_dirs=include_dirs,
                    extra_compile_args=extra_compile_args,
                    extra_link_args=extra_link_args,
                    language="c++",
                    )
    return module

ext_modules = [buildExtension()]

cmdclass = {'install_data':smart_install_data,
            'build_py':smart_build_py,
            'build_ext': build_ext,
            }

description = "Quantitative X-Ray Fluorescence Analysis Support Library"
long_description = """
Tools to evaluate the expected X-ray fluorescence measured when a sample is excitated by an X-ray beam. Secondary and tertiary excitation effects taken into account.
"""

packages = ['fisx', 'fisx.tests']
setup(
    name='fisx',
    version=__version__,
    author="V. Armando Sole",
    author_email="sole@esrf.fr",
    description=description,
    long_description=long_description,
    license="MIT",
    url="https://github.com/vasole/fisx",
    packages=packages,
    ext_modules=ext_modules,
    data_files=data_files,
    cmdclass=cmdclass,
)
