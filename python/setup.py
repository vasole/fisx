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
DATA_DIR = os.path.join('fisx', 'fisx_data')
DOC_DIR = os.path.join('fisx', 'fisx_data')

from distutils.command.build_py import build_py
class smart_build_py(build_py):
    def run (self):
        toReturn = build_py.run(self)
        global DATA_DIR
        global DOC_DIR
        global INSTALL_DIR
        defaultDataPath = os.path.join('fisx', 'fisx_data')
        defaultDocPath = os.path.join('fisx', 'fisx_data')
        if (DATA_DIR == defaultDataPath) or\
           (DOC_DIR == defaultDocPath):
            #default, just make sure the complete path is there
            install_cmd = self.get_finalized_command('install')
            INSTALL_DIR = getattr(install_cmd, 'install_lib')

        #packager should have given the complete path
        #in other cases
        if DATA_DIR == defaultDataPath:
            DATA_DIR = os.path.join(INSTALL_DIR, DATA_DIR)
        if DOC_DIR == defaultDocPath:
            #default, just make sure the complete path is there
            DOC_DIR = os.path.join(INSTALL_DIR, DOC_DIR)
        target = os.path.join(self.build_lib, "fisx", "DataDir.py")
        fid = open(target,'r')
        content = fid.readlines()
        fid.close()
        fid = open(target,'w')
        for line in content:
            lineToBeWritten = line
            if lineToBeWritten.startswith("DATA_DIR"):
                lineToBeWritten = "DATA_DIR = r'%s'\n" % DATA_DIR
            if line.startswith("DOC_DIR"):
                lineToBeWritten = "DOC_DIR = r'%s'\n" % DOC_DIR
            fid.write(lineToBeWritten)
        fid.close()
        return toReturn

from distutils.command.install_data import install_data
class smart_install_data(install_data):
    def run(self):
        global INSTALL_DIR
        global DATA_DIR
        global DOC_DIR
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
data_files = [(DATA_DIR, fileList)]

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

packages = ['fisx', 'fisx.tests']
setup(
    name='fisx',
    author="V. Armando Sole",
    author_email="sole@esrf.fr",
    license="MIT",
    packages=packages,
    ext_modules=ext_modules,
    data_files=data_files,
    cmdclass=cmdclass,
)
