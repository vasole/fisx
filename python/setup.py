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

packages = ['fisx']
setup(
    name='fisx',
    author="V. Armando Sole",
    license="MIT",
    packages=packages,
    ext_modules=ext_modules,
    cmdclass={'build_ext': build_ext},
)
