#!/usr/bin/env python
# -*- coding: utf-8 -*-
#/*##########################################################################
#
# The fisx library for X-Ray Fluorescence
#
# Copyright (c) 2014-2018 European Synchrotron Radiation Facility
#
# This file is part of the fisx X-ray developed by V.A. Sole
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in
# all copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
# THE SOFTWARE.
#
#############################################################################*/
import glob
import os
import sys
if "bdist_wheel" in sys.argv:
    from setuptools import setup, Extension
    from setuptools.command.build_py import build_py
    from distutils.command.install_data import install_data
    from setuptools.command.sdist import sdist
else:
    try:
        from setuptools import setup, Extension
        from setuptools.command.build_py import build_py
        from distutils.command.install_data import install_data
        from setuptools.command.sdist import sdist
    except ImportError:
        from distutils.core import setup, Extension
        from distutils.sysconfig import get_python_lib
        from distutils.command.build_py import build_py
        from distutils.command.install_data import install_data
        from distutils.command.sdist import sdist


# check if cython is not to be used despite being present
def use_cython():
    """
    Check if cython is disabled from the command line or the environment.
    """
    if "WITH_CYTHON" in os.environ:
        if os.environ["WITH_CYTHON"] in ["False", "0", 0]:
            print("No Cython requested by environment")
            return False

    if ("--no-cython" in sys.argv):
        sys.argv.remove("--no-cython")
        os.environ["WITH_CYTHON"] = "False"
        print("No Cython requested by command line")
        return False
    return True

if use_cython():
    try:
        from Cython.Distutils import build_ext
        from Cython.Compiler.Version import version
        if (version < "0.28.3") and (sys.version_info > (3, 7)):
            build_ext = None
        elif version < "0.17":
            build_ext = None
    except:
        build_ext = None
else:
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


def get_version():
    """Returns current version number from version.py file"""
    import version
    return version.strictversion

__version__ = get_version()

print("fisx X-Ray Fluorescence Toolkit %s\n" % __version__)

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
            if lineToBeWritten.startswith("FISX_DATA_DIR"):
                lineToBeWritten = "FISX_DATA_DIR = r'%s'\n" % FISX_DATA_DIR
            if line.startswith("FISX_DOC_DIR"):
                lineToBeWritten = "FISX_DOC_DIR = r'%s'\n" % FISX_DOC_DIR
            fid.write(lineToBeWritten)
        fid.close()
        return toReturn

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

topLevel = os.path.dirname(__file__)
if not len(topLevel):
    topLevel = "."
fileList = glob.glob(os.path.join(topLevel, "fisx_data", "*.dat"))
fileList.append(os.path.join(topLevel, "changelog.txt"))
fileList.append(os.path.join(topLevel, "LICENSE"))
fileList.append(os.path.join(topLevel, "README.rst"))
fileList.append(os.path.join(topLevel, "TODO"))
data_files = [(FISX_DATA_DIR, fileList)]

# actual build stuff
FORCE = False
cython_dir = os.path.join(topLevel, "python", "cython")
if build_ext:
    #make sure everything is totally clean?
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
        pyx = sorted(pyx, key=str.lower)
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
    inSrc = os.path.join(cython_dir, 'default', '_fisx.cpp')
    inFile = open(inSrc, 'rb')
    inLines = inFile.readlines()
    inFile.close()
    outSrc = os.path.join(cython_dir, '_fisx.cpp')
    if os.path.exists(outSrc):
        outFile = open(outSrc, 'rb')
        outLines = outFile.readlines()
        outFile.close()
        if outLines != inLines:
            os.remove(outSrc)
    if not os.path.exists(outSrc):
        outFile = open(outSrc, 'wb')
        outFile.writelines(inLines)
        outFile.close()
    src = [outSrc]

src += glob.glob(os.path.join(topLevel,
                'src', 'fisx_*.cpp'))

include_dirs = [numpy.get_include(),
                os.path.join(topLevel, "src")]

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


class sdist_debian(sdist):
    """
    Tailor made sdist for debian
    * remove auto-generated doc
    * remove cython generated .c files
    * remove cython generated .c files
    * remove .bat files
    * include .l man files
    """
    @staticmethod
    def get_debian_name():
        name = "python-fisx_%s" % (__version__)
        return name

    def prune_file_list(self):
        sdist.prune_file_list(self)

        # this is for Cython files specifically: remove C & html files
        search_root = topLevel
        for root, _, files in os.walk(search_root):
            for afile in files:
                if os.path.splitext(afile)[1].lower() == ".pyx":
                    base_file = os.path.join(root, afile)[len(search_root) + 1:-4]
                    self.filelist.exclude_pattern(pattern=base_file + ".c")
                    self.filelist.exclude_pattern(pattern=base_file + ".cpp")
                    self.filelist.exclude_pattern(pattern=base_file + ".html")

    def make_distribution(self):
        self.prune_file_list()
        sdist.make_distribution(self)
        dest = self.archive_files[0]
        dirname, basename = os.path.split(dest)
        base, ext = os.path.splitext(basename)
        while ext in [".zip", ".tar", ".bz2", ".gz", ".Z", ".lz", ".orig"]:
            base, ext = os.path.splitext(base)
        if ext:
            dest = "".join((base, ext))
        else:
            dest = base
        # sp = dest.split("-")
        # base = sp[:-1]
        # nr = sp[-1]
        debian_arch = os.path.join(dirname, self.get_debian_name() + ".orig.tar.gz")
        os.rename(self.archive_files[0], debian_arch)
        self.archive_files = [debian_arch]
        print("Building debian .orig.tar.gz in %s" % self.archive_files[0])

cmdclass = {'install_data': smart_install_data,
            'build_py': smart_build_py,
            'debian_src': sdist_debian}
if build_ext:
    cmdclass['build_ext'] = build_ext

description = "Quantitative X-Ray Fluorescence Analysis Support Library"
with open("README.rst") as f:
    long_description = f.read()

# tell distutils where to find the packages
package_dir = {"": "python"}
packages = ['fisx', 'fisx.tests']

classifiers = ["Development Status :: 5 - Production/Stable",
               "Programming Language :: C++",
               "Programming Language :: Python :: 2",
               "Programming Language :: Python :: 3",
               "Programming Language :: Cython",
               "Environment :: Console",
               "Intended Audience :: Developers",
               "Intended Audience :: End Users/Desktop",
               "Intended Audience :: Science/Research",
               "License :: OSI Approved :: MIT License",
               "Topic :: Software Development :: Libraries :: Python Modules",
               "Operating System :: Microsoft :: Windows",
               "Operating System :: Unix",
               "Operating System :: MacOS :: MacOS X",
               "Operating System :: POSIX",
               "Topic :: Scientific/Engineering :: Chemistry",
               "Topic :: Scientific/Engineering :: Physics"
               ]

install_requires = ["numpy"]

# cython is not mandatory to build because generated code is supplied
setup_requires = ["numpy"]

if __name__ == "__main__":
    setup(
        name='fisx',
        version=__version__,
        author="V. Armando Solé",
        author_email="sole@esrf.fr",
        description=description,
        long_description=long_description,
        license="MIT",
        url="https://github.com/vasole/fisx",
        download_url="https://github.com/vasole/fisx/archive/v%s.tar.gz" % __version__, 
        package_dir=package_dir,
        packages=packages,
        ext_modules=ext_modules,
        data_files=data_files,
        cmdclass=cmdclass,
        classifiers=classifiers,
        install_requires=install_requires,
        setup_requires=setup_requires,
    )
