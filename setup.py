'''
Setup for Aurora package. Basic call (install in editable mode)
pip install -e .

To install with a different Fortran compiler, use e.g.
python3 setup.py build --fcompiler=gnu95
or
python3 setup.py build --fcompiler=intelem 

It should be possible to pass any f2py flags via the command line, e.g. using
python3 setup.py build --fcompiler=intelem --opt="-fast"

'''

import setuptools
import os, sys, subprocess
from numpy.distutils.core import setup, Extension

package_name='aurorafusion'

with open('README.md', 'r') as fh:
    long_description = fh.read()

wrapper = Extension(name='_aurora', 
                    sources=['aurora/main.f90',
                             'aurora/grids.f90',
                             'aurora/impden.f90',
                             'aurora/math.f90'])

# use local makefile and avoid numpy's Extension class...
#cmd = 'make clean; make aurora'
#result = subprocess.call(cmd, shell=True) 


# load version number from .version file
aurora_dir = os.path.dirname(os.path.abspath(__file__))
with open(aurora_dir+'/.version') as vfile:
    version = vfile.read().strip()


setup(name=package_name,
      version=version,
      description=long_description,
      long_description_content_type='text/markdown',
      url='https://github.com/fsciortino/Aurora',
      author='F. Sciortino',
      author_email='sciortino@psfc.mit.edu',
      packages=['aurora'], #setuptools.find_packages(),
      requires=['numpy','scipy','matplotlib','xarray',
                'omfit_commonclasses','omfit_eqdsk','omfit_gapy'],
      ext_modules=[wrapper],
      classifiers=['Programming Language :: Python :: 3',
                   'Operating System :: OS Independent',
                   ],
      )

# move shared-object library to ./aurora
#filename = [filename for filename in os.listdir('.') if filename.startswith('_aurora')]
#print(filename)
#os.rename(filename, './aurora/'+filename)
