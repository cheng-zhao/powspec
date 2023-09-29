from distutils.extension import Extension
from distutils.core import setup
from Cython.Build import cythonize
import numpy
import glob, os

from numpy.compat import py3k
try:
      os.remove("src/pypowspec.c")
except: pass
includes = [numpy.get_include(), '/usr/include', 'etc', 'io', 'lib', 'math', 'src']
sources = glob.glob(f"io/*.c") + glob.glob(f"lib/*.c") + glob.glob(f"math/*.c") + glob.glob(f"src/*.c")
pypowspec = Extension("pypowspec",
                  sources=['src/pypowspec.pyx'] + sources,
                  include_dirs=[f"{os.environ.get('FFTW_DIR')}/../include", f"{os.environ.get('FFTW_DIR')}/include"] + includes,
                  library_dirs=[f"{os.environ.get('FFTW_DIR')}", f"{os.environ.get('FFTW_DIR')}/lib"],
                  language='c',
                  extra_compile_args=["-DOMP", "-fopenmp", "-lfftw3_omp"],
                  extra_link_args=["-fopenmp", "-lfftw3_omp"]
             )


setup(name='pypowspec',
      ext_modules=cythonize([pypowspec], gdb_debug=True),
      packages=['pypowspec'])