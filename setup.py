from distutils.core import setup
from Cython.Build import cythonize
import numpy as np

setup(name="window_functions", ext_modules=cythonize('window_functions.pyx'),include_dirs=[np.get_include()])