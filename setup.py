from distutils.core import setup
from distutils.extension import Extension
from Cython.Build import cythonize
import numpy
import sys
# HACK for now
# ext_modules = [Extension(..., include_path=[numpy.get_include()])]

sys.argv = ['setup.py', 'build_ext', '--inplace']

# CFLAGS
# -Wno-unused-function
extensions = [
    Extension(
        "organismal/test", 
        ["organismal/test.pyx"],
        extra_compile_args = ['-Wno-unused-function', '-std=c++11'],
        depends = ['organismal/func.h'],
        # include_path = [numpy.get_include()],
        # include_dirs = ['./organismal'],
        # libraries = [...],
        # library_dirs = [...]),
    )
]
setup(
    name='organismal',
    ext_modules=cythonize(
        extensions,
        aliases={ 'NUMPY_PATH': numpy.get_include() },
    )
)
