from distutils.core import setup
from distutils.extension import Extension
from Cython.Build import cythonize
import numpy
import sys

# HACK for now
sys.argv = ['setup.py', 'build_ext', '--inplace']

# export LDFLAGS=-lc++

# CFLAGS
# -Wno-unused-function
extensions = [
    Extension(
        "organismal/pubsub2_ext", 
        ["organismal/pubsub2_ext.pyx", "organismal/pubsub2_c.cpp", "organismal/scheme_cooperative.cpp"],
        extra_compile_args = [
            # '-ffast-math',
            '-Wno-unused-function', 
            '-stdlib=libc++',
            '-std=c++11', 
            '-mmacosx-version-min=10.8',
        ],
        depends = ['organismal/pubsub2_c.h'],
        language = 'c++',
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
        nthreads=4,
    )
)
