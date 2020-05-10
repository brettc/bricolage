import os
import glob

from setuptools import setup, Extension
from Cython.Build import cythonize

ext_lib_path = "src"
# include_dir = os.path.join(ext_lib_path, 'include')
include_dir = ext_lib_path
sources = glob.glob(os.path.join(ext_lib_path, "*.cpp"))
print(sources)
inc_python = "/home/ubuntu/miniconda3/envs/bricolage/include"

import numpy
inc_numpy = numpy.get_include()
print(inc_numpy)
# inc_numpy = "/home/ubuntu/miniconda3/envs/bricolage/lib/python2.7/site-packages/numpy/core/include"

# setup(
#             name="My hello app",
#                 ext_modules=cythonize('hello.pyx', compiler_directives={'embedsignature': True}),
#                 )

# Use as macros = [('<DEFINITION>', '<VALUE>')]
# where value can be None
macros = None

ext_libraries = [
    [
        "grn",
        {
            "sources": sources,
            "include_dirs": [include_dir, inc_python],
            "macros": macros,
        },
    ]
]

extensions = [
    Extension(
        "test",
        sources=["bricolage/test.pyx"],
        language="c++",
        include_dirs=['./src'],
        libraries=["grn"],
    ),
    # Extension(
    #     "core_ext",
    #     sources=["bricolage/core_ext.pyx"],
    #     language="c++",
    #     include_dirs=['./src'],
    #     libraries=["grn"],
    # ),
    # Extension(
    #     "analysis_ext",
    #     sources=["bricolage/analysis_ext.pyx"],
    #     language="c++",
    #     include_dirs=['./src', inc_python, '.', inc_numpy],
    #     libraries=["grn"],
    # ),
    # Extension(
    #     "logic2_ext",
    #     sources=["bricolage/logic2_ext.pyx"],
    #     language="c++",
    #     include_dirs=['./src', inc_python, '.', inc_numpy],
    #     libraries=["grn"],
    # ),
]

setup(ext_modules=cythonize(extensions), libraries=ext_libraries)
