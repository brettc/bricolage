from distutils.core import setup
import distutils.ccompiler as cc
from distutils.cmd import Command
from distutils.extension import Extension
from Cython.Build import cythonize
import numpy
import os, os.path
import platform
import ctypes
import sys

# HACK for now
sys.argv = ['setup.py', 'build_ext', '--inplace']

class build_pubsub(Command):
    description = """build chipmunk to a shared library"""
    
    user_options = [('compiler=', 'c', 'specify the compiler type. It must understand GCC arguments')
                    ,('release', 'r', 'build chipmunk without debug asserts')
                    ]
    
    boolean_options = ['release']
    
    help_options = [
        ('help-compiler', None, "list available compilers", cc.show_compilers)
        ]

    compiler = None  
        
    def initialize_options (self):
        self.compiler = None
        self.release = True
        
    def finalize_options (self):
        pass
    
    def compile_chipmunk(self):
        if self.release:
            print("compiling chipmunk in Release mode (No debug output or asserts)" )
        else:
            print("compiling chipmunk in Debug mode (Defualt, prints debug output and asserts)")
        
        compiler = cc.new_compiler(compiler=self.compiler)

        # source_folders = ['organismal', os.path.join('chipmunk_src','constraints')]
        sources = ['organismal/pubsub2_c.cpp']
        # for folder in source_folders:
        #     for fn in os.listdir(folder):
        #         fn_path = os.path.join(folder, fn)
        #         if fn_path[-1] == 'c':
        #             sources.append(fn_path)
        #         elif fn_path[-1] == 'o':
        #             os.remove(fn_path)
                    
        include_folders = [
            './organismal',
            '/Users/Brett/anaconda/include',
            '/Users/Brett/anaconda/lib/python2.7/site-packages/numpy/core/include',
            '/Users/Brett/anaconda/include/python2.7'
        ]
        
        compiler_preargs = ['-Wno-unknown-pragmas', '-fPIC', '-Wno-unused-function', '-stdlib=libc++', '-std=c++11', '-mmacosx-version-min=10.8']

        # cc -bundle -undefined dynamic_lookup -L/Users/Brett/anaconda/lib -arch x86_64
        # -arch x86_64 build/temp.macosx-10.5-x86_64-2.7/organismal/pubsub2_ext.o
        # build/temp.macosx-10.5-x86_64-2.7/organismal/pubsub2_c.o
        # -L/Users/Brett/anaconda/lib -o
        # /Users/Brett/Dropbox/Code/organismal/organismal/pubsub2_ext.so
        
        if self.release:
            compiler_preargs.append('-DNDEBUG')
        
        # check if we are on a 64bit python
        arch = ctypes.sizeof(ctypes.c_voidp) * 8
        
        if arch == 64 and platform.system() == 'Linux':
            compiler_preargs += ['-m64', '-O3']
        elif arch == 32 and platform.system() == 'Linux':
            compiler_preargs += ['-m32', '-O3']
        elif platform.system() == 'Darwin':
            #No -O3 on OSX. There's a bug in the clang compiler when using O3.
            # compiler_preargs += ['-arch', 'i386', '-arch', 'x86_64']
            compiler_preargs += ['-O3', '-arch', 'x86_64']
        
        if platform.system() in ('Windows', 'Microsoft'):
            # Compile with stddecl instead of cdecl (-mrtd). 
            # Using cdecl cause a missing bytes issue in some cases
            # Because -mrtd and -fomit-frame-pointer (which is included in -O)
            # gives problem with function pointer to the sdtlib free function
            # we also have to use -fno-omit-frame-pointer
            compiler_preargs += ['-mrtd', '-O3', '-shared', '-fno-omit-frame-pointer'] 
        if arch == 64 and platform.system() in ('Windows', 'Microsoft'):
            compiler_preargs += ['-m64']
        if arch == 32 and platform.system() in ('Windows', 'Microsoft'):
            compiler_preargs += ['-m32']
            
        for x in compiler.executables:
            args = getattr(compiler, x)
            try:
                args.remove('-mno-cygwin') #Not available on newer versions of gcc 
                args.remove('-mdll')
            except:
                pass
        
        objs = compiler.compile(sources, include_dirs=include_folders, extra_preargs=compiler_preargs)
        
        libname = 'organismal'
        if arch == 64 and platform.system() != 'Darwin':
            libname += '64'
        if platform.system() == 'Darwin':
            libname = compiler.library_filename(libname, lib_type='dylib')
            # compiler.set_executable('linker_so', ['cc', '-dynamiclib', '-arch', 'i386', '-arch', 'x86_64'])
            # compiler.set_executable('linker_so', ['cc', '-dynamiclib', '-arch', 'x86_64', '-mmacosx-version-min=10.8', '-stdlib=libc++', '-L/Users/brett/anaconda/lib', '-lpython2.7', '-lstdc++'])
            compiler.set_executable('linker_so', ['cc', '-dynamiclib', '-arch', 'x86_64', '-mmacosx-version-min=10.8', '-std=c++11', '-stdlib=libc++', '-L/Users/brett/anaconda/lib', '-lpython2.7', '-lstdc++'])
        else:
            libname = compiler.library_filename(libname, lib_type='shared')
        linker_preargs = []
        if platform.system() == 'Linux' and platform.machine() == 'x86_64':
            linker_preargs += ['-fPIC']
        if platform.system() in ('Windows', 'Microsoft'):
            # link with stddecl instead of cdecl
            linker_preargs += ['-mrtd'] 
            # remove link against msvcr*. this is a bit ugly maybe.. :)
            compiler.dll_libraries = [lib for lib in compiler.dll_libraries if not lib.startswith("msvcr")]
        compiler.link(cc.CCompiler.SHARED_LIBRARY, objs, libname, output_dir = './organismal', extra_preargs=linker_preargs)
        
    def run(self):
        self.compile_chipmunk()

# export LDFLAGS=-lc++

# CFLAGS
# -Wno-unused-function
extensions = [
    Extension(
        "organismal/pubsub2_ext", 
        # ["organismal/pubsub2_ext.pyx"],
        # ["organismal/pubsub2_ext.pyx", "organismal/pubsub2_c.cpp"],
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
        # libraries = ['organismal'],
        # library_dirs = ['./organismal'],
    )
]
setup(
    name='organismal',
    ext_modules=cythonize(
        extensions,
        aliases={ 'NUMPY_PATH': numpy.get_include() },
        nthreads=4,
    ),
    cmdclass={'xxx' : build_pubsub}
)
