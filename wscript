# vim: ft=python
top = '.'
out = 'build'

import os
import numpy
import os.path
import glob

from waflib.Tools.ccroot import USELIB_VARS
USELIB_VARS['pyext'] = USELIB_VARS['cxxshlib'] = \
        set(['LIBPATH', 'STLIBPATH', 'LINKFLAGS', 'LINKDEPS'])
        # set(['LIB', 'STLIB', 'LIBPATH', 'STLIBPATH', 'LINKFLAGS', 'LINKDEPS'])

def options(ctx):
    ctx.load('compiler_cxx')
    ctx.load('python')
    ctx.load('cython')

def configure(ctx):
    ctx.load('compiler_cxx')
    ctx.load('python')
    ctx.check_python_headers()
    ctx.load('cython')

    # This appears to be a bug in waf -- it should know this from the cxx
    # setting. I don't see why it doesn't yet
    ctx.env.append_value('CXXFLAGS', ['-O3', '-Wno-unknown-pragmas',
                                      '-Wno-unused-function',
                                      '-stdlib=libc++', 
                                      '-std=c++11',
                                      '-mmacosx-version-min=10.8'])
    ctx.env.PREFIX = '.'

def build(ctx):
    # The Box2D library
    # ctx.env.ARCH = ['x86_64']
    # ctx(features = 'cxx cxxshlib',
    #     source   = ctx.path.ant_glob('box2d/**/*.cpp'),
    #     target   = 'cybox2d/box2d',
    #     includes = '.')

    # ctx.install_as('cybox2d/box2d.dylib', 'cybox2d/libbox2d.dylib')
    ctx(features = 'cxx cxxshlib pyext',
        # source   = ['organismal/pubsub2_ext.pyx'],
        source   = ['organismal/pubsub2_ext.pyx', 'organismal/pubsub2_c.cpp'],
        depends_on=['organismal/pubsub2_c.h'],
        target   = "./organismal/pubsub2_ext",
        includes = [
            '.', 
            './organismal', 
            '/Users/Brett/anaconda/lib/python2.7/site-packages/numpy/core/include',
        ]
        # includes = '/usr/loc'
        # use      = 'cybox2d/box2d'
       )
    ctx.install_files('./organismal', ['./organismal/pubsub2_ext.so'])

    # Extension classes
    # cybox_core(ctx, '_core')
    # cybox_core(ctx, '_core', ['callbacks.cpp']) 
    # cybox_core(ctx, '_bodies')
    # cybox_core(ctx, '_joints')
    # cybox_core(ctx, '_fixtures')

    # cybox_ext(ctx, 'cell/_type1')
    #
    # cybox_ext(ctx, 'ui/bodies')
    # cybox_ext(ctx, 'ui/connect')
    # ui_source = glob.glob('ui/*.pyx')
    # for c in ui_source:
        # cybox_ext(ctx, c)
        # print os.path.splitext(c)[0]
        # cybox_ext(ctx, os.path.splitext(c)[0])
    # cybox_ext(ctx, 'ui/sim1c')
    # cybox_ext(ctx, 'ui/sim2c'c

    # numpy_ext(ctx, 'controller/_shapespace')
    # numpy_ext(ctx, 'controller/_controller')
