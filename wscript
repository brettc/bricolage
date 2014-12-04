top = '.'
out = 'build'

import os
import numpy
import os.path
import glob

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
    ctx.env.append_unique('CYTHONFLAGS', '--cplus')
    ctx.env.append_value('CXXFLAGS', ['-O3', '-Wno-unknown-pragmas',
                                      '-Wno-unused-function',
                                      '-stdlib=libc++', 
                                      '-std=c++11',
                                      '-mmacosx-version-min=10.8'])
    # ctx.env.append_value('LIBFLAGS', ['-O3', '-Wno-unknown-pragmas',
    #                                   '-Wno-unused-function',
    #                                   '-stdlib=libc++', 
    #                                   '-std=c++11',
    #                                   '-mmacosx-version-min=10.8'])

def cybox_core(ctx, nm, additional=[]):
    source = 'cybox2d/%s.pyx' % nm
    if additional:
        addit = [source] + ['cybox2d/%s' % a for a in additional]
        source = ' '.join(addit)

    ctx(features = 'cxx cxxshlib pyext',
        source   = source,
        target   = "./organismal/%s" % nm,
        includes = '. ./organismal',
        # use      = 'cybox2d/box2d'
       )

def cybox_ext(ctx, nm):
    source = '%s.pyx' % nm

    # The build happens in the same folder as the file. So we need to append
    # an include that references the parent. NOT sure why it doesn't pick this
    # up from the "includes"
    ctx.env.append_unique('CYTHONFLAGS', '-I../')
    ctx.env.append_unique('CCFLAGS', '-mmacosx-version-min=10.8')
    ctx(features = 'cxx cxxshlib pyext',
        source   = source,
        target   = nm,
        includes = '. cybox2d',
        use      = 'cybox2d/box2d'
       )
    # start_dir = bld.path.find_dir('src/bar')
    # ctx.install_files('.', ['%s.so' % nm], relative_trick=True) 

def numpy_ext(ctx, nm):
    source = '%s.pyx' % nm
    ctx(features = 'cxx cxxshlib pyext',
        source   = source,
        target   = nm,
        includes = numpy.get_include()
       )

def build(ctx):
    # The Box2D library
    # ctx.env.ARCH = ['x86_64']
    # ctx(features = 'cxx cxxshlib',
    #     source   = ctx.path.ant_glob('box2d/**/*.cpp'),
    #     target   = 'cybox2d/box2d',
    #     includes = '.')

    # ctx.install_as('cybox2d/box2d.dylib', 'cybox2d/libbox2d.dylib')
    ctx(features = 'cxx cxxshlib pyext',
        source   = ['organismal/pubsub2_ext.pyx', 'organismal/pubsub2_c.cpp'],
        target   = "./organismal/pubsub2_ext",
        includes = [
            '.', 
            './organismal', 
            '/Users/Brett/anaconda/lib/python2.7/site-packages/numpy/core/include',
        ]
        # includes = '/usr/loc'
        # use      = 'cybox2d/box2d'
       )

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
