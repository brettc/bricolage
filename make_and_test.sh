# TODO: rewrite this in python for portability
# or maybe figure out how to do it all with waf
# First, build everything
./waf build

# Change to the build folder, then copy the built stuff back into the
# development folder
cd build
find . -name \*.so -o -name \*.dylib | rsync -q -a -vv --files-from=- . ../
cd ..

# Now run some tests
py.test -v tests
