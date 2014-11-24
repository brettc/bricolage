from distutils.ccompiler import new_compiler

# Create compiler with default options
c = new_compiler()
workdir = "./organismal"

# Optionally add include directories etc.
# c.add_include_dir("/usr/include/c++/4.2.1")
c.add_include_dir("/usr/include/c++/4.2.1/tr1")
c.add_include_dir("/usr/local/include")

print c
print dir(c)

# Compile into .o files
objects = c.compile(["./organismal/pubsub2_c.cpp"],
                    cc_args='-std=c++11'
                    )

# Create static or shared library
# c.create_static_lib(objects, "foo", output_dir=workdir)
# c.link_shared_lib(objects, "foo", output_dir=workdir)
