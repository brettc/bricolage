"""
Cribbed and changed from this article and github site:

https://github.com/hynek/attrs/blob/master/setup.py
https://hynek.me/articles/sharing-your-labor-of-love-pypi-quick-and-dirty/#fn4
"""
import codecs
import os
import re

from setuptools import setup, find_packages


###############################################################################


NAME = "bricolage"
PACKAGES = find_packages(where=".")
META_PATH = os.path.join("bricolage", "__init__.py")
KEYWORDS = ["genetic", "evolution", "regulation"]
CLASSIFIERS = [
    "Development Status :: 4 - Alpha",
    "Intended Audience :: Science/Research",
    "Natural Language :: English",
    "License :: OSI Approved :: GNU General Public License v3 (GPLv3)",
    "Operating System :: OS Independent",
    "Programming Language :: C++",
    "Programming Language :: Python",
    "Programming Language :: Python :: 2.7",
    "Topic :: Scientific/Engineering",
]


INSTALL_REQUIRES = [
    "pathlib",
    "pandas",
    "networkx", 'numpy'
]

###############################################################################

HERE = os.path.abspath(os.path.dirname(__file__))


def read(*parts):
    """
    Build an absolute path from *parts* and and return the contents of the
    resulting file.  Assume UTF-8 encoding.
    """
    with codecs.open(os.path.join(HERE, *parts), "rb", "utf-8") as f:
        return f.read()


META_FILE = read(META_PATH)


def find_meta(meta):
    """
    Extract __*meta*__ from META_FILE.
    """
    meta_match = re.search(
        r"^__{meta}__ = ['\"]([^'\"]*)['\"]".format(meta=meta),
        META_FILE, re.M
    )
    if meta_match:
        return meta_match.group(1)
    raise RuntimeError("Unable to find __{meta}__ string.".format(meta=meta))


if __name__ == "__main__":
    setup(
        name=NAME,
        description=find_meta("description"),
        license=find_meta("license"),
        url=find_meta("uri"),
        version=find_meta("version"),
        author=find_meta("author"),
        author_email=find_meta("email"),
        maintainer=find_meta("author"),
        maintainer_email=find_meta("email"),
        keywords=KEYWORDS,
        long_description=read("README.rst"),
        packages=PACKAGES,
        zip_safe=False,
        classifiers=CLASSIFIERS,
        install_requires=INSTALL_REQUIRES,
    )
