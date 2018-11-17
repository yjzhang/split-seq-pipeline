from setuptools import setup
import os

HOME = os.path.expanduser('~')

DEFAULT_INSTALL_DIR = os.path.join(HOME, 'split_seq')

setup(
    # basic stuff here
    scripts = [
            'split-seq'
    ],
)
