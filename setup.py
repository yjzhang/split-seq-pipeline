from setuptools import setup, find_packages
from setuptools.command.install import install
import os

HOME = os.path.expanduser('~')

DEFAULT_INSTALL_DIR = os.path.join(HOME, 'split_seq')



class CustomInstallCommand(install):
    """Customized setuptools install command - prints a friendly greeting."""

    def run(self):
        print('test')
        install.run(self)

setup(
    name='split_seq',
    version='0.0.1',
    description='split seq tools',
    # basic stuff here
    scripts = [
            'split-seq'
    ],
    packages = ['split_seq'],
    install_requires=[
        'numpy',
        'pysam',
    ],
)
