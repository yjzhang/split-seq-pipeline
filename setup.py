from setuptools import setup, find_packages
from distutils.command.install import install
import os

HOME = os.path.expanduser('~')

# this is where the dependent libraries will be installed.
DEFAULT_INSTALL_DIR = os.path.join(HOME, 'split_seq_reqs')



class CustomInstall(install):
    """Customized setuptools install command - prints a friendly greeting."""

    def run(self):
        import subprocess
        subprocess.call('sh install_dependencies.sh {0}'.format(DEFAULT_INSTALL_DIR), shell=True)
        print('test')
        super.run()

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
    cmdclass={'install': CustomInstall},
)
