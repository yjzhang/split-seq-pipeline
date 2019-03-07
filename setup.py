from setuptools import setup, find_packages
from distutils.command.install import install
import os

HOME = os.path.expanduser('~')

# this is where the non-python dependencies will be installed.
DEFAULT_INSTALL_DIR = os.path.join(HOME, 'split_seq_reqs')



class CustomInstall(install):
    """Customized setuptools install command - prints a friendly greeting."""

    def run(self):
        print(DEFAULT_INSTALL_DIR)
        try:
            os.makedirs(DEFAULT_INSTALL_DIR)
        except:
            pass
        import subprocess
        subprocess.call('sh install_dependencies.sh', shell=True)
        super().run()


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
        'pandas',
        'scipy',
        'matplotlib',
    ],
    zip_safe=False,
    package_data={'split_seq': ['barcodes/*.csv']},
    include_package_data=True,
    # TODO: commenting this out; the user can run install_dependencies on their own.
    #cmdclass={'install': CustomInstall},
)
