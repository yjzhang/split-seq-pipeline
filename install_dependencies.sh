#!/bin/bash

# this tries to install all the necessary dependencies...
# first argument is directory where stuff is to be installed
DEFAULT_INSTALL_DIR=~/split_seq_reqs

install_dir=$DEFAULT_INSTALL_DIR
if [ -z "$1" ]
then
    echo "Default install dir: $install_dir"
else
    install_dir=$1
fi

echo "Installing dependencies to $install_dir"

# TODO: allow for user install dir as a custom arg
cd $install_dir

mkdir bin

# append to path?
#if echo $PATH | grep -q $1; then
#    echo "PATH contains dir"
#else
#    export PATH=$1/bin:$PATH
#    echo "\nexport PATH=$1/bin:$PATH" >> ~/.bashrc
#    . ~/.bashrc
#fi

# install STAR
echo "Installing STAR..."
# TODO: test if STAR already exists

if [ -x "$(command -v STAR)" ]; then
    echo "STAR is installed"
else
    wget https://github.com/alexdobin/STAR/archive/2.6.1c.tar.gz
    tar -xzf 2.6.1c.tar.gz
    cd STAR-2.6.1c
    cd source
    make STAR
    cp STAR ../../bin
    cd $install_dir
    rm 2.6.1c.tar.gz
    echo "STAR is installed"
fi

cd $install_dir

# install samtools

if [ -x "$(command -v samtools)" ]; then
    echo "samtools is installed"
else
    wget https://github.com/samtools/samtools/releases/download/1.9/samtools-1.9.tar.bz2
    tar -xf samtools-1.9.tar.bz2
    cd samtools-1.9
    ./configure --prefix=$install_dir/
    make
    make install
    cp samtools ../bin
    cd $install_dir
    rm samtools-1.9.tar.bz2
    echo "samtools is installed"
fi

