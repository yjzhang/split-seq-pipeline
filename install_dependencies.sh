#!/bin/sh

# this tries to install all the necessary dependencies...
# first argument is directory where stuff is to be installed

install_dir=$1
cd $install_dir

mkdir bin

# append to path?
#if echo $PATH | grep -q $1; then
#    echo "PATH contains dir"
#else
#    export PATH=$1/bin:$PATH
#    echo "\nexport PATH=$1/bin:$PATH" >> ~/.bashrc
#fi

# install STAR
echo "Installing STAR..."
# TODO: test if STAR already exists

if [ -x "$(command -v STAR)" ]; then
    wget https://github.com/alexdobin/STAR/archive/2.6.1c.tar.gz
    tar -xzf 2.6.1c.tar.gz
    cd STAR-2.6.1c
    cd source
    make STAR
    cp STAR ../../bin
fi
echo "STAR is installed"

cd $install_dir
j
# install samtools

if [ -x "$(command -v samtools)" ]; then
    wget https://github.com/samtools/samtools/releases/download/1.9/samtools-1.9.tar.bz2
    tar -xf samtools-1.9.tar.bz2
    cd samtools-1.9
    ./configure --prefix=$install_dir/bin
    make
    make install
    cp samtools ../bin
fi

echo "samtools is installed"

