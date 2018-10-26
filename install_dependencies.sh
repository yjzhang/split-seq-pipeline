#!/bin/sh

# this tries to install all the necessary dependencies...
# first argument is directory where stuff is to be installed

install_dir=$1
cd $install_dir

mkdir bin

# append to path?
if echo $PATH | grep -q $1; then
    echo "PATH contains dir"
else
    export PATH=$1/bin:$PATH
    echo "\nexport PATH=$1/bin:$PATH" >> ~/.bashrc
fi

# install STAR
echo "Installing STAR..."
# TODO: test if STAR already exists

wget https://github.com/alexdobin/STAR/archive/2.6.1c.tar.gz
tar -xzf 2.6.1c.tar.gz
cd STAR-2.6.1c
cd source
make STAR
cp STAR ../../bin

cd $install_dir

# install starcode
echo "Installing starcode..."

wget https://github.com/gui11aume/starcode/archive/1.3.tar.gz
tar -xzf 1.3.tar.gz
cd starcode-1.3
cd starcode
make
cp starcode ../bin

cd $install_dir

# install samtools

wget https://github.com/samtools/samtools/releases/download/1.9/samtools-1.9.tar.bz2
tar -xf samtools-1.9.tar.bz2
cd samtools-1.9
./configure --prefix=$install_dir/bin
make
make install
cp samtools ../bin

cd $install_dir

# install picard
wget https://github.com/broadinstitute/picard/releases/download/2.18.14/picard.jar
mv picard.jar bin/


