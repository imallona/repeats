#!/bin/bash
##
## scTE
##
## https://www.nature.com/articles/s41467-021-21808-x#Sec2
##
## 3 May 2021

NTHREADS=20
cd ~/virtenvs

# requires python >= 3.7


sudo apt install build-essential zlib1g-dev libncurses5-dev libgdbm-dev \
     libnss3-dev libssl-dev libsqlite3-dev libreadline-dev libffi-dev curl libbz2-dev -y

mkdir -p ~/soft/python3
cd $_
wget https://www.python.org/ftp/python/3.9.2/Python-3.9.2.tgz
tar xzvf *tgz
cd Python-3.9.2/
./configure  --enable-optimizations  --prefix=/home/imallona/soft/python3/Python-3.9.2
make -j $NTHREADS

make install

~/soft/python3/Python-3.9.2/bin/python3 -m pip install --upgrade pip
~/soft/python3/Python-3.9.2/bin/pip3 install virtualenv

cd ~/virtenvs

~/soft/python3/Python-3.9.2/bin/virtualenv scTE

source ~/virtenvs/scTE/bin/activate

cd scTE

git clone https://github.com/JiekaiLab/scTE

cd scTE

python setup.py install

# add snakemake, too

which pip
pip install snakemake
