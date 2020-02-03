#!/bin/bash
##
## Software management till upgrading to conda
##
## 03 Feb 2020
## Izaskun Mallona
## GPLv3


mkdir -p ~/soft/python
cd $_
wget https://www.python.org/ftp/python/3.6.9/Python-3.6.9.tgz
tar xzvf Python-3.6.9.tgz
cd Python-3.6.9
./configure prefix=~/soft/python/Python-3.6.9
make prefix=~/soft/python/Python-3.6.9
make install prefix=~/soft/python/Python-3.6.9
## to ensure pip with that python
# ~/soft/python/Python-3.6.9/bin/python -m ensurepip --default-pip
~/soft/python/Python-3.6.9/bin/pip3 install virtualenv


mkdir -p ~/virtenvs
cd $_
~/soft/python/Python-3.6.9/bin/virtualenv snakemake
cd $_
source bin/activate
python --version
wget https://github.com/snakemake/snakemake/archive/v5.10.0.tar.gz
tar xzvf v5.10.0.tar.gz
cd snakemake-5.10.0
python setup.py install
deactivate


mkdir -p ~/soft/bowtie
cd $_
wget https://sourceforge.net/projects/bowtie-bio/files/bowtie/1.2.3/bowtie-src-x86_64.zip
unzip bowtie-src-x86_64.zip
cd bowtie-1.2.3
sudo apt-get install libtbb-dev
make ## --prefix='/home/imallona/soft/bowtie/bowtie-1.2.3'


mkdir -p ~/soft/star
cd $_
wget https://github.com/alexdobin/STAR/archive/2.7.3a.tar.gz
tar -xzf 2.7.3a.tar.gz
cd STAR-2.7.3a
cd source
make STAR


source ~/virtenvs/snakemake/bin/activate
pip install biopython.convert
deactivate
