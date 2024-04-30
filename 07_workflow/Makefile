install_python:
	mkdir -p $HOME/soft/python
	cd $HOME/soft/python
	wget "https://www.python.org/ftp/python/3.12.3/Python-3.12.3.tgz"
	tar zxvf Python-3.12.3.tgz
	cd Python-3.12.3
	./configure --enable-optimizations --prefix=$HOME/soft/python/Python-3.12.3
	make -j 20
	make altinstall
install:
	export PATH=$HOME/soft/python/Python-3.12.3$PATH
	pipx install pipenv --python="/home/imallona/soft/python/Python-3.12.3/python"
	pipenv install -r requirements.txt --python="/home/imallona/soft/python/Python-3.12.3/python"
activate:
	pipenv shell
run:
	snakemake -p -s --cores 1 --configfile config.yaml
dry:
	snakemake -p --cores 1 -n -p -F --configfile config.yaml
dag:
	snakemake --dag --configfile config.yaml
graph:
	snakemake --dag --configfile config.yaml | dot -Tpng > workflow_dag.png
clean:
	rm -rf ./out
