#!/bin/bash
# Set up environment ready for make2dgp
# usage:  ./install <install_dir>

if [ -z "$1" ]; then
	echo "export PATH=\$PATH:${PWD}" >> ~/.profile
	echo "export GP2DSOURCELOC='${PWD}'" >> ~/.profile
else
	mkdir $1
	if [ $? == 0 ]; then
		fullpath=`readlink -f $1`
		cp *.f95 $fullpath
		cp make2dgp makemovie.sh Makefile $fullpath
		echo "export PATH=\$PATH:$fullpath" >> ~/.profile
		echo "export GP2DSOURCELOC='$fullpath'" >> ~/.profile
	fi
fi

source ~/.profile
