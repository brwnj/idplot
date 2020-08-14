#!/bin/bash

basedir=$(pwd)

git clone https://gitlab.com/lamhm/3seq.git
cd 3seq && make -j 2
cp ./3seq /usr/local/bin

cd $basedir
rm -rf 3seq
