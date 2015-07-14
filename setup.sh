#!/bin/bash
mkdir savefiles;
mkdir wqaa;
mkdir dncuts_eigensolver;
cd wqaa;
git clone https://github.com/gvacaliuc/wqaa.git;
cd ..;
cd dncuts_eigensolver;
git clone https://github.com/gvacaliuc/dncuts-eigensolver.git;
cd ..;
echo "Directory has been intialized.  You may now run 'main.py'";
