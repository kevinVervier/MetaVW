# Install GDL (in gdl-1.2/GDL) #
#------------------------------#
cd ext
tar -zxf gdl-1.2.tar.gz
cd gdl-1.2
sh autogen.sh
sh autogen.sh
./configure --prefix=$PWD/GDL/
make
make install
cd ../..

# Install tools #
#---------------#
make
