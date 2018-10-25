

# install tools #
#---------------#
cd tools
# install external tool GDL (in gdl-1.2.1/GDL)
cd ext
tar -zxf gdl-1.2.1.tar.gz
cd gdl-1.2.1
sh autogen.sh
./configure --prefix=$PWD/GDL/
make
make install
cd ../..
# install dedicated tools 
make
# return 
cd ..


# add tools directory to path #
#-----------------------------#
export PATH=$PWD/tools:${PATH}

echo "***************************************************************"
echo "******* you should add the following in your .bashrc **********"
echo "export PATH=$PWD/tools:\${PATH}"
echo "***************************************************************"


# create alias for vw #
#---------------------#
cd tools/ext
ln -s vw-7.20150623 vw
cd ../..
# add to path
export PATH=$PWD/tools/ext:${PATH}


echo "***************************************************************"
echo "******* you should add the following in your .bashrc **********"
echo "export PATH=$PWD/tools/ext:\${PATH}"
echo "***************************************************************"


