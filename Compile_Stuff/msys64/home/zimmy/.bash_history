pacman -Syu
pacman -S mingw-w64-x86_64-gsl
pacman -S git
pacman -S mingw-w64-ucrt-x86_64-zlib
git clone git://git.code.sf.net/p/matio/matio
cd matio
git submodule update --init
export PATH="/c/Users/zimmy/Code/Matlab_to_C/Compile_Stuff/TDM-GCC-64/bin:$PATH"
pacman -S libtool
pacman -S automake
pacman -S autoconf m4
pacman -S make
./autogen.sh
./configure
make
make check
make install
