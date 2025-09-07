pacman -Syu
pacman -S gsl
pacman -S mingw-w64-x86_64-gsl
pacman -S git
git clone git://git.code.sf.net/p/matio/matio
cd matio
git submodule update --init
export PATH="/c/TDM-GCC-64/bin:$PATH"
pacman -S libtool
pacman -S automake\
pacman -S autoconf m4
pacman -S make
./autogen.sh
./configure
make
export PATH="/c/Users/zimmy/Code/Mat_2_C/Compiler_Stuff/TDM-GCC-64/bin:$PATH"
make
./autogen.sh
./configure
make
export PATH="/c/Users/zimmy/Code/Mat_2_C/Compile_Stuff/TDM-GCC-64/bin:$PATH"
./configure
make
make check
make install
pacman -S mingw-w64-ucrt-x86_64-zlib
cd matio
./autogen
./autogen.sh
./configure
make
make check
make install
lear
cd matio
make clean
./autogen.sh
./configure
make
export PATH="/c/Users/zimmy/Code/Mat_2_C/Compile_Stuff/TDM-GCC-64/bin:$PATH"
make
make install
cd matio
./configure
export PATH="/c/Users/zimmy/Code/Mat_2_C/Compile_Stuff/TDM-GCC-64/bin:$PATH"
./configure
ldd /path/to/libmatio-13.dll
ldd C:/Users/zimmy/Code/Mat_2_C/Compile_Stuff/msys64/ucrt64/bin/libmatio-13.dll
cd ~
pacman -S mingw-w64-ucrt-x86_64-zlib
cd matio
make clean
./autogen
./autogen.sh
./configure --with-zlib=/ucrt64
make
make check
make install
