# How to fix missing glibc2.14 error

Encountered this on Comet which has old software and an ancient version of `blast`. First download glibc2.14 from this link:
`http://ftp.gnu.org/gnu/libc/glibc-2.14.tar.gz` and create a directory `/home/isegota/fakeroot`.

Untargz `tar -xvfz glibc-2.14.tar.gz`, go into the folder. Run:
```sh
mkdir build
cd build
../configure --prefix=/home/isegota/fakeroot
make -j4
# This may take few minutes to complete.
make install
```

How to expose this library to an executable in Linux:

```sh
LD_LIBRARY_PATH=/home/isegota/fakeroot/glibc-2.14/lib
export LD_LIBRARY_PATH
```

Add this to `~/.bash_profile` to have it available automatically in each session.