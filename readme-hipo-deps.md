# Dependencies install ubuntu 

## Default: static

This default setup works for the GitHub runner `ubuntu-latest`. Please replace `${{runner.workspace}}` with the path to your working directory.

1. Clone GKLib
```
git clone https://github.com/KarypisLab/GKlib.git
```

2. Clone METIS
```
git clone https://github.com/KarypisLab/METIS.git
```

3. Create installs dir 
```
mkdir installs
```

4. Install GKlib
```
cd GKlib
make config prefix=${{runner.workspace}}/installs
make
make install
```

5. Install METIS
```
cd METIS
make config prefix=${{runner.workspace}}/installs
make
make install
```

6. Check METIS and GKlib
```
cd installs
ls 
ls lib
```

7. Install BLAS or OpenBLAS
```
sudo apt update
sudo apt install libblas-dev 
```
or 
```
sudo apt update
sudo apt install libopenblas-dev 
```

8. Configure HiGHS
```
cmake -S. -B build -DHIPO=ON -DMETIS_ROOT=${{runner.workspace}}/installs -DGKLIB_ROOT=${{runner.workspace}}/installs
```

--------------

## Alternative: shared

On some systems, the following error was encountered: 

```
relocation against symbol `gk_cur_jbufs` can not be used when making a shared object, please recompile with -fPIC
``` 
Possibly, METIS is trying to compile a shared lib and to link to a static version of GKlib. It can be resolved by making all libraries shared. Delete everything in `installs` and the `HiGHS/build` dir.


4. Install GKlib shared
```
cd GKlib
make config shared=1 prefix=${{runner.workspace}}/installs
make
make install
```

Check if the shared library with no numbers in the file name is in `${{runner.workspace}}/installs/lib/`. If not, make a link to it with 
```
ln ${{runner.workspace}}/installs/lib/libGKlib.so.0 ${{runner.workspace}}/installs/lib/libGKlib.so
```

5. Install METIS shared
```
cd METIS
make config shared=1 gklib_path=${{runner.workspace}}/installs prefix=${{runner.workspace}}/installs
make
make install
```

Go to step 8.

--------

### Notes

Make sure, that in the HiGHS CMake output, the path where the GKlib library is found matches the path where the GKlib include header is. If not, there may be a system GKlib library which is being picked up by mistake. I is working on resolving this. 