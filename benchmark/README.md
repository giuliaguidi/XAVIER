* Clone repository recursively
```
git submodule update --init --recursive
```
* OR Update repository recursively
```
git submodule update --recursive --remote
```
* Build [libgaba](https://github.com/ocxtal/libgaba)
```
cd xavier/benchmark/tools/libgaba
make
```
* Build [parasail](https://github.com/jeffdaily/parasail#compiling-and-installing)
```
cd xavier/benchmark/tools/parasail
autoreconf -fi
./configure
make 
make install
export LD_LIBRARY_PATH=/where/you/installed/lib
```
