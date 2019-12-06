## State-of-the-Art Benchmark

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
* Build [SSW](https://github.com/mengyao/Complete-Striped-Smith-Waterman-Library)
```
cd tools/Complete-Striped-Smith-Waterman-Library/src/
make
```

## Python Semi-Global X-Drop Alignment Benchmark

Command line input:
```
-i = input sequence file 1 (">" new sequence)
-j = input sequence file 2 (">" new sequence)
-o = output file name
-x = x-drop
-m = match
-d = mismatch
-g = gap
```
Run:
```
python3 align.py -i test1.seq -j test2.seq -x 50 -o test.out -m 1 -d -1 -g -1
```

