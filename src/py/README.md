search for 6x6 magic squares with magic product `2^10 * 3^4 * 5^3 * 7^2`
```
pypy3 rowsearch.py 10 4 3 2
```

search for 6x6 magic squares with magic product `2^10 * 3^4 * 5^3 * 7^2` and magic sum `0` or `327`
```
pypy3 rowsearch.py 10 4 3 2 --sum 0 327
```

search for 6x6 magic squares with magic product `2^10 * 3^4 * 5^3 * 7^2` in a deterministic order, without parallelism
```
pypy3 rowsearch.py 10 4 3 2 --deterministic --serial
```

dump enumerations to a file
```
pypy3 enumeration.py 10 4 3 2 --file enumeration.txt
```

post-process arrangement.c output, normal verbosity
```
pypy3 postprocess.py out.txt -v
```

post-process arrangement.c output, extra verbosity, assert that there are 17 solutions
```
pypy3 postprocess.py out.txt  -vv --expected-nsols 17
```
