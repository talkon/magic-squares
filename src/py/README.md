search for 6x6 magic squares with magic product `2^10 * 3^4 * 5^3 * 7^2`
```
pypy3 rowsearch.py 10 4 3 2
```

search for 6x6 magic squares with magic product `2^10 * 3^4 * 5^3 * 7^2` and magic sum `327`
```
pypy3 rowsearch.py 10 4 3 2 --sum 327
```

search for 6x6 magic squares with magic product `2^10 * 3^4 * 5^3 * 7^2` in a deterministic order
```
pypy3 rowsearch.py 10 4 3 2 --deterministic
```
