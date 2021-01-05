### PDS-Project-2

Run any of the three versions accordingly:
```
printf "\033c"&./run_v0.sh " -rand -n1000 -m1 -k100 -d10 -v0 "
```
```
printf "\033c"&./run_v1.sh " -rand -n1000 -m1 -k100 -d10 -v1 "
```
```
printf "\033c"&./run_v2.sh " -rand -n1000 -m1 -k100 -d10 -v2 "
```

*(m is not used on V1 and V2. A separate Query set can be created for V0 by setting `\#define V0_USE_X_AS_Y false` in `main.c` [true by default] ).*



Use a local .mtx file:

```
printf "\033c"&./run_v0.sh " -x <Path-To-Corpus-File> -y <Path-To-Query-File> -v0 "
```

If the Query set has more dimensions than the Corpus set, they will be removed import.

Optional Flags for File Imports (Use with caution! Bugs might change the results if configured wrong):

```
-transx	Transpose Corpus Dataset before use (any abnormalities in matrix sizes will be corrected).
-transy Transpose Query Dataset before use (-//-).
-n		Set the maximum number of Corpus Points (any additional will be removed from the set).
-m		Set the maximum number of Query Points (-//-).
-d		Set the number of dimensions (-//-).
-k		Set number of nearest neighbours.
```

Any of the flags `-v0`, `-v1`, `-v2` can be used with file imports. V1 and V2 use only the Corpus Dataset, but a query dataset is required as input.



Run Test with randomized data or using an .mtx file:

```
printf "\033c"&make all&./tester.o "-rand -n1000 -m1 -k100 -d14 "&make clean>/dev/null
```

*(In order for this test to work with randomized data, `#define RAND_SEED` should be false in `auxlib.h` [false by default]).*

