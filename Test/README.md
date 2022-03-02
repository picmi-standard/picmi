# Testing
PICMI is currently not thoroughly tested on its own.
Feel free to contribute!

[python unittest docs](https://docs.python.org/3.8/library/unittest.html)

## Path Workaround
To allow `import picmistandard` to refer to the actual source of this repository,
this directory contains a symbolic link named `picmistandard` to the actual source.

By supplying an appropriate `PYTHONPATH` this module is loaded.

## Unittests
Unittests are launched from the `__main__` function from the unittest directory.
This tests the currently available module `picmistandard`.

The file structure follows the source 1-to-1.

To test the development version run:

```
PYTHONPATH=.:$PYTHONPATH python -m unit
```

## E2E
Execute the example as end-to-end test by launching `./launch_e2e_test.sh` from this directory.
Note that it requires the python module `fbpic` to be available.
