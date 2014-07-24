# Notes For Me

## Getting Cython Numpy Stuff to Work

* PYTHON: These are from import numpy

```python
import numpy
dtype = numpy.int32
```

Data type | Description
---       | ---
bool_     | Boolean (True or False) stored as a byte
int_      | Default integer type (same as C long; normally either int64 or int32)
intc      | Identical to C int (normally int32 or int64)
intp      | Integer used for indexing (same as C ssize_t; normally either int32 or int64)
int8      | Byte (-128 to 127)
int16     | Integer (-32768 to 32767)
int32     | Integer (-2147483648 to 2147483647)
int64     | Integer (-9223372036854775808 to 9223372036854775807)
uint8     | Unsigned integer (0 to 255)
uint16    | Unsigned integer (0 to 65535)
uint32    | Unsigned integer (0 to 4294967295)
uint64    | Unsigned integer (0 to 18446744073709551615)
float_    | Shorthand for float64.
float16   | Half precision float: sign bit, 5 bits exponent, 10 bits mantissa
float32   | Single precision float: sign bit, 8 bits exponent, 23 bits mantissa
float64   | Double precision float: sign bit, 11 bits exponent, 52 bits mantissa

* C/C++: These are from numpy/npy_common.h

```c
#include "numpy/npy_common.h"
npy_byte blarg;
```

declaration           | typedef
---                   | ---
typedef signed char   | npy_byte;
typedef unsigned char | npy_ubyte;
typedef signed int    | npy_int;
typedef unsigned int  | npy_uint;
typedef signed long   | npy_long;
typedef unsigned long | npy_ulong;
