#fit: a module to fitting in a gnuplot-like manner

## usage

```python
import numpy as np
import matplotlib.pyplot as plt
import scipy

import fit

## define function
def testfunc(x, a, w, c):
    return a*np.sin(w*x)+c

## make test data
xs     = np.linspace(-5, 5, 50)
noise  = scipy.randn(len(xs))*0.1
ys     = 1.5*np.sin(1.7*xs) + 1.9 + noise

## fit
result = fit.fit(testfunc, xs, ys, 0.1,
                 via={'a': 1, 'w':1.5, 'c':2})

# string expression with `x` as independent variable can be used
# result = fit.fit('a*sin(w*x)+c', xs, ys, 0.1,
#                  via={'a': 1, 'w':1.5, 'c':2})

## plot
plt.errorbar(xs, ys, 0.1, fmt='.')
plt.plot(xs, testfunc(xs, **result.values()))
plt.grid()
plt.show()
```

to produce output

```
[[Fit Statistics]]
    # function evals   = 23
    # data points      = 50
    # variables        = 3
    chi-square         = 52.6636603007
    reduced chi-square = 1.12050341065
    p-value            = 0.264220474313
[[Variables]]
    a:   1.53839019 +/- 0.020580 (1.34%) (init= 1)
    w:   1.70295979 +/- 0.005136 (0.30%) (init= 1.5)
    c:   1.91410839 +/- 0.014969 (0.78%) (init= 2)
[[Correlations]] (unreported correlations are <  0.100)
```

## dependencies:
- numpy
- scipy
- sympy (if you use expression)
- lmfit >= 0.9.0
- ad

## fit data other than 1-D real array
the output of functino should be 1-D array of real number.
fit() can take `convert` keyword argument to wrap function output,
so if you wish to fit complex data as 2-d vector can do something like:
```python
import numpy as np

def complex_to_cartesian2darray(x):
    x = np.atleast_1d(x)
    shape = x.shape
    return np.concatenate([np.real(x), np.imag(x)], axis=len(shape)-1)

def cartesian2darray_to_complex(x):
    assert x.shape[-1] % 2 == 0
    size = x.shape[-1] / 2
    return x[...,:size] + 1j*x[...,size:]

def fit_complex(func, x, z, error, params):
    """fit a complex-valued function and data, by converting complex value to 2-d vector.
    give error as two-element tuple, whose first element is for real, second is for imag."""

    convert = complex_to_cartesian2darray

    if dataerror is None:
        errRI = None
    else:
        err_R = error[0]
        err_I = error[1]
        errRI = convert(np.broadcast_to(err_R, iqs.shape) +
                        1j*np.broadcast_to(err_I, iqs.shape))
    r = fit.fit(func, freqs, convert(iqs), errRI, params, silent=True, convert=convert)
    return r
```

## todo
- change module name?
