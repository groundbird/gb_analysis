import collections
import inspect
from copy import deepcopy
import inspect

import scipy.stats
import numpy as np

import lmfit

from .expr_with_args import Expr_with_args, EWA_to_func, EWA_to_gradfunc, EWA_from_string
from .fitresult import FitResult

try:
    import sympy
except ImportError:
    sympy_installed = False
else:
    sympy_installed = True

def fit(function, xs, ys, err=None, via=None, range=None, silent=True, convert=None, **kws):
    return _fit(function, xs, ys, err, via, range, silent, convert, **kws)

def _fit(function, xs, ys, err=None, via=None, range=None, silent=True, convert=None, **kws):
    """
    hoge
    """

    ## get function
    if type(function) == Expr_with_args:
        func = EWA_to_func(function)
    elif callable(function):
        func = function
    elif isinstance(function, str):
        function = EWA_from_string(function)
        func = EWA_to_func(function)
    else:
        raise ValueError("can't get function of %s" % (function,))

    # get its argument
    funcsignature = inspect.signature(func)

    ## prepare params argument
    params = lmfit.Parameters()

    for k in list(funcsignature.parameters.values())[1:]: # skip 'x'
        if not silent : print(f'ADD PARAM -> {k.name:20s} : {k.default}')
        if k.default is not k.empty:
            params.add(k.name, value=k.default)
        else:
            params.add(k.name, value=1)

    if via is None:
        pass
    elif isinstance(via, lmfit.Parameters):
        if not silent:
            print("INSERT MANUAL PARAMETERS")
            print(via)
        params = via
    elif isinstance(via, dict):
        for k, v in list(via.items()):
            if not silent: print(f'ADD PARAM DEFAULT -> {k:20s} : {v}')
            if isinstance(v, (int, float, complex)):
                params[k].value=v
            elif isinstance(v, lmfit.Parameter):
                params[k] = v
            else:
                raise RuntimeError(f'{type(v)} for {k}: Not Implemented parameter value type')
    else:
        raise RuntimeError(f'{type(via)}: Not Implemented via type')

    ## prepare Dfun if function is an Expr_with_args object
    if type(function) == Expr_with_args:
        usenames = []
        for k in list(funcsignature.parameters.keys())[1:]:
            if params[k].vary == True:
                usenames.append(k)
        kws['Dfun'] = EWA_to_gradfunc(function, usenames)

    # if conversion function is specified, wrap func and Dfun
    if convert is not None:
        if callable(convert):
            convert_name = convert.__name__
            funcstr = f'''def wrapped({' ,'.join(list(funcsignature.parameters.keys()))}):
                ret = convert(func({', '.join(list(funcsignature.parameters.keys()))}))
                return ret'''
            Dfunstr = f'''def wrapped_Dfun({' ,'.join(list(funcsignature.parameters.keys()))}):
                ret = convert(Dfun({', '.join(list(funcsignature.parameters.keys()))}))
                return ret'''
            # make namespace for exec
            ns = globals()
            ns.update(locals())
            exec(funcstr, ns)
            if 'Dfun' in kws:
                # update kws['Dfun'] to wrapped one
                ns['Dfun'] = kws['Dfun']
                exec(Dfunstr, ns)
                kws['Dfun'] = wrapped_Dfun
            func = wrapped
        else:
            raise RuntimeError("can't call convert function %s" % convert)

    if 'Dfun' in kws:
        # https://github.com/lmfit/lmfit-py/blob/master/examples/example_derivfunc.py
        Dfun_orig = kws['Dfun']
        if err is None:
            def Dfun_pars(p, x, data=None):
                v = dict((k, v) for (k, v) in list(p.valuesdict().items()) if k in list(funcsignature.parameters.keys())[1:])
                ret = Dfun_orig(x, **v)
                return ret
        else:
            def Dfun_pars(p, x, data=None):
                v = dict((k, v) for (k, v) in list(p.valuesdict().items()) if k in list(funcsignature.parameters.keys())[1:])
                ret = Dfun_orig(x, **v)/err
                return ret
        kws['Dfun']      = Dfun_pars
        kws['col_deriv'] = 1

    # assert funcargs[0][0] == 'x'   # for now

    if err is None:
        # treat all error as 1
        def residue(p, xs):
            v = dict((k, v) for (k, v) in list(p.valuesdict().items()) if k in list(funcsignature.parameters.keys())[1:])
            return (func(xs,**v) - ys)
    else:
        # w/ err
        def residue(p, xs):
            v = dict((k, v) for (k, v) in list(p.valuesdict().items()) if k in list(funcsignature.parameters.keys())[1:])
            return (func(xs,**v) - ys)/err

    try:
        minimized = lmfit.minimize(residue, params, args=(xs,), **kws)
    except NameError as e: ## for EWA with un derivable expression
        if e.args != ("global name 'Derivative' is not defined",):
            raise
        else:
            ## retry without Dfun
            del kws["Dfun"]
            minimized = lmfit.minimize(residue, params, args=(xs,), **kws)

    minimized = lmfit.minimize(residue, params, args=(xs,), **kws)

    result    = FitResult(minimized, err, function)

    if not silent:
        result.report()

    return result
