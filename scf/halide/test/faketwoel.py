#!/usr/bin/env python3

import sys
sys.path.append('../tools')
import numpy as np
import halide as hl
from loops import zones, halide_pipeline

def call_twoel(N=8, seed=2, tracing=False, tracing_g=False):
    inputs = [
        { "name": "delo2"  , "d": 0, "value": 0.001 },
        { "name": "delta"  , "d": 0, "value": 0.001 },
        { "name": "rdelta" , "d": 0, "value": 0.001 },
        { "name": "expnt"  , "d": 1, "value": 0.00001 },
        { "name": "rnorm"  , "d": 1 },
        { "name": "x"      , "d": 1 },
        { "name": "y"      , "d": 1 },
        { "name": "z"      , "d": 1 },
        { "name": "fm"     , "d": 2, "shape": [1002, 5] },
        { "name": "g_fock" , "d": 2 },
        { "name": "g_dens" , "d": 2 },
        { "name": "g_trace", "d": 4, "value": 0.0 },
    ]

    outputs = [
        { "name": "rv"    , "d": 1, "shape": [1] },
        { "name": "g_fock", "d": 2 },
    ]

    inputs  = { x["name"]: x for x in inputs }
    outputs = { x["name"]: x for x in outputs }

    buffers = []
    buffers_by_name = {}
    # be consistent
    np.random.seed(seed)
    for key in inputs:
        param = inputs[key]
        if param['d'] == 0:
            thing = 0.2
        else:
            shape = [N]*param['d']
            if 'shape' in param:
                shape = param['shape']
            thing = hl.Buffer(hl.Float(64), shape, name=key)
            if 'value' in param:
                if param['value'] != 0.0:
                    for pos in np.ndindex(*shape):
                        thing[pos] = param['value']
            else:
                values = np.random.rand(*shape) - 0.5
                for pos in np.ndindex(*shape):
                    thing[pos] = values[pos]
        buffers.append(thing)
        buffers_by_name[key] = thing

    # get JIT pipeline
    twoel, outfuncs, inparams = halide_pipeline(zones[4:5], tracing=tracing, tracing_g=tracing_g)
    twoel.compile_jit()
    # set parameter values
    for param in inparams:
        name = param.name()
        if name in buffers_by_name:
            thing = buffers_by_name[name]
        elif name.endswith("_in") and name[:-3] in buffers_by_name:
            name = name[:-3]
            thing = buffers_by_name[name]
        else:
            raise KeyError(name)
        param.set(thing)
    # evaluate it
    rv, g_fock_out = twoel.realize(N, N)
    rv = rv[0]
    g_fock_out = np.array(g_fock_out)
    print("rv", rv)
    print("g_fock_out", g_fock_out)

if __name__ == "__main__":
    call_twoel(8, 2, True, True)
