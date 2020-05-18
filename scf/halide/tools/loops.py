#!/usr/bin/env python3

'''generate halide code to take advantage of the symmetry of g()'''

import sys

original_zones = [
    # GOOD
    {
        "name": "original",
        "iterators": [ "i", "j", "k", "l" ],
        "conditions": [],
        "updates": [
            ("i","j","k","l",1.0,"i","j","k","l"),
            ("i","j","k","l",-.5,"i","k","j","l"),
        ]
    }
]

zones = [
    # GOOD
    {
        "name": "ij_diagonal_kl_diagonal_pairs_diagonal",
        "iterators": [ "i" ],
        "conditions": [ ("i", "==", "j"), ("k", "==", "l"), ("i", "==", "k") ],
        "updates": [
            ("i","i","i","i",0.5,"i","i","i","i"),
        ]
    },
    {
        "name": "ij_diagonal_kl_diagonal_pairs_lower",
        "iterators": [ "i", "k" ],
        "conditions": [ ("i", "==", "j"), ("k", "==", "l"), ("ij", "<", "kl") ],
        "updates": [
            ("i","i","k","k",1.0,"i","i","k","k"),
            ("i","i","k","k",1.0,"k","k","i","i"),
            ("i","i","k","k",-.5,"i","k","i","k"),
            ("i","i","k","k",-.5,"k","i","k","i"),
        ]
    },
    {
        "name": "ij_diagonal_kl_lower",
        "iterators": [ "i", "j", "k", "l" ],
        "conditions": [ ("i", "==", "j"), ("k", "<", "l") ],
        "updates": [
            ("i","i","k","l",1.0,"i","i","k","l"),
            ("i","i","l","k",1.0,"i","i","l","k"),
            ("k","l","i","j",1.0,"k","l","i","j"),
            ("l","k","i","j",1.0,"l","k","i","j"),
            ("i","i","k","l",-.5,"i","k","i","l"),
            ("i","i","l","k",-.5,"i","l","i","k"),
            ("k","l","i","j",-.5,"k","i","l","j"),
            ("l","k","i","j",-.5,"l","i","k","j"),
        ]
    },
    {
        "name": "ij_lower_kl_lower_pairs_lower",
        "iterators": [ "i", "j", "k", "l" ],
        "conditions": [ ("i", "<", "j"), ("k", "<", "l"), ("ij", "<", "kl") ],
        "updates": [
            ("i","j","k","l",1.0,"i","j","k","l"),
            ("i","j","l","k",1.0,"i","j","l","k"),
            ("j","i","k","l",1.0,"j","i","k","l"),
            ("j","i","l","k",1.0,"j","i","l","k"),
            ("i","j","k","l",1.0,"k","l","i","j"),
            ("i","j","l","k",1.0,"k","l","j","i"),
            ("j","i","k","l",1.0,"l","k","i","j"),
            ("j","i","l","k",1.0,"l","k","j","i"),
            ("i","j","k","l",-.5,"k","i","l","j"),
            ("i","j","l","k",-.5,"k","j","l","i"),
            ("j","i","k","l",-.5,"l","i","k","j"),
            ("j","i","l","k",-.5,"l","j","k","i"),
            ("i","j","k","l",-.5,"i","k","j","l"),
            ("i","j","l","k",-.5,"i","l","j","k"),
            ("j","i","k","l",-.5,"j","k","i","l"),
            ("j","i","l","k",-.5,"j","l","i","k"),
        ]
    },
    {
        "name": "ij_lower_kl_lower_pair_diagonal",
        "iterators": [ "i", "j" ],
        "conditions": [ ("i", "<", "j"), ("ij", "==", "kl") ],
        "updates": [
            ("i","j","i","j",1.0,"i","j","i","j"),
            ("i","j","i","j",1.0,"i","j","j","i"),
            ("i","j","i","j",1.0,"j","i","i","j"),
            ("i","j","i","j",1.0,"j","i","j","i"),
            ("i","j","i","j",-.5,"i","i","j","j"),
            ("i","j","i","j",-.5,"i","j","j","i"),
            ("i","j","i","j",-.5,"j","i","i","j"),
            ("i","j","i","j",-.5,"j","j","i","i"),
        ]
    },
]

# for zone in zones:
#     updates = []
#     for gi, gj, gk, gl, coeff, ox, oy, dx, dy in zone['updates']:
#         updates.append({"coefficient": coeff, "g": (gi, gj, gk, gl), "out": (ox, oy), "dens": (dx, dy)})
#     zone['updates'] = updates
#     print(zone)
# exit(0)

def test_symmetry(N=8):
    # SIMPLE LOOPS VERSION

    def tally_score(i, j, our_score, quiet=False):
        if not quiet:
            print("g_fock[%d,%d] = old g_fock[%d,%d]"%(i, j, i, j))
        total_points = 0
        good_points = 0
        for key in sorted(our_score.keys()):
            total_points += 1
            si, sj, sk, sl, dx, dy = key
            record = our_score[key]
            good_coeff = record["good_coeff"]
            good_names = record["good_names"]
            symm_coeff = record["symm_coeff"]
            symm_names = record["symm_names"]
            prefix = "  "
            if good_coeff == symm_coeff:
                good_points += 1
            else:
                prefix = "!!"
            if not quiet:
                print("%s + g(%d,%d,%d,%d) * g_dens(%d,%d) * %4.1f : symmetric=%4.1f (%s → %s)"%(prefix, si, sj, sk, sl, dx, dy, good_coeff, symm_coeff, good_names, symm_names))
        if not quiet:
            print("    .... that's %d good points out of %d total."%(good_points,total_points))
        return good_points, total_points

    def sort_iters(*args):
        """find the lower triangle, or the 4-dimensional symmetry equivalant of that.  Re-orders indices to prefer i < j, i < k, k < l."""
        if len(args) == 4:
            i, j, k, l = args
            if i > j:
                j, i = i, j
            if k > l:
                l, k = k, l
            if (i, j) > (k, l):
                k, l, i, j = i, j, k, l
            return i, j, k, l
        elif len(args) == 2:
            i, j = args
            if i > j:
                j, i = i, j
            return i, j
        else:
            raise Exception("unhandled number of indexes")

    def history_to_log_dict(history):
        """input: array of log entry tuples. output: dictionary of coefficients by index tuple"""
        log = {}
        for name, gi, gj, gk, gl, coeff, dx, dy in history:
            si, sj, sk, sl = sort_iters(gi, gj, gk, gl)
            dx, dy = sort_iters(dx, dy)
            key = (si, sj, sk, sl, dx, dy)
            if key not in log:
                log[key] = (0.0, [])
            saved_coeff, saved_names = log[key]
            saved_coeff += coeff
            if name != "":
                saved_names.append(name)
            log[key] = saved_coeff, saved_names
        return log

    class FockSimulator:
        """this simulates a fock() call to determine how many operations should occur at each position in the output space.  It is used to check the results of symmetry decomposition, later."""
        def __init__(self, N):
            self.N = N
            self.matrix = []
            for i in range(N):
                row = []
                for j in range(N):
                    log = Log()
                    row.append(log)
                self.matrix.append(row)

        def __getitem__(self, t):
            if type(t) != tuple:
                raise TypeError()
            x, y = t
            return self.matrix[x][y]

        def print_log(self):
            N = self.N
            for i in range(N):
                for j in range(N):
                    print("g_fock[%d,%d] = "%(i,j))
                    print("    old g_fock[%d,%d]"%(i,j))
                    history = self[i,j].history
                    def sort_key(t):
                        name, gi, gj, gk, gl, coeff, dx, dy = t
                        gi, gj, gk, gl = sort_iters(gi, gj, gk, gl)
                        dx, dy = sort_iters(dx, dy)
                        rv = gi*N*N*N*N*N + gj*N*N*N*N + gk*N*N*N + gl*N*N + dx*N + dy
                        return rv
                    history.sort(key=sort_key)
                    for entry in history:
                            name, gi, gj, gk, gl, coeff, dx, dy = entry
                            gi, gj, gk, gl = sort_iters(gi, gj, gk, gl)
                            dx, dy = sort_iters(dx, dy)
                            sign = '+'
                            if coeff < 0:
                                sign = '-'
                                coeff = -coeff
                            if name != "":
                                name = "  # " + name
                            print("    %c %.1f  g(%d,%d,%d,%d)  g_dens(%d,%d)%s"%(sign, coeff, gi, gj, gk, gl, dx, dy, name))

    class Log():
        """Represents an element of the output matrix, logs anything that the algorithm writes to that element."""
        def __init__(self):
            self.history = []

        def log(self, name, coeff, gi, gj, gk, gl, dx, dy):
            entry = [name, gi, gj, gk, gl, coeff, dx, dy]
            self.history.append(entry)

    def match_conditions(N, conditions, i, j, k, l):
        """return True if the conditions are met for a given position in an iteration space"""
        real = { "i": i, "j": j, "k": k, "l": l, "ij": i*N+j, "kl": k*N+l }
        conds = {
            "==": lambda a, b: real[a] == real[b],
            "!=": lambda a, b: real[a] != real[b],
            "<" : lambda a, b: real[a] <  real[b],
            ">" : lambda a, b: real[a] >  real[b],
            "<=": lambda a, b: real[a] <= real[b],
            ">=": lambda a, b: real[a] >= real[b],
        }
        for a, cond, b in conditions:
            if cond in conds:
                cond = conds[cond]
                if cond(a, b) is False:
                    return False
            else:
                raise Exception("unknown cond '%s'"%cond)
        return True
    direct_fock = FockSimulator(N)
    for i in range(N):
        for j in range(N):
            for k in range(N):
                for l in range(N):
                    direct_fock[i,j].log('%d,%d,%d,%d forward'%(i,j,k,l), 1.0, i, j, k, l, k, l)
                    direct_fock[i,k].log('%d,%d,%d,%d reverse'%(i,j,k,l), -.5, i, j, k, l, j, l)


    # SYMMETRIC LOOPS VERSION
    symmetry_fock = FockSimulator(N)
    def apply_updates(name, updates, i, j, k, l):
        real = { "i": i, "j": j, "k": k, "l": l }
        for m in range(len(updates)):
            gi, gj, gk, gl, coeff, ox, oy, dx, dy = updates[m]
            #print("gi=%d gj=%d gk=%d gl=%d coeff=%f ox=%d oy=%d dx=%d dy=%d #%d"%(real[gi], real[gj], real[gk], real[gl], coeff, real[ox], real[oy], real[dx], real[dy], m))
            symmetry_fock[real[ox], real[oy]].log("%d,%d,%d,%d "%(i,j,k,l) + name + "#%d"%m, coeff, i, j, k, l, real[dx], real[dy])

    for i in range(N):
        for j in range(N):
            for k in range(N):
                for l in range(N):
                    for zone in zones:
                        if match_conditions(N, zone['conditions'], i, j, k, l):
                            apply_updates(zone['name'], zone['updates'], i, j, k, l)


    # Compare the two.

    scores = {}
    for i in range(N):
        for j in range(N):
            good_log = history_to_log_dict(direct_fock[i,j].history)
            symm_log = history_to_log_dict(symmetry_fock[i,j].history)
            score = {}
            for key in good_log:
                si, sj, sk, sl, dx, dy = key
                good_coeff, good_names = good_log[key]
                score[key] = { "good_coeff": good_coeff, "good_names": good_names, "symm_coeff": 0.0, "symm_names": [] }

            for key in symm_log:
                if key not in score:
                    score[key] = { "good_coeff": 0.0, "good_names": [] }
                coeff, names = symm_log[key]
                score[key]["symm_coeff"] = coeff
                score[key]["symm_names"] = names

            for key in score:

                if "good_coeff" not in score[key]:
                    score[key]["good_coeff"] = 0.0
                good_coeff = score[key]["good_coeff"]

                if "symm_coeff" not in score[key]:
                    score[key]["symm_coeff"] = 0.0

            good, total = tally_score(i, j, score, quiet=True)
            scores[i,j] = (good, total)

    print("--------")
    print("good/total by element:")
    failure_count = 0
    for i in range(N):
        for j in range(N):
            good, total = scores[i,j]
            print("  %3d/%3d"%(good, total), end="")
            if good != total:
                failure_count += 1
        print()
    if failure_count:
        print("%d failures found; symmetry optimizations do not equal the original computation."%failure_count)
        sys.exit(1)

    print("--------")


def halide_gen(zones):
    import halide as hl
    from math import pi

    def halide_apply_where_clause(r, a, cond, b):
        # This looks weird because the Halide python bindings don't have a way to overload 'or' and 'and'.
        # So they overload the bitwise operators and interpret them as logical operators.
        # But the precedence is weird, hence the extra parens.
        if type(a) == list:
            a = a[0] * a[1].extent() + a[1]
            b = b[0] * b[1].extent() + b[1]
        if cond == '==':
            r.where(a == b)
        if cond == '!=':
            r.where(a != b)
        if cond == '<':
            r.where(a < b)
        if cond == '<=':
            r.where(a <= b)
        if cond == '>':
            r.where(a > b)
        if cond == '>=':
            r.where(a >= b)
        print("resulting where clause:", r)

    # input scalars
    delo2  = hl.Param(hl.Float(64), "delo2")
    delta  = hl.Param(hl.Float(64), "delta")
    rdelta = hl.Param(hl.Float(64), "rdelta")

    # input vectors
    expnt_in = hl.ImageParam(hl.Float(64), 1, "expnt_in")
    rnorm_in = hl.ImageParam(hl.Float(64), 1, "rnorm_in")
    x_in     = hl.ImageParam(hl.Float(64), 1, "x_in")
    y_in     = hl.ImageParam(hl.Float(64), 1, "y_in")
    z_in     = hl.ImageParam(hl.Float(64), 1, "z_in")

    # input matrices
    fm_in        = hl.ImageParam(hl.Float(64), 2, "fm_in")
    g_fock_in_in = hl.ImageParam(hl.Float(64), 2, "g_fock_in")
    g_dens_in    = hl.ImageParam(hl.Float(64), 2, "g_dens_in")


    # output scalars
    rv = hl.Func("rv")

    # output matrix
    g_fock_out = hl.Func("g_fock_out")

    # our function's API
    all_inputs = [ delo2, delta, rdelta, expnt_in, rnorm_in, x_in, y_in, z_in, fm_in, g_fock_in_in, g_dens_in ]
    all_outputs = [ rv, g_fock_out ]

    # iterators
    i = hl.Var("i")
    j = hl.Var("j")
    k = hl.Var("k")
    l = hl.Var("l")

    nbfn = g_fock_in_in.height()

    # clamp all inputs, to prevent out-of-bounds errors from odd tile sizes and such

    expnt     = hl.BoundaryConditions.constant_exterior(expnt_in    , 0)
    rnorm     = hl.BoundaryConditions.constant_exterior(rnorm_in    , 0)
    x         = hl.BoundaryConditions.constant_exterior(x_in        , 0)
    y         = hl.BoundaryConditions.constant_exterior(y_in        , 0)
    z         = hl.BoundaryConditions.constant_exterior(z_in        , 0)
    fm        = hl.BoundaryConditions.constant_exterior(fm_in       , 0)
    g_fock_in = hl.BoundaryConditions.constant_exterior(g_fock_in_in, 0)
    g_dens    = hl.BoundaryConditions.constant_exterior(g_dens_in   , 0)

    # define g()
    dx = hl.Func("dx")
    dy = hl.Func("dy")
    dz = hl.Func("dz")
    r2 = hl.Func("r2")
    expnt2 = hl.Func("expnt2")
    expnt_inv = hl.Func("expnt_inv")

    dx[i,j] = x[i] - x[j]
    dy[i,j] = y[i] - y[j]
    dz[i,j] = z[i] - z[j]

    r2[i,j] = dx[i,j] * dx[i,j] + dy[i,j] * dy[i,j] + dz[i,j] * dz[i,j]

    expnt2[i,j]     = expnt[i] + expnt[j]
    expnt_inv[i,j] = hl.f64(1.0) / expnt2[i,j]

    fac2   = hl.Func("fac2")
    ex_arg = hl.Func("ex_arg")
    ex     = hl.Func("ex")
    denom  = hl.Func("denom")
    fac4d  = hl.Func("fac4d")
    fac2[i,j] = expnt[i] * expnt[j] * expnt_inv[i,j]
    ex_arg[i,j,k,l] = -fac2[i,j] * r2[i,j] - fac2[k,l] * r2[k,l]
    ex[i,j,k,l] = hl.select(ex_arg[i,j,k,l] < hl.f64(-37.0), hl.f64(0.0), hl.exp(ex_arg[i,j,k,l]))
    denom[i,j,k,l]  = expnt2[i,j] * expnt2[k,l] * hl.sqrt(expnt2[i,j] + expnt2[k,l])
    fac4d[i,j,k,l]  = expnt2[i,j] * expnt2[k,l] /        (expnt2[i,j] + expnt2[k,l])

    x2   = hl.Func("x2")
    y2   = hl.Func("y2")
    z2   = hl.Func("z2")
    rpq2 = hl.Func("rpq2")
    x2[i,j] = (x[i] * expnt[i] + x[j] * expnt[j]) * expnt_inv[i,j]
    y2[i,j] = (y[i] * expnt[i] + y[j] * expnt[j]) * expnt_inv[i,j]
    z2[i,j] = (z[i] * expnt[i] + z[j] * expnt[j]) * expnt_inv[i,j]
    rpq2[i,j,k,l] = (
          (x2[i,j] - x2[k,l]) * (x2[i,j] - x2[k,l])
        + (y2[i,j] - y2[k,l]) * (y2[i,j] - y2[k,l])
        + (z2[i,j] - z2[k,l]) * (z2[i,j] - z2[k,l]))

    f0t   = hl.Func("f0t")
    f0n   = hl.Func("f0n")
    f0x   = hl.Func("f0x")
    f0val = hl.Func("f0val")
    f0t[i,j,k,l] = fac4d[i,j,k,l] * rpq2[i,j,k,l]
    f0n[i,j,k,l] = hl.clamp(hl.i32((f0t[i,j,k,l] + delo2) * rdelta), fm_in.dim(0).min(), fm_in.dim(0).max())
    f0x[i,j,k,l] = delta * f0n[i,j,k,l] - f0t[i,j,k,l]
    f0val[i,j,k,l] = hl.select(f0t[i,j,k,l] >= hl.f64(28.0),
         hl.f64(0.88622692545276) / hl.sqrt(f0t[i,j,k,l]),
                                           fm[f0n[i,j,k,l],0]
         + f0x[i,j,k,l] *                 (fm[f0n[i,j,k,l],1]
         + f0x[i,j,k,l] * hl.f64(0.5) *   (fm[f0n[i,j,k,l],2]
         + f0x[i,j,k,l] * hl.f64(1./3.) * (fm[f0n[i,j,k,l],3]
         + f0x[i,j,k,l] * hl.f64(0.25) *   fm[f0n[i,j,k,l],4]))))

    g = hl.Func("g")
    g[i,j,k,l] = (hl.f64(2.00) * hl.f64(pow(pi, 2.50)) / denom[i,j,k,l]) * ex[i,j,k,l] * f0val[i,j,k,l] * rnorm[i] * rnorm[j] * rnorm[k] * rnorm[l]

    # define g_fock()
    g_fock_components = [ g_fock_in ]
    zone_funcs = {}
    for zone in zones:
        # each symmetry zone has its own iteration space, implemented as an RDom with a where() clause.
        name = "g_" + zone['name']
        # the iterators that the RDom will have
        iter_vars = {
            "i": "i",
            "j": "j",
            "k": "k",
            "l": "l",
        }
        # equalities between iterators reduce the iteration space
        for a, cond, b in zone['conditions']:
            if cond == '==':
                if a in [ 'ij', 'kl' ]:
                    for _a, _b in zip(a, b):
                        iter_vars[_b] = _a
                else:
                    iter_vars[b] = a
        print("zone %s"%name)
        for _ in range(1):
            for a in iter_vars:
                # handle 4 rounds of aliasing
                b = iter_vars[a]
                if b in iter_vars and iter_vars[b] != b:
                    iter_vars[a] = iter_vars[b]
        print("final iter mapping:", iter_vars)
        distinct_iters = sorted(set(iter_vars.values()))
        distinct_iters = [ x for x in distinct_iters if len(x) == 1 ]
        #print("distinct iters:", distinct_iters)

        # TODO: this is redundant with the real version below, consolidate this better
        updates = {}
        for gi, gj, gk, gl, coeff, oi, oj, di, dj in zone['updates']:
            okey = (oi, oj)
            dkey = (di, dj)
            gkey = (gi, gj, gk, gl)
            updatekey = (oi, oj, gi, gj, gk, gl)
            if updatekey not in updates:
                updates[updatekey] = {}
            if dkey not in updates[updatekey]:
                updates[updatekey][dkey] = 0.0
            updates[updatekey][dkey] += coeff

        piece_count = 0
        for updatekey in updates:
            for dkey in updates[updatekey].keys():
                piece_count += 1

        print("piece_count:", piece_count)
        rdom_iters = [(0, piece_count)] + ([(0, nbfn)] * len(distinct_iters))
        print("rdom iters:", rdom_iters)
        r = hl.RDom(rdom_iters, name+"_dom")
        # set local variables for RVars
        expanded_iters = {}
        distinct_iters = [r[i] for i in range(len(r))]
        assigned_already = {}
        ru = distinct_iters.pop(0)
        for a, b in iter_vars.items():
            if a in [ 'ij', 'kl' ]:
                continue
            if b in assigned_already:
                expanded_iters[a] = assigned_already[b]
            else:
                iterator = distinct_iters.pop(0)
                expanded_iters[a] = iterator
                assigned_already[b] = iterator
        gi = expanded_iters['i']
        gj = expanded_iters['j']
        gk = expanded_iters['k']
        gl = expanded_iters['l']

        # set up the correct iteration space through dynamic RDom.where() calls
        expanded_iters['ij'] = [ expanded_iters["i"], expanded_iters["j"] ]
        expanded_iters['kl'] = [ expanded_iters["k"], expanded_iters["l"] ]

        print("generating where clause conditions for %s"%name)
        for a, cond, b in zone['conditions']:
            a = expanded_iters[a]
            b = expanded_iters[b]
            if cond == '==':
                continue
            halide_apply_where_clause(r, a, cond, b)


        # consolidate multiple coefficients for the same (reduced) output and dens indexes
        updates = {}
        for gi, gj, gk, gl, coeff, oi, oj, di, dj in zone['updates']:
            oi = iter_vars[oi]
            oj = iter_vars[oj]
            di = iter_vars[di]
            dj = iter_vars[dj]
            gi = iter_vars[gi]
            gj = iter_vars[gj]
            gk = iter_vars[gk]
            gl = iter_vars[gl]
            okey = (oi, oj)
            dkey = (di, dj)
            gkey = (gi, gj, gk, gl)
            updatekey = (oi, oj, gi, gj, gk, gl)
            if updatekey not in updates:
                updates[updatekey] = {}
            if dkey not in updates[updatekey]:
                updates[updatekey][dkey] = 0.0
            updates[updatekey][dkey] += coeff

        oiv = []
        ojv = []
        div = []
        djv = []
        coeffs = []
        piece_count = 0
        for updatekey in sorted(updates.keys()):
            oi, oj, gi, gj, gk, gl = updatekey
            oi = expanded_iters[oi]
            oj = expanded_iters[oj]
            gi = expanded_iters[gi]
            gj = expanded_iters[gj]
            gk = expanded_iters[gk]
            gl = expanded_iters[gl]
            for dkey in updates[updatekey].keys():
                di, dj = dkey
                di = expanded_iters[di]
                dj = expanded_iters[dj]
                coeff = updates[updatekey][dkey]
                oiv.append(oi)
                ojv.append(oj)
                div.append(di)
                djv.append(dj)
                coeffs.append(coeff)
                print("%s[%s, %s] += g[%s, %s, %s, %s] * g_dens[%s, %s] * %f"%(name, oi.name(), oj.name(), gi.name(), gj.name(), gk.name(), gl.name(), di.name(), dj.name(), coeff))
                piece_count += 1

        # generate this symmetry zone

        def maybe_mux(s):
            if len(set(s)) == 1:
                return s[0]
            else:
                return hl.mux(hl.Expr(ru), s)
        gg = g[gi, gj, gk, gl]
        zone_func = hl.Func(name)
        zone_func[i,j] = hl.f64(0.0)
        expr = gg * g_dens[maybe_mux(div), maybe_mux(djv)] * maybe_mux(coeffs)
        zone_func[maybe_mux(oiv), maybe_mux(ojv)] += expr
        print("%s[%s, %s] += %s"%(name, maybe_mux(oiv), maybe_mux(ojv), expr))

        #if name == 'g_pairwise_symmetry':
        g_fock_components.append(zone_func)
        zone_funcs[name] = { "func": zone_func, "zone": zone, "updates": updates, "iters": expanded_iters, "rdom": r, "unroll": ru }

    expr = None
    for zone in g_fock_components:
        if expr is None:
            expr = zone[i,j]
        else:
            expr += zone[i,j]
    g_fock = hl.Func("g_fock")
    g_fock[i,j] = expr

    g_fock_out[i, j] = g_fock[i,j]

    rv[i] = hl.f64(0.0)
    r_rv = hl.RDom([(0, nbfn), (0, nbfn)])
    rv[0] += g_fock[r_rv] * g_dens[r_rv]
    rv[0] *= hl.f64(0.5)



    # scheduling
    # each zone should be one for-loop, and g should be calculated once for each inner loop iteration.
    # unroll the updates within a zone, parallelize and vectorize as much as possible
    for zone_name in zone_funcs:
        record = zone_funcs[zone_name]
        zone_func = record['func']
        zone = record['zone']
        updates = record['updates']
        expanded_iters = record['iters']
        r = record['rdom']
        distinct_iters = [r[i] for i in range(len(r))]
        ru = record['unroll']
        gi = expanded_iters['i']
        gj = expanded_iters['j']
        gk = expanded_iters['k']
        gl = expanded_iters['l']
        first = None
        print("----")
        print(zone['name'], "before scheduling")
        zone_func.print_loop_nest()

        # schedule the zone buffer initialization
        zone_func.compute_root().parallel(j).vectorize(i, 8)
        # schedule the actual work
        print("zone_func:", zone_func)
        print("distinct_iters:", distinct_iters)
        zone_func.compute_root().update(0).reorder(*distinct_iters).unroll(ru).atomic().vectorize(distinct_iters[1], 8)
        if len(distinct_iters) > 2:
            zone_func.update(0).parallel(distinct_iters[-1])

        print(zone['name'], "after scheduling")
        zone_func.print_loop_nest()

    # schedule the sum
    # consider parallelizing over wider strips of j
    g_fock.compute_root().reorder(i, j).parallel(j).vectorize(i, 8)
    # schedule the copy to output
    # consider parallelizing over wider strips of j
    g_fock_out.parallel(j).vectorize(i, 8)
    # schedule the rv reduction
    rv_intm = rv.update(0).reorder(r_rv.x, r_rv.y).rfactor([(r_rv.x, k), (r_rv.y, l)])
    rv_intm.parallel(l).vectorize(k)

    g_fock_out.print_loop_nest()

    # code generation
    p = hl.Pipeline(all_outputs)
    p.compile_to(
        {
            hl.Output.c_header: "twoel.h",
            hl.Output.c_source: "twoel.cpp",
            hl.Output.static_library: "twoel.a",
            hl.Output.stmt_html: "twoel.html",
        }, all_inputs, "twoel"
    )

    return all_outputs, all_inputs


if __name__ == "__main__":
    test_symmetry()
    outputs, inputs = halide_gen(zones)
    print({"outputs": outputs})
    print({"inputs": inputs})
