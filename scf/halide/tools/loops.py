#!/usr/bin/env python3

'''generate halide code to takes advantage of the symmetry of g()'''

import sys

N = 8

#!/usr/bin/env python3

class Fock:
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

def match_conditions(N, conditions, i, j, k, l):
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


zones = [
#i < j, i < k, k < l: twoel_i_j_k_l_all_different
# All different  i < j , k < l , i < k , j < l    (1,2,3,4)   (i,j,k,l), (j,i,k,l), (i,j,l,k), (j,i,l,k), (k,l,i,j), (k,l,j,i), (l,k,i,j), (l,k,j,i)
    {
        "name": "no_symmetry",
        "conditions": [
            ("i", "<", "j"), ("i", "<=", "k"), ("k", "<", "l"), ("ij", "<", "kl")
        ],
        "updates": [
            ("i","j","k","l",1.0,"i","j","k","l"), ("j","i","k","l",1.0,"j","i","k","l"), ("i","j","l","k",1.0,"i","j","l","k"), ("j","i","l","k",1.0,"j","i","l","k"),
            ("k","l","i","j",1.0,"k","l","i","j"), ("k","l","j","i",1.0,"k","l","j","i"), ("l","k","i","j",1.0,"l","k","i","j"), ("l","k","j","i",1.0,"l","k","j","i"),
            ("i","j","k","l",-.5,"i","k","j","l"), ("j","i","k","l",-.5,"j","k","i","l"), ("i","j","l","k",-.5,"i","l","j","k"), ("j","i","l","k",-.5,"j","l","i","k"),
            ("k","l","i","j",-.5,"k","i","l","j"), ("k","l","j","i",-.5,"k","j","l","i"), ("l","k","i","j",-.5,"l","i","k","j"), ("l","k","j","i",-.5,"l","j","k","i")
        ]
    },
#i==j: twoel_i_eq_j
# Diagonal a     i == j, k < l , i < k , j < l    (1,1,3,4)   (i,j,k,l),            (i,j,l,k),            (k,l,i,j),            (l,k,i,j)
    {
        "name": "single_symmetry",
        "conditions": [
            ("i", "==", "j"), ("k", "<", "l")
        ],
        "updates": [
            ("i","j","k","l",1.0,"i","j","k","l"),                    ("i","j","l","k",1.0,"i","j","l","k"),
            ("k","l","i","j",1.0,"k","l","i","j"),                    ("l","k","i","j",1.0,"l","k","i","j"),
            ("i","j","k","l",-.5,"i","k","j","l"),                    ("i","j","l","k",-.5,"i","l","j","k"),
            ("k","l","i","j",-.5,"k","i","l","j"),                    ("l","k","i","j",-.5,"l","i","k","j")
        ]
    },
#i==k, j==l: twoel_ij_eq_kl
# Diagonal c     i < j ,       , i == k, j == l   (1,2,1,2)   (i,j,k,l), (j,i,k,l), (i,j,l,k), (j,i,l,k)
    {
        "name": "pairwise_symmetry",
        "conditions": [
            ("i", "<", "j"), ("ij", "==", "kl")
        ],
        "updates": [
            ("i","j","k","l",1.0,"i","j","k","l"), ("j","i","k","l",1.0,"j","i","k","l"), ("i","j","l","k",1.0,"i","j","l","k"), ("j","i","l","k",1.0,"j","i","l","k"),
            ("i","j","k","l",-.5,"i","k","j","l"), ("j","i","k","l",-.5,"j","k","i","l"), ("i","j","l","k",-.5,"i","l","j","k"), ("j","i","l","k",-.5,"j","l","i","k"),
        ]
    },
#i==j, k==l: twoel_i_eq_j_and_k_eq_l
# Diagonal ab    i == j, k == l, i < k            (1,1,2,2)   (i,j,k,l),                                  (k,l,i,j)
    {
        "name": "double_symmetry",
        "conditions": [
            ("i", "==", "j"), ("k", "==", "l"), ("i", "<", "k")
        ],
        "updates": [
            ("i","j","k","l",1.0,"i","j","k","l"),
            ("k","l","i","j",1.0,"k","l","i","j"),
            ("i","j","k","l",-.5,"i","k","j","l"),
            ("k","l","i","j",-.5,"k","i","l","j")
        ]
    },
#i==j==k==l: twoel_i_eq_j_eq_k_eq_l
# Diagonal abc   i == j,       , i == k, j == l   (1,1,1,1)   (i,j,k,l)
    {
        "name": "triple_symmetry",
        "conditions": [
            ("i", "==", "j"), ("i", "==", "k"), ("j", "==", "l")
        ],
        "updates": [
            ("i","j","k","l",1.0,"i","j","k","l"),
            ("i","j","k","l",-.5,"i","k","j","l"),
        ]
    },
]


# SIMPLE LOOPS VERSION
direct_fock = Fock(N)
for i in range(N):
    for j in range(N):
        for k in range(N):
            for l in range(N):
                direct_fock[i,j].log('%d,%d,%d,%d forward'%(i,j,k,l), 1.0, i, j, k, l, k, l)
                direct_fock[i,k].log('%d,%d,%d,%d reverse'%(i,j,k,l), -.5, i, j, k, l, j, l)


# SYMMETRIC LOOPS VERSION
symmetry_fock = Fock(N)
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

def print_comparison(i, j, our_score):
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
        print("%s + g(%d,%d,%d,%d) * g_dens(%d,%d) * %4.1f : symmetric=%4.1f (%s → %s)"%(prefix, si, sj, sk, sl, dx, dy, good_coeff, symm_coeff, good_names, symm_names))
    print("    .... that's %d good points out of %d total."%(good_points,total_points))
    return good_points, total_points

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
            symm_coeff = score[key]["symm_coeff"]

        good, total = print_comparison(i, j, score)
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
def fix_pair_conditions(a):
    if len(a) == 2:
        return '%s * nbfn + g%s'%(a[0], a[1])
    return a

for zone in reversed(zones):
    name = zone['name']
    print("        RDom g_%s_dom(0, nbfn, 0, nbfn, 0, nbfn, 0, nbfn);"%name)
    print("        gi = g_%s_dom[0]; gj = g_%s_dom[1]; gk = g_%s_dom[2]; gl = g_%s_dom[3];"%(name,name,name,name))
    where = []
    for a, cond, b in zone['conditions']:
        where.append("g%s %s g%s"%(fix_pair_conditions(a), cond, fix_pair_conditions(b)))
    where = " && ".join(where)
    print("        g_%s_dom.where(%s);"%(name,where))
    print("        Expr g_%s = g(gi, gj, gk, gl);"%name)
    updates = {}
    for gi, gj, gk, gl, coeff, oi, oj, di, dj in zone['updates']:
        okey = "g%s,g%s"%(oi, oj)
        dkey = "g%s,g%s"%(di, dj)
        if okey not in updates:
            updates[okey] = {}
        if dkey not in updates[okey]:
            updates[okey][dkey] = 0.0
        updates[okey][dkey] += coeff
    for okey in sorted(updates.keys()):
        values = []
        for dkey in sorted(updates[okey].keys()):
            coeff = updates[okey][dkey]
            if coeff == 1.0:
                coeff = ''
            else:
                coeff = "Expr(%.1f) * "%coeff
            values.append("%sg_dens(%s)"%(coeff,dkey))
        print("        g_fock(%s) += g_%s * (%s);"%(okey,name," + ".join(values)))
