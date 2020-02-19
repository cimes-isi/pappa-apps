/* Copyright (C) 2013 ET International, Inc. */

#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <math.h>
#include "tce.h"

#define CHECK_R1
#define CHECK_R2

#define EPSILON 1e-15

/* The extra nesting is to convince cpp to generate a function name like cc2_t1,
 * rather than MODEL_t1 or errors like "stray ‘##’ in program" and
 * "unknown type name ‘cc2’" */
#define _FUNC2(model,suffix) model ## _ ## suffix
#define _FUNC1(model,suffix) _FUNC2(model,suffix)
#define  FUNC(suffix) _FUNC1(MODEL,suffix)

/* Similarly, the extra nesting convinces cpp to generate a string like "cc2",
 * rather than "macro" or "MODEL" */
#define _STRINGIFY1(macro) #macro
#define STRINGIFY(macro) _STRINGIFY1(macro)

void FUNC(t1)(double* d_f1,double* d_i0,double* d_t1,double* d_t2,double* d_v2,int* k_f1_offset,int* k_i0_offset,int* k_t1_offset,int* k_t2_offset,int* k_v2_offset);
void FUNC(t2)(double* d_f1,double* d_i0,double* d_t1,double* d_t2,double* d_v2,int* k_f1_offset,int* k_i0_offset,int* k_t1_offset,int* k_t2_offset,int* k_v2_offset);

#define ENTRIES(a) (sizeof(a)/sizeof(a[0]))
static char *scalar_names[] = (char*[]) {
    "noa",
    "nob",
    "nva",
    "nvb",
    "noab",
    "nvab",
    "io_v2",
    "intorb",
    "size_f1",
    "size_t1",
    "size_t2",
    "size_v2",
    "restricted",
};
enum {
    si_noa = 0,
    si_nob,
    si_nva,
    si_nvb,
    si_noab,
    si_nvab,
    si_io_v2,
    si_intorb,
    si_size_f1,
    si_size_t1,
    si_size_t2,
    si_size_v2,
    si_restricted,
    NUM_SCALARS
};
static char *vector_names[] = (char*[]) {
    "k_sym",
    "k_b2am",
    "k_spin",
    "k_alpha",
    "k_range",
    "k_offset",
    "k_f1_offset",
    "k_t1_offset",
    "k_t2_offset",
    "k_v2_offset",
    "k_sym_alpha",
    "k_spin_alpha",
    "k_range_alpha",
    "k_v2_alpha_offset",
};
enum {
    vi_sym = 0,
    vi_b2am,
    vi_spin,
    vi_alpha,
    vi_range,
    vi_offset,
    vi_f1_offset,
    vi_t1_offset,
    vi_t2_offset,
    vi_v2_offset,
    vi_sym_alpha,
    vi_spin_alpha,
    vi_range_alpha,
    vi_v2_alpha_offset,
    NUM_VECTORS
};
static char *tensor_names[] = (char*[]) {
    "f1",
    "v2",
    "t1",
    "t2",
    "r1",
    "r2",
};
enum {/* Note: this determines the order of command line arguments. */
    ti_f1 = 0,
    ti_v2,
    ti_t1,
    ti_t2,
    ti_r1,
    ti_r2,
    NUM_TENSORS
};

int *k_sym;
int *k_b2am;
int *k_spin;
int *k_alpha;
int *k_range;
int *k_offset;
int *k_sym_alpha;
int *k_spin_alpha;
int *k_range_alpha;
int *k_v2_alpha_offset;
int noab, nvab, noa, nob, nva, nvb, io_v2;
bool restricted, intorb;

double FACTORIAL[20];
int *scalars, *vector_sizes, **vectors, *tensor_sizes;
double **tensors;

void tce_hash_v2_prepare(const int *hash, int M);

int main(int argc, char **argv) {
    FILE *f;
    int i, j, rv = 0;
    char line[256];
    scalars      = calloc(NUM_SCALARS, sizeof(int));
    vector_sizes = calloc(NUM_VECTORS, sizeof(int));
    vectors      = calloc(NUM_VECTORS, sizeof(int*));
    tensor_sizes = calloc(NUM_TENSORS, sizeof(int));
    tensors      = calloc(NUM_TENSORS, sizeof(double*));
    if(!scalars || !vector_sizes || !vectors || !tensor_sizes)
        die("malloc failure");
    if(argc < 8)
        die("Usage: %s <params> <f1> <v2> <t1> <t2> <r1> <r2>", argv[0]);
    f = fopen(argv[1], "r");
    if(!f)
        die("Cannot open params file %s", argv[1]);
    if(!fgets(line, sizeof(line), f))
        die("Could not read line of params file");
    /* parse scalar constants */
    uint64_t mask = 0;
    while(!strchr(line, '[')) {
        /* line[] contains something like "size_v2=   7970724", parse it. */
        for(i = 0; i < NUM_SCALARS; i++) {
            if(!strncmp(line, scalar_names[i], strlen(scalar_names[i])) && (
                    line[strlen(scalar_names[i])] == ' '
                 || line[strlen(scalar_names[i])] == '='))
                break;
        }
        if(i == NUM_SCALARS)
            die("Could not find scalar variable name in line \"%s\"", line);
        mask |= (1ULL<<i);
        char *equals = strchr(line,'=');
        if(!equals)
            die("Could not find equal sign in line \"%s\"", line);
        scalars[i] = strtol(equals+1,NULL,0);

        if(!fgets(line, sizeof(line), f))
            die("Could not read line of params file");
    }
    for(i = 0; i < NUM_SCALARS; i++) {
        if(!(mask & (1ULL<<i)))
            printf("scalar \"%s\" not found in params file, defaulting to 0.\n", scalar_names[i]);
    }
    /* parse vector constants */
    mask = 0;
    while(1) {
        /* line[] contains something like "k_sym[18]=", parse it. */
        for(i = 0; i < NUM_VECTORS; i++) {
            if(!strncmp(line, vector_names[i], strlen(vector_names[i]))
            && line[strlen(vector_names[i])] == '[')
                break;
        }
        if(i == NUM_VECTORS)
            die("Could not find vector variable name in line \"%s\"", line);
        mask |= (1ULL<<i);
        char *bracket = strchr(line,'[');
        if(!bracket)
            die("Could not find square brackets in line \"%s\"", line);
        vector_sizes[i] = strtol(bracket+1,NULL,0);
        vectors[i] = malloc(sizeof(vectors[i][0])*vector_sizes[i]);
        if(!vectors[i])
            die("Could not malloc buffer of size %zd for vector %s (%d)",
                sizeof(vectors[i][0])*vector_sizes[i], vector_names[i], i);
        for(j = 0; j < vector_sizes[i]; j++) {
            if(!fgets(line, sizeof(line), f))
                die("Could not read element %d of vector %s (%d)",
                    j, vector_names[i], i);
            /* line[] contains something like "     1     0", parse it. */
            char *ptr;
            int index, value;
            index = strtol(line, &ptr, 0);
            value = strtol(ptr , NULL, 0);
            if(j != index-1)
                die("index mismatch in %s: %d != %d", vector_names[i],j+1,index);

            vectors[i][j] = value;
        }

        if(!fgets(line, sizeof(line), f)) {
            fclose(f);
            break;
        }
    }
    for(i = 0; i < NUM_VECTORS; i++) {
        if(!(mask & (1ULL<<i)))
            printf("vector \"%s\" not found in params file, defaulting to NULL.\n", vector_names[i]);
    }
    /* read tensors */
    tensor_sizes[0] = scalars[si_size_f1]; /* size_f1 */
    tensor_sizes[1] = scalars[si_size_v2]; /* size_v2 */
    tensor_sizes[2] = scalars[si_size_t1]; /* size_t1 */
    tensor_sizes[3] = scalars[si_size_t2]; /* size_t2 */
    tensor_sizes[4] = scalars[si_size_t1]; /* r1 has the same size as t1 */
    tensor_sizes[5] = scalars[si_size_t2]; /* r2 has the same size as t2 */
    for(i = 0; i < 6; i++) {
        f = fopen(argv[2+i], "r");
        if(!f)
            die("Cannot open %s file %s", tensor_names[i], argv[2+i]);
        tensors[i] = malloc(sizeof(tensors[i][0]) * tensor_sizes[i]);
        if(!tensors[i])
            die("Could not malloc buffer of size %zd for tensor %s (%d)",
                sizeof(tensors[i][0])*tensor_sizes[i], tensor_names[i], i);
        int rv = fread(tensors[i], sizeof(tensors[i][0]), tensor_sizes[i], f);
        if(rv != tensor_sizes[i])
            die("reading tensor %s (%d), fread() returned %d, wanted %d",
                    tensor_names[i], i, rv, tensor_sizes[i]);
        fclose(f);
    }
    /* Set up global scalars */
    noa        = scalars[si_noa];
    nob        = scalars[si_nob];
    nva        = scalars[si_nva];
    nvb        = scalars[si_nvb];
    noab       = scalars[si_noab];
    nvab       = scalars[si_nvab];
    io_v2      = scalars[si_io_v2];
    intorb     = scalars[si_intorb];
    restricted = scalars[si_restricted];
    /* Set up global vectors */
    k_sym             = vectors[vi_sym];
    k_b2am            = vectors[vi_b2am];
    k_spin            = vectors[vi_spin];
    k_alpha           = vectors[vi_alpha];
    k_range           = vectors[vi_range];
    k_offset          = vectors[vi_offset];
    k_sym_alpha       = vectors[vi_sym_alpha];
    k_spin_alpha      = vectors[vi_spin_alpha];
    k_range_alpha     = vectors[vi_range_alpha];
    k_v2_alpha_offset = vectors[vi_v2_alpha_offset];
    FACTORIAL[0] = 1;
    for(i = 1; i < ENTRIES(FACTORIAL); i++)
        FACTORIAL[i] = FACTORIAL[i-1] * i;

    if(intorb)
        tce_hash_v2_prepare(k_v2_alpha_offset, 10);

#ifdef CHECK_R1
    printf("Calling " STRINGIFY(MODEL) "_t1 (size=%d)\n", scalars[si_size_t1]);
    /* Call cc2_t1 or ccsd_t1 or whatever */
    double *my_r1 = calloc(scalars[si_size_t1], sizeof(double));
    FUNC(t1)
           (tensors[ti_f1],my_r1,tensors[ti_t1],tensors[ti_t2],tensors[ti_v2],
            vectors[vi_f1_offset],vectors[vi_t1_offset],vectors[vi_t1_offset],
            vectors[vi_t2_offset],vectors[vi_v2_offset]);
#endif

#ifdef CHECK_R2
    /* Call cc2_t2 or ccsd_t2 or whatever */
    printf("Calling " STRINGIFY(MODEL) "_t2 (size=%d)\n", scalars[si_size_t2]);
    double *my_r2 = calloc(scalars[si_size_t2], sizeof(double));
    FUNC(t2)(tensors[ti_f1],my_r2,tensors[ti_t1],tensors[ti_t2],tensors[ti_v2],
            vectors[vi_f1_offset],vectors[vi_t2_offset],vectors[vi_t1_offset],
            vectors[vi_t2_offset],vectors[vi_v2_offset]);
#endif /* CHECK_R2 */

#ifdef CHECK_R1
    /* Check r1 */
    printf("Checking " STRINGIFY(MODEL) "_t1\n");
    for(i = 0; i < scalars[si_size_t1]; i++) {
        if(isnan(my_r1[i]))
            die("got NAN!");
        if(fabs(tensors[ti_r1][i] - my_r1[i]) > EPSILON) {
            rv++;
            if(rv == 20) {
                debug("[more t1 errors suppressed.]");
                continue;
            } else if(rv < 20)
                debug(STRINGIFY(MODEL) "_t1 output index %d differs: (expected, got) %e, %e",
                      i, tensors[ti_r1][i], my_r1[i]);
        }
    }
#endif /* CHECK_R1 */

    if(rv >= 20)
        rv = 19;

#ifdef CHECK_R2
    /* Check r2 */
    printf("Checking " STRINGIFY(MODEL) "_t2\n");
    for(i = 0; i < scalars[si_size_t2]; i++) {
        if(isnan(my_r2[i]))
            die("got NAN!");
        if(fabs(tensors[ti_r2][i] - my_r2[i]) > EPSILON) {
            rv++;
            if(rv == 20) {
                debug("[more t2 errors suppressed.]");
                continue;
            } else if(rv < 20)
                debug(STRINGIFY(MODEL) "_t2 output index %d differs: (expected, got) %e, %e",
                        i, tensors[ti_r2][i], my_r2[i]);
        }
    }
#endif /* CHECK_R2 */

    if(rv > 50)
        rv = 50;

    /* Free memory */
    for(i = 0; i < NUM_TENSORS; i++)
        if(tensors[i])
            free(tensors[i]);
    for(i = 0; i < NUM_VECTORS; i++)
        if(vectors[i])
            free(vectors[i]);
#ifdef CHECK_R1
    free(my_r1);
#endif
#ifdef CHECK_R2
    free(my_r2);
#endif /* CHECK_R2 */
    free(vector_sizes);
    free(tensor_sizes);
    free(scalars);
    free(vectors);
    free(tensors);

    printf("%d errors detected in total.\n", rv);
    return rv;
}

/* Straight translations of src/tce/tce_restricted.F, minus the int_mb nonsense */
void tce_restricted_2(int a1b, int a2b, int *b1b, int *b2b) {
    if(restricted && k_spin[a1b]+k_spin[a2b] == 4) {
        *b1b = k_alpha[a1b] - 1;
        *b2b = k_alpha[a2b] - 1;
    } else {
        *b1b = a1b;
        *b2b = a2b;
    }
}
void tce_restricted_4(int a1b, int a2b, int a3b, int a4b, int *b1b, int *b2b, int *b3b, int *b4b) {
    if(restricted && k_spin[a1b]+k_spin[a2b]+k_spin[a3b]+k_spin[a4b] == 8) {
        *b1b = k_alpha[a1b] - 1;
        *b2b = k_alpha[a2b] - 1;
        *b3b = k_alpha[a3b] - 1;
        *b4b = k_alpha[a4b] - 1;
    } else {
        *b1b = a1b;
        *b2b = a2b;
        *b3b = a3b;
        *b4b = a4b;
    }
}


/* Rearrange the 2d input array, multiplying values by the amplitude factor.
 * [a,b] defines the length, along each dimension, of the input arrays.
 * [i,j] defines the target index order.  [0,1] is a direct copy; [1,0] transposes. */
void tce_sort_2(double *in, double *out, int a, int b, int i, int j, double factor) {
    /* more-or-less direct translation of tce_sort2.F */
    int pos[2],len[2],ia,ib,j1,j2;
    ia = 0;
    if(i == 0 && j == 1) {
        /* no permutation, just copy and apply the scale factor */
        for(ia = 0; ia < a * b; ia++)
            out[ia] = in[ia] * factor;
        return;
    }
    len[0] = a;
    len[1] = b;
    for(j1 = 0; j1 < a; j1++) {
        pos[0] = j1;
        for(j2 = 0; j2 < b; j2++) {
            pos[1] = j2;
            ib = pos[j]+len[j]*pos[i];
            out[ib] = in[ia++] * factor;
        }
    }
}

/* Rearrange the 4d input array, multiplying values by the amplitude factor.
 * [a,b,c,d] defines the length, along each dimension, of the input arrays.
 * [i,j,k,l] defines the target index order.  [0,1,2,3] is a direct copy. */
void tce_sort_4(double *in, double *out, int a, int b, int c, int d, int i, int j, int k, int l, double factor) {
    /* a simple extrapolation of TCE_SORT_2.  The original TCE has two
     * tce_sort_4*.F sources, it's not clear which one is used, and they are
     * both much more heavily optimized than this. */
    int pos[4],len[4],ia,ib,j1,j2,j3,j4;
    ia = 0;
    len[0] = a;
    len[1] = b;
    len[2] = c;
    len[3] = d;
    for(j1 = 0; j1 < a; j1++) {
        pos[0] = j1;
        for(j2 = 0; j2 < b; j2++) {
            pos[1] = j2;
            for(j3 = 0; j3 < c; j3++) {
                pos[2] = j3;
                for(j4 = 0; j4 < d; j4++) {
                    pos[3] = j4;
                    ib = len[j]*len[k]*len[l]*pos[i]
                       +        len[k]*len[l]*pos[j]
                       +               len[l]*pos[k]
                       +                      pos[l];
                    out[ib] = in[ia++] * factor;
                }
            }
        }
    }
}
void tce_sortacc_4(double *in, double *out, int a, int b, int c, int d, int i, int j, int k, int l, double factor) {
    /* Same as tce_sort_4, above, but uses += instead of =. */
    int pos[4],len[4],ia,ib,j1,j2,j3,j4;
    ia = 0;
    len[0] = a;
    len[1] = b;
    len[2] = c;
    len[3] = d;
    for(j1 = 0; j1 < a; j1++) {
        pos[0] = j1;
        for(j2 = 0; j2 < b; j2++) {
            pos[1] = j2;
            for(j3 = 0; j3 < c; j3++) {
                pos[2] = j3;
                for(j4 = 0; j4 < d; j4++) {
                    pos[3] = j4;
                    ib = len[j]*len[k]*len[l]*pos[i]
                       +        len[k]*len[l]*pos[j]
                       +               len[l]*pos[k]
                       +                      pos[l];
                    out[ib] += in[ia++] * factor;
                }
            }
        }
    }
}

#define index_pair(i,j) ((i*(i+1))/2+j)
#define index_point(i,j,n) ((n*(n+1))/2-((n-i)*(n-i+1))/2-(n-j))


static int hash_comparator(const void *a_, const void *b_) {
    int *a = (int*)a_, *b = (int*)b_;
    return *a - *b;
}
/* Despite being called a "hash", it's a simple key-value lookup table.
 * hash points at an array of 2N+1 integers, with the following layout:
 * [ N | keys (N of these) | values (N of these) ] */
int tce_hash(int *hash, int key) {
    /* essentially the same as the one in tce_hash.F */
    int N = hash[0];

    int *elem = bsearch(&key, &hash[1], N, sizeof(key), hash_comparator);
    if(!elem) {
        int i;
        for(i = 0; i < NUM_VECTORS; i++)
            if(hash == vectors[i])
                break;
        char *name;
        if(i < NUM_VECTORS)
            name = vector_names[i];
        else
            name = "unknown";
        die("key %d not found in %s hash", key, name);
    }
    int o = elem - &hash[1];
    return hash[N+o+1];
}
struct tce_hash_v2_lookup {
    int i, j, k, l, pos;
    long offset;
} *tce_hash_v2_lookup;
int tce_hash_v2_count;
void tce_hash_v2_prepare(const int *hash, int M) {
    /* Build a lookup table such that tce_hash_v2() never has to iterate more
     * than M times.  (Said another way, it never has to move more than M
     * entries from the starting point. This lookup table replaces the
     * k_v2_alpha_offset array provided by nwchem. */
    int N = noa+nva;
    int entries = 0, skip = 0;
    long offset = 0;
    int i, j, k, l, pos, pos_u = N*N*N*N;
    i = j = k = l = pos = 0;
    printf("preparing tce_hash_v2 lookup table...\n");
    while(pos < pos_u) {
        if(k_spin_alpha[i]+k_spin_alpha[j] == k_spin_alpha[k]+k_spin_alpha[l]) {
            if((k_sym_alpha[i] ^ k_sym_alpha[j] ^ k_sym_alpha[k] ^ k_sym_alpha[l]) == irrep_v) {
                if(index_pair(j,i) >= index_pair(l,k)) {
//                    debug("block %d: %d,%d,%d,%d", entries, i,j,k,l);
                    entries++;
                }
            }
        }
        l++;
        pos++;
        if(l == N) {
            k++;
            if(k == N) {
                k = 0;
                j++;
                if(j == N) {
                    i++;
                    j = i;
                    pos += j * N * N;
                }
            }
            l = k;
            pos += l;
        }
    }
    if(entries % M)
        entries += entries % M;
    entries /= M;
    entries++;
    printf("%d entries allocated...", entries);
    tce_hash_v2_lookup = malloc(sizeof(struct tce_hash_v2_lookup)*entries);
    memset(tce_hash_v2_lookup, 0, sizeof(struct tce_hash_v2_lookup)*entries);
    struct tce_hash_v2_lookup *this = tce_hash_v2_lookup;
    i = j = k = l = pos = 0;
    while(pos < pos_u) {
        if(k_spin_alpha[i]+k_spin_alpha[j] == k_spin_alpha[k]+k_spin_alpha[l]) {
            if((k_sym_alpha[i] ^ k_sym_alpha[j] ^ k_sym_alpha[k] ^ k_sym_alpha[l]) == irrep_v) {
                if(index_pair(j,i) >= index_pair(l,k)) {
                    if(!(skip++ % M)) {
                        this->i = i;
                        this->j = j;
                        this->k = k;
                        this->l = l;
                        this->pos = pos;
                        this->offset = offset;
                        this++;
                    }
                    offset += k_range_alpha[i]*
                              k_range_alpha[j]*
                              k_range_alpha[k]*
                              k_range_alpha[l];
                }
            }
        }
        l++;
        pos++;
        if(l == N) {
            k++;
            if(k == N) {
                k = 0;
                j++;
                if(j == N) {
                    i++;
                    j = i;
                    pos += j * N * N;
                }
            }
            l = k;
            pos += l;
        }
    }
    this->i = i;
    this->j = j;
    this->k = k;
    this->l = l;
    this->pos = pos;
    this->offset = offset;
    this++;
    tce_hash_v2_count = this - tce_hash_v2_lookup;
    printf("%d used.\n", tce_hash_v2_count);
}
long tce_hash_v2(const int *hash, int key) {
    /* Equivalent to the subroutine in tce_hash.F.  Adapted to use C style
     * zero-based indexing, and rewritten loop iterator logic. */
    int N = noa+nva;
    long offset;
    int middle, min, max,
        i, j, k, l,
        i_start, j_start, k_start, l_start,
        pos, pos_l, pos_u;
    /* The purpose of this function is as follows.  Given an [i,j,k,l]
     * permutation index (key), generated by:
     * key = i*(noa+nva)^3 + j*(noa+nva)^2 + k*(noa+nva) + l
     * And a method of calculating the size for each position in the [i,j,k,l]
     * space:
     * size = k_range_alpha[i]*k_range_alpha[j]*k_range_alpha[k]*k_range_alpha[l]
     * And given some symmetry constraints (funky if-statements) to determine
     * whether the ijkl position is actually stored in the v2 array, find
     * the offset in the v2 array for the data block associated with the key.
     * (Or die with an error message if the symmetry constraints failed.)
     *
     * The tce_hash_v2_lookup[] table provides some starting points for this
     * search.  The code starts from the closest predecessor row, and loops
     * until it gets to the right index, adding offsets for valid blocks it
     * encounters along the way.
     */
    min = middle = 0;
    max = tce_hash_v2_count - 1;
    while(min < max) {
        middle = (max - min) / 2 + min + 1;
        if(tce_hash_v2_lookup[middle].pos <= key)
            min = middle;
        else
            max = middle - 1;
    }
    middle = min;
    if(middle == tce_hash_v2_count-1) {
        if(tce_hash_v2_lookup[middle+1].pos == key)
            return tce_hash_v2_lookup[middle+1].offset;
        die("tce_hash_v2: key %d not found", key);
    }
    die_if(tce_hash_v2_lookup[middle].pos > key || tce_hash_v2_lookup[middle+1].pos <= key);
    // middle is fixed
    i_start = tce_hash_v2_lookup[middle].i;
    j_start = tce_hash_v2_lookup[middle].j;
    k_start = tce_hash_v2_lookup[middle].k;
    l_start = tce_hash_v2_lookup[middle].l;

    offset      = tce_hash_v2_lookup[middle].offset;
    pos = pos_l = tce_hash_v2_lookup[middle].pos;
    pos_u       = tce_hash_v2_lookup[middle+1].pos;
    i = i_start;
    j = j_start;
    k = k_start;
    l = l_start;
    die_if(l + N*(k + N*(j + N*i)) != pos_l);

    while(pos < pos_u) {
        if(k_spin_alpha[i]+k_spin_alpha[j] == k_spin_alpha[k]+k_spin_alpha[l]) {
            if((k_sym_alpha[i] ^ k_sym_alpha[j] ^ k_sym_alpha[k] ^ k_sym_alpha[l]) == irrep_v) {
                if(index_pair(j,i) >= index_pair(l,k)) {
                    if(pos == key) break;
                    long size = k_range_alpha[i]*
                                k_range_alpha[j]*
                                k_range_alpha[k]*
                                k_range_alpha[l];

                    offset += size;
                }
            }
        }

        /* advance ijkl and pos1, taking advantage of j>=i and l>=k to skip some rounds */
        /* A j>=l skip is also possible, but is not currently implemented */
        l++;
        pos++;
        if(l == N) {
            k++;
            if(k == N) {
                k = 0;
                j++;
                if(j == N) {
                    i++;
                    if(i == N)
                        die("tce_hash_v2: failed to find value for key %d (end of permutation range)", key);
                    j = i;
                    pos += j * N * N;
                }
            }
            l = k;
            pos += l;
        }
    }
    if(pos >= pos_u)
        die("tce_hash_v2: failed to find value for key %d (end of bucket range)", key);
    return offset;
}

void tce_get_block_ind_i(double *in, double *out, int dim, int x2, int x1, int x4, int x3);
void tce_get_hash_block(double *in, double *out, int dim, int *hash, int key) {
    int offset;
    offset = tce_hash(hash, key);
    memcpy(out, &in[offset], dim*sizeof(*out));
    return;
}

void tce_add_hash_block(double *out, double *in, int dim, int *hash, int key) {
    int i, offset = tce_hash(hash, key);
    for(i = 0; i < dim; i++)
        out[offset+i] += in[i];
}

void write_tensor(char *filename, double *data, size_t len) {
    size_t record_size = 2*1024*1024/8;
    double d;
    FILE *f = fopen(filename, "w");
    fwrite(data, sizeof(*data), len, f);
    while(len % record_size) {
        d = 0;
        fwrite(&d, sizeof(d), 1, f);
        len++;
    }
    fclose(f);
}

void tce_get_block_ind_i(double *in, double *out, int dim, int x2, int x1, int x4, int x3) {
    /* Straight translation of the subroutine in get_block_ind.F, minus the unused "indexc" argument. */
    int i,j,k,l,ispin,size[4],sorder[4],porder[4],offset_alpha,key_alpha;
    int ix1,ix2,ix3,ix4,irow,icol;
    bool do_first,do_second;
    int first_h,second_h;
    double *tmp;

//    debug("tce_get_block_ind_i: x1=%d x2=%d x3=%d x4=%d", x1, x2, x3, x4);
    size[0] = k_range[x1];
    size[1] = k_range[x2];
    size[2] = k_range[x3];
    size[3] = k_range[x4];

    memset(out,0,dim*sizeof(*out));

    // v^{x3<x4}_{x1<x2} => ( x3 x1 | x4 x2 ) - ( x3 x2 | x4 x1 )
    do_first = do_second = 0;

    ispin = k_spin[x3] + k_spin[x4] + k_spin[x1] + k_spin[x2];
    if(ispin == 4) do_first = do_second = 1;
    if(ispin == 8) do_first = do_second = 1;
    if(k_spin[x3] == 1 && k_spin[x4] == 2 && k_spin[x1] == 1 && k_spin[x2] == 2)
        do_first = 1;
    if(k_spin[x3] == 2 && k_spin[x4] == 1 && k_spin[x1] == 2 && k_spin[x2] == 1)
        do_first = 1;
    if(k_spin[x3] == 2 && k_spin[x4] == 1 && k_spin[x1] == 1 && k_spin[x2] == 2)
        do_second = 1;
    if(k_spin[x3] == 1 && k_spin[x4] == 2 && k_spin[x1] == 2 && k_spin[x2] == 1)
        do_second = 1;

    //  do_first &&  do_second: v^{x3<x4}_{x1<x2} => ( x3 x1 | x4 x2 ) - ( x3 x2 | x4 x1 )
    //  do_first && !do_second: v^{x3<x4}_{x1<x2} => ( x3 x1 | x4 x2 ) -         0
    // !do_first &&  do_second: v^{x3<x4}_{x1<x2} =>         0         - ( x3 x2 | x4 x1 )

    if (do_first) {
        bool l31s,l42s,lp31p42;
        // first half
        i = max(k_b2am[x1], k_b2am[x3]);
        j = min(k_b2am[x1], k_b2am[x3]);
        k = max(k_b2am[x2], k_b2am[x4]);
        l = min(k_b2am[x2], k_b2am[x4]);
        irow = index_pair(i,j);
        icol = index_pair(k,l);
        if(irow >= icol) {
//            debug("...which translates to %d,%d,%d,%d", j-1, i-1, l-1, k-1);
            first_h = (((j-1)*(noa+nva) + i-1)*(noa+nva) + l-1)*(noa+nva) + k-1;
        } else {
//            debug("...which translates to %d,%d,%d,%d", l-1, k-1, j-1, i-1);
            first_h = (((l-1)*(noa+nva) + k-1)*(noa+nva) + j-1)*(noa+nva) + i-1;
        }
        key_alpha = first_h;
        // defining the order
        ix1 = k_b2am[x1];
        ix2 = k_b2am[x2];
        ix3 = k_b2am[x3];
        ix4 = k_b2am[x4];
        l31s = (ix3 < ix1);
        l42s = (ix4 < ix2);
        irow = l31s ? index_pair(ix1,ix3) : index_pair(ix3,ix1);
        icol = l42s ? index_pair(ix2,ix4) : index_pair(ix4,ix2);
        lp31p42 = (irow < icol);

        offset_alpha = tce_hash_v2(k_v2_alpha_offset, key_alpha);
        tmp = tce_double_malloc(dim);
        memcpy(tmp, in+offset_alpha, dim*sizeof(*out));

        if(     !lp31p42 && !l31s && !l42s) {
            // --- ( g3 g1 | g4 g2 )
            memcpy(sorder,(int[]){1,0,3,2},sizeof(sorder));
            memcpy(porder,(int[]){3,1,2,0},sizeof(porder));
        } else if(!lp31p42 && !l31s && l42s) {
            // --- ( g3 g1 | g2 g4 )
            memcpy(sorder,(int[]){3,1,0,2},sizeof(sorder));
            memcpy(porder,(int[]){3,0,2,1},sizeof(porder));
        } else if(!lp31p42 && l31s && !l42s) {
            // --- ( g1 g3 | g4 g2 )
            memcpy(sorder,(int[]){1,3,2,0},sizeof(sorder));
            memcpy(porder,(int[]){2,1,3,0},sizeof(porder));
        } else if(!lp31p42 && l31s && l42s) {
            // --- ( g1 g3 | g2 g4 )
            memcpy(sorder,(int[]){3,1,2,0},sizeof(sorder));
            memcpy(porder,(int[]){2,0,3,1},sizeof(porder));
        } else if(lp31p42 && !l31s && !l42s) {
            // --- ( g4 g2 | g3 g1 )
            memcpy(sorder,(int[]){0,2,1,3},sizeof(sorder));
            memcpy(porder,(int[]){1,3,0,2},sizeof(porder));
        } else if(lp31p42 && !l31s && l42s) {
            // --- ( g2 g4 | g3 g1 )
            memcpy(sorder,(int[]){0,2,3,1},sizeof(sorder));
            memcpy(porder,(int[]){1,2,0,3},sizeof(porder));
        } else if(lp31p42 && l31s && !l42s) {
            // --- ( g4 g2 | g1 g3 )
            memcpy(sorder,(int[]){2,0,1,3},sizeof(sorder));
            memcpy(porder,(int[]){0,3,1,2},sizeof(porder));
        } else if(lp31p42 && l31s && l42s) {
            // --- ( g2 g4 | g1 g3 )
            memcpy(sorder,(int[]){2,0,3,1},sizeof(sorder));
            memcpy(porder,(int[]){0,2,1,3},sizeof(porder));
        } else
            die("unhandled permutation");
//        debug("...in order {%d,%d,%d,%d}", porder[0], porder[1], porder[2], porder[3]);
        tce_sortacc_4(tmp,out,
                      size[sorder[0]],size[sorder[1]],size[sorder[2]],size[sorder[3]],
                           porder[0] ,     porder[1] ,     porder[2] ,     porder[3] ,
                      1.0);

        free(tmp);
    } // !spin cases

    if(do_second) {
        bool l32s,l41s,lp32p41;
        // second half

        i = max(k_b2am[x2], k_b2am[x3]);
        j = min(k_b2am[x2], k_b2am[x3]);
        k = max(k_b2am[x1], k_b2am[x4]);
        l = min(k_b2am[x1], k_b2am[x4]);
        irow = index_pair(i,j);
        icol = index_pair(k,l);
        if(irow >= icol) {
//            debug("...minus %d,%d,%d,%d", j-1, i-1, l-1, k-1);
            second_h = (((j-1)*(noa+nva) + i-1)*(noa+nva) + l-1)*(noa+nva) + k-1;
        } else {
//            debug("...minus %d,%d,%d,%d", l-1, k-1, j-1, i-1);
            second_h = (((l-1)*(noa+nva) + k-1)*(noa+nva) + j-1)*(noa+nva) + i-1;
        }
        key_alpha = second_h;
        // defining the order
        ix2 = k_b2am[x2];
        ix1 = k_b2am[x1];
        ix3 = k_b2am[x3];
        ix4 = k_b2am[x4];
        l32s = (ix3 < ix2);
        l41s = (ix4 < ix1);
        irow = l32s ? index_pair(ix2,ix3) : index_pair(ix3,ix2);
        icol = l41s ? index_pair(ix1,ix4) : index_pair(ix4,ix1);
        lp32p41 = (irow < icol);

        offset_alpha = tce_hash_v2(k_v2_alpha_offset, key_alpha);
        tmp = tce_double_malloc(dim);
        /* call ga_get(d_v2orb,offset_alpha+1,offset_alpha+size,1,1,dbl_mb(k_a),1) */
        memcpy(tmp, in+offset_alpha, dim*sizeof(*out));

        if(!lp32p41 && !l32s && !l41s) {
            // --- ( g3 g2 | g4 g1 )
            memcpy(sorder,(int[]){0,3,1,2},sizeof(sorder)); memcpy(porder,(int[]){3,1,0,2},sizeof(porder));
        } else if(!lp32p41 && !l32s && l41s) {
            // --- ( g3 g2 | g1 g4 )
            memcpy(sorder,(int[]){1,0,3,2},sizeof(sorder)); memcpy(porder,(int[]){3,0,1,2},sizeof(porder));
        } else if(!lp32p41 && l32s && !l41s) {
            // --- ( g2 g3 | g4 g1 )
            memcpy(sorder,(int[]){0,3,2,1},sizeof(sorder)); memcpy(porder,(int[]){2,1,0,3},sizeof(porder));
        } else if(!lp32p41 && l32s && l41s) {
            // --- ( g2 g3 | g1 g4 )
            memcpy(sorder,(int[]){0,3,2,1},sizeof(sorder)); memcpy(porder,(int[]){2,0,1,3},sizeof(porder));
        } else if(lp32p41 && !l32s && !l41s) {
            // --- ( g4 g1 | g3 g2 )
            memcpy(sorder,(int[]){1,2,0,3},sizeof(sorder)); memcpy(porder,(int[]){1,3,2,0},sizeof(porder));
        } else if(lp32p41 && l32s && !l41s) {
            // --- ( g4 g1 | g2 g3 )
            memcpy(sorder,(int[]){2,1,0,3},sizeof(sorder)); memcpy(porder,(int[]){0,3,2,1},sizeof(porder));
        } else if(lp32p41 && !l32s && l41s) {
            // --- ( g1 g4 | g3 g2 )
            memcpy(sorder,(int[]){1,0,3,2},sizeof(sorder)); memcpy(porder,(int[]){1,2,3,0},sizeof(porder));
        } else if(lp32p41 && l32s && l41s) {
            // --- ( g1 g4 | g2 g3 )
            memcpy(sorder,(int[]){2,1,0,3},sizeof(sorder)); memcpy(porder,(int[]){0,2,3,1},sizeof(porder));
        } else
            die("unhandled permutation");
//        debug("...in order {%d,%d,%d,%d}", porder[0], porder[1], porder[2], porder[3]);
        tce_sortacc_4(tmp,out,
                      size[sorder[0]],size[sorder[1]],size[sorder[2]],size[sorder[3]],
                           porder[0] ,     porder[1] ,     porder[2] ,     porder[3] ,
                      -1.0);

        free(tmp);
    } // !spin cases
//    dump_block(out, size[0]*size[1]*size[2]*size[3]);
}
void dump_block(double *data, int len) {
    int i;
    for(i = 0; i < len; i++)
        debug("  [%d] = %g", i, data[i]);
}
void dump_tensor(double *data, int *hash, int dims, int *sizes, int *starts) {
    int blockid, i, iter[dims], N = hash[0];
    int *hash_keys = &hash[1], *hash_values = &hash[N+1];
    for(blockid = 0; blockid < N; blockid++) {
        int len = 1;
        int key = hash_keys[blockid];
        int offset = hash_values[blockid];
        for(i = dims - 1; i >= 0; i--) {
            int tmp = key % sizes[i];
            key -= tmp;
            key /= sizes[i];
            iter[i] = tmp + starts[i];
        }
        for(i = 0; i < dims; i++)
            len *= k_range[iter[i]];
        printf("block len %d at ", len);
        for(i = 0; i < dims; i++)
            printf("%s%d", i ? "," : "", iter[i]);
        printf(" (index %d):\n", hash_keys[blockid]);
        dump_block(data+offset,len);
    }
}
