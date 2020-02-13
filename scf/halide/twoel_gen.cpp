#include <Halide.h>
#include <stdio.h>
using namespace Halide;

#define pi 3.141592653589793
#define tol2e 0.000001

class twoel_generator : public Halide::Generator<twoel_generator> {
public:

    /* input scalars */
    Input<double> delo2{"delo2"};
    Input<double> delta{"delta"};
    Input<double> rdelta{"rdelta"};

    /* output scalars */
    Output<double> rv{"rv"};

    /* input matrices */
    Input<Buffer<double>> fm_in{"fm_in", 2};
    Input<Buffer<double>> g_fock_in_in{"g_fock_in", 2};
    Input<Buffer<double>> g_dens_in{"g_dens_in", 2};

    /* partial g() precomputation data */
    Input<Buffer<double>> g_precalc{"g_precalc", 3};

    /* output matrix */
    Output<Buffer<double>> g_fock_out{"g_fock_out", 2};

    Var i{"i"}, j{"j"}, k{"k"}, l{"l"};

    void generate() {
        /* clamp all inputs, to prevent out-of-bounds errors from odd tile sizes and such */
        Func fm, g_fock_in, g_dens;

        fm        = BoundaryConditions::constant_exterior(fm_in       , 0);
        g_fock_in = BoundaryConditions::constant_exterior(g_fock_in_in, 0);
        g_dens    = BoundaryConditions::constant_exterior(g_dens_in   , 0);

        Func ex{"ex"}, denom{"denom"}, fac4d{"fac4d"};

        ex(i,j,k,l) = g_precalc(i,j,0) * g_precalc(k,l,0);
        denom(i,j,k,l)  = g_precalc(i,j,1) * g_precalc(k,l,1) * sqrt(g_precalc(i,j,1) + g_precalc(k,l,1));
        fac4d(i,j,k,l)  = g_precalc(i,j,1) * g_precalc(k,l,1) /     (g_precalc(i,j,1) + g_precalc(k,l,1));

        Func rnorm2{"rnorm2"}, x2{"x2"}, y2{"y2"}, z2{"z2"}, rpq2{"rpq2"};
        x2(i,j)     = g_precalc(i,j,2);
        y2(i,j)     = g_precalc(i,j,3);
        z2(i,j)     = g_precalc(i,j,4);
        rnorm2(i,j) = g_precalc(i,j,5);
        rpq2(i,j,k,l) =
              (x2(i,j) - x2(k,l)) * (x2(i,j) - x2(k,l))
            + (y2(i,j) - y2(k,l)) * (y2(i,j) - y2(k,l))
            + (z2(i,j) - z2(k,l)) * (z2(i,j) - z2(k,l));

        Func f0t{"f0t"}, f0n{"f0n"}, f0x{"f0x"}, f0val{"f0val"};
        f0t(i,j,k,l) = fac4d(i,j,k,l) * rpq2(i,j,k,l);
        f0n(i,j,k,l) = clamp(cast<int>((f0t(i,j,k,l) + delo2) * rdelta), fm_in.dim(0).min(), fm_in.dim(0).max());
        f0x(i,j,k,l) = delta * f0n(i,j,k,l) - f0t(i,j,k,l);
        f0val(i,j,k,l) = select(f0t(i,j,k,l) >= Expr(28.0),
             Expr(0.88622692545276) / sqrt(f0t(i,j,k,l)),
                                             fm(f0n(i,j,k,l),0)
             + f0x(i,j,k,l) *               (fm(f0n(i,j,k,l),1)
             + f0x(i,j,k,l) * Expr(0.5) *   (fm(f0n(i,j,k,l),2)
             + f0x(i,j,k,l) * Expr(1./3.) * (fm(f0n(i,j,k,l),3)
             + f0x(i,j,k,l) * Expr(0.25) *   fm(f0n(i,j,k,l),4)))));

        Func gg{"gg"};
        gg(i, j, k, l) = (Expr(2.00) * Expr(pow(pi, 2.50)) / denom(i,j,k,l)) * ex(i,j,k,l) * f0val(i,j,k,l) * rnorm2(i,j) * rnorm2(k,l);

        RDom r(0,g_fock_in_in.width(), 0, g_fock_in_in.height());

        Func g_fock{"g_fock"};
        g_fock(i,j) = g_fock_in(i,j)
                    + sum(gg(i,j,r.x,r.y) * g_dens(r.x,r.y))
                    - sum(gg(i,r.x,j,r.y) * g_dens(r.x,r.y)) * Expr(0.5);

        g_fock_out(i,j) = g_fock(i,j);
        rv() = sum(g_fock(r.x,r.y) * g_dens(r.x,r.y)) * Expr(0.5);

#ifdef TRACING
        gg.trace_realizations();
        gg.trace_loads();
        gg.trace_stores();
        g_dens_in.trace_loads();
        fm_in.trace_loads();
        x_in.trace_loads();
        y_in.trace_loads();
        z_in.trace_loads();
        expnt_in.trace_loads();
        rnorm_in.trace_loads();
        g_fock_in_in.trace_loads();
        g_fock_out.trace_stores();
#endif /* TRACING */
    }

    void schedule() {
        if(auto_schedule) {
            /* input scalars */
            delta.set_estimate(0.014);
            delo2.set_estimate(0.007);
            rdelta.set_estimate(71.4285);
            /* input matrices */
            g_fock_in_in.dim(0).set_estimate(0, 1024);
            g_fock_in_in.dim(1).set_estimate(0, 1024);
            g_precalc.dim(0).set_estimate(0, 1024);
            g_precalc.dim(1).set_estimate(0, 1024);
            g_precalc.dim(2).set_estimate(6, 6);
            g_dens_in.dim(0).set_estimate(0, 1024);
            g_dens_in.dim(1).set_estimate(0, 1024);
            fm_in.dim(0).set_estimate(0, 2001);
            fm_in.dim(1).set_estimate(0, 5);
            /* output matrix */
            g_fock_out.set_estimate(i, 0, 1024);
            g_fock_out.set_estimate(j, 0, 1024);
        } else {
            /* insert schedule here */
        }
    }
};

HALIDE_REGISTER_GENERATOR(twoel_generator, twoel);
