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

    /* input vectors */
    Input<Buffer<double>> expnt_unclamped{"expnt", 1};
    Input<Buffer<double>> rnorm_unclamped{"rnorm", 1};
    Input<Buffer<double>> x_unclamped{"x", 1};
    Input<Buffer<double>> y_unclamped{"y", 1};
    Input<Buffer<double>> z_unclamped{"z", 1};

    /* input matrices */
    Input<Buffer<double>> fm_unclamped{"fm", 2};
    Input<Buffer<double>> g_fock_in_unclamped{"g_fock_in", 2};
    Input<Buffer<double>> g_dens_unclamped{"g_dens", 2};

    /* output matrix */
    Output<Buffer<double>> g_fock_out{"g_fock_out", 2};

    Var i{"i"}, j{"j"}, k{"k"}, l{"l"};

    void generate() {
        /* clamp all inputs, to prevent out-of-bounds errors from odd tile sizes and such */
        Func expnt{"delo2_clamped"}, rnorm{"delo2_clamped"},
             x{"delo2_clamped"}, y{"delo2_clamped"}, z{"delo2_clamped"},
             fm{"delo2_clamped"}, g_fock_in{"delo2_clamped"}, g_dens{"delo2_clamped"};

        expnt     = BoundaryConditions::constant_exterior(expnt_unclamped    , 0);
        rnorm     = BoundaryConditions::constant_exterior(rnorm_unclamped    , 0);
        x         = BoundaryConditions::constant_exterior(x_unclamped        , 0);
        y         = BoundaryConditions::constant_exterior(y_unclamped        , 0);
        z         = BoundaryConditions::constant_exterior(z_unclamped        , 0);
        fm        = BoundaryConditions::constant_exterior(fm_unclamped       , 0);
        g_fock_in = BoundaryConditions::constant_exterior(g_fock_in_unclamped, 0);
        g_dens    = BoundaryConditions::constant_exterior(g_dens_unclamped   , 0);

        Func dx{"dx"}, dy{"dy"}, dz{"dz"}, r2{"r2"}, expnt2{"expnt2"}, expnt_inv{"expnt_inv"};

        dx(i,j) = x(i) - x(j);
        dy(i,j) = y(i) - y(j);
        dz(i,j) = z(i) - z(j);

        r2(i,j) = dx(i,j) * dx(i,j) + dy(i,j) * dy(i,j) + dz(i,j) * dz(i,j);

        expnt2(i,j)     = expnt(i) + expnt(j);
        expnt_inv(i,j) = Expr(1.0) / expnt2(i,j);

        Func fac2{"fac2"}, ex_arg{"ex_arg"}, ex{"ex"}, denom{"denom"}, fac4d{"fac4d"};
        fac2(i,j) = expnt(i) * expnt(j) * expnt_inv(i,j);
        ex_arg(i,j,k,l) = -fac2(i,j) * r2(i,j) - fac2(k,l) * r2(k,l);
        ex(i,j,k,l) = select(ex_arg(i,j,k,l) < Expr(-37.0), Expr(0.0), exp(ex_arg(i,j,k,l)));
        denom(i,j,k,l)  = expnt2(i,j) * expnt2(k,l) * sqrt(expnt2(i,j) + expnt2(k,l));
        fac4d(i,j,k,l)  = expnt2(i,j) * expnt2(k,l) /     (expnt2(i,j) + expnt2(k,l));

        Func x2{"x2"}, y2{"y2"}, z2{"z2"}, rpq2{"rpq2"};
        x2(i,j)     = (x(i) * expnt(i) + x(j) * expnt(j)) * expnt_inv(i,j);
        y2(i,j)     = (y(i) * expnt(i) + y(j) * expnt(j)) * expnt_inv(i,j);
        z2(i,j)     = (z(i) * expnt(i) + z(j) * expnt(j)) * expnt_inv(i,j);
        rpq2(i,j,k,l) =
              (x2(i,j) - x2(k,l)) * (x2(i,j) - x2(k,l))
            + (y2(i,j) - y2(k,l)) * (y2(i,j) - y2(k,l))
            + (z2(i,j) - z2(k,l)) * (z2(i,j) - z2(k,l));

        Func f0t{"f0t"}, f0n{"f0n"}, f0x{"f0x"}, f0val{"f0val"};
        f0t(i,j,k,l) = fac4d(i,j,k,l) * rpq2(i,j,k,l);
        f0n(i,j,k,l) = clamp(cast<int>((f0t(i,j,k,l) + delo2) * rdelta), fm_unclamped.dim(0).min(), fm_unclamped.dim(0).max());
        f0x(i,j,k,l) = delta * f0n(i,j,k,l) - f0t(i,j,k,l);
        f0val(i,j,k,l) = select(f0t(i,j,k,l) >= Expr(28.0),
             Expr(0.88622692545276) / sqrt(f0t(i,j,k,l)),
                                             fm(f0n(i,j,k,l),0)
             + f0x(i,j,k,l) *               (fm(f0n(i,j,k,l),1)
             + f0x(i,j,k,l) * Expr(0.5) *   (fm(f0n(i,j,k,l),2)
             + f0x(i,j,k,l) * Expr(1./3.) * (fm(f0n(i,j,k,l),3)
             + f0x(i,j,k,l) * Expr(0.25) *   fm(f0n(i,j,k,l),4)))));

        Func gg{"gg"};
        gg(i, j, k, l) = (Expr(2.00) * Expr(pow(pi, 2.50)) / denom(i,j,k,l)) * ex(i,j,k,l) * f0val(i,j,k,l) * rnorm(i) * rnorm(j) * rnorm(k) * rnorm(l);

        RDom r(0,g_fock_in_unclamped.width(), 0, g_fock_in_unclamped.height());

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
        g_dens_unclamped.trace_loads();
        fm_unclamped.trace_loads();
        x_unclamped.trace_loads();
        y_unclamped.trace_loads();
        z_unclamped.trace_loads();
        expnt_unclamped.trace_loads();
        rnorm_unclamped.trace_loads();
        g_fock_in_unclamped.trace_loads();
        g_fock_out.trace_stores();
#endif /* TRACING */
    }

    void schedule() {
        if(auto_schedule) {
            /* input scalars */
            delta.set_estimate(0.014);
            delo2.set_estimate(0.007);
            rdelta.set_estimate(71.4285);
            /* input vectors */
            x_unclamped.dim(0).set_estimate(0, 1024);
            y_unclamped.dim(0).set_estimate(0, 1024);
            z_unclamped.dim(0).set_estimate(0, 1024);
            expnt_unclamped.dim(0).set_estimate(0, 1024);
            rnorm_unclamped.dim(0).set_estimate(0, 1024);
            /* input matrices */
            g_fock_in_unclamped.dim(0).set_estimate(0, 1024);
            g_fock_in_unclamped.dim(1).set_estimate(0, 1024);
            g_dens_unclamped.dim(0).set_estimate(0, 1024);
            g_dens_unclamped.dim(1).set_estimate(0, 1024);
            fm_unclamped.dim(0).set_estimate(0, 2001);
            fm_unclamped.dim(1).set_estimate(0, 5);
            /* output matrix */
            g_fock_out.set_estimate(i, 0, 1024);
            g_fock_out.set_estimate(j, 0, 1024);
        } else {
            /* insert schedule here */
        }
    }
};

HALIDE_REGISTER_GENERATOR(twoel_generator, twoel);
