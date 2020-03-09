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
    Input<Buffer<double>> expnt_in{"expnt_in", 1};
    Input<Buffer<double>> rnorm_in{"rnorm_in", 1};
    Input<Buffer<double>> x_in{"x_in", 1};
    Input<Buffer<double>> y_in{"y_in", 1};
    Input<Buffer<double>> z_in{"z_in", 1};

    /* input matrices */
    Input<Buffer<double>> fm_in{"fm_in", 2};
    Input<Buffer<double>> g_fock_in_in{"g_fock_in", 2};
    Input<Buffer<double>> g_dens_in{"g_dens_in", 2};

    /* output matrix */
    Output<Buffer<double>> g_fock_out{"g_fock_out", 2};

    Var i{"i"}, j{"j"}, k{"k"}, l{"l"};

    void generate() {
        Expr nbfn = g_fock_in_in.height();

        /* clamp all inputs, to prevent out-of-bounds errors from odd tile sizes and such */
        Func expnt, rnorm, x, y, z, fm, g_fock_in, g_dens;

        expnt     = BoundaryConditions::constant_exterior(expnt_in    , 0);
        rnorm     = BoundaryConditions::constant_exterior(rnorm_in    , 0);
        x         = BoundaryConditions::constant_exterior(x_in        , 0);
        y         = BoundaryConditions::constant_exterior(y_in        , 0);
        z         = BoundaryConditions::constant_exterior(z_in        , 0);
        fm        = BoundaryConditions::constant_exterior(fm_in       , 0);
        g_fock_in = BoundaryConditions::constant_exterior(g_fock_in_in, 0);
        g_dens    = BoundaryConditions::constant_exterior(g_dens_in   , 0);

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
        f0n(i,j,k,l) = clamp(cast<int>((f0t(i,j,k,l) + delo2) * rdelta), fm_in.dim(0).min(), fm_in.dim(0).max());
        f0x(i,j,k,l) = delta * f0n(i,j,k,l) - f0t(i,j,k,l);
        f0val(i,j,k,l) = select(f0t(i,j,k,l) >= Expr(28.0),
             Expr(0.88622692545276) / sqrt(f0t(i,j,k,l)),
                                             fm(f0n(i,j,k,l),0)
             + f0x(i,j,k,l) *               (fm(f0n(i,j,k,l),1)
             + f0x(i,j,k,l) * Expr(0.5) *   (fm(f0n(i,j,k,l),2)
             + f0x(i,j,k,l) * Expr(1./3.) * (fm(f0n(i,j,k,l),3)
             + f0x(i,j,k,l) * Expr(0.25) *   fm(f0n(i,j,k,l),4)))));

        Func g_fock{"g_fock"};
        g_fock(i,j) = g_fock_in(i,j);

        Func g{"g"};
        g(i,j,k,l) = (Expr(2.00) * Expr(pow(pi, 2.50)) / denom(i,j,k,l)) * ex(i,j,k,l) * f0val(i,j,k,l) * rnorm(i) * rnorm(j) * rnorm(k) * rnorm(l);
        RVar gi, gj, gk, gl;


        // symmetry domain A: 3-way symmetry.  i == j == k == l.  1D iteration space
        RDom g_triple_symmetry_dom(0, nbfn);
        gi = g_triple_symmetry_dom[0];
        Expr g_triple_symmetry = g(gi, gi, gi, gi);
        g_fock(gi,gi) += g_triple_symmetry * (Expr(0.5) * g_dens(gi,gi));

        // symmetry domain B: 2-way symmetry.  i == j, k == l, i < k.  2D iteration space
        RDom g_double_symmetry_dom(0, nbfn, 0, nbfn);
        gi = g_double_symmetry_dom[0]; gk = g_double_symmetry_dom[1];
        g_double_symmetry_dom.where(gi < gk);
        Expr g_double_symmetry = g(gi, gi, gk, gk);
        g_fock(gi,gi) += g_double_symmetry * (g_dens(gk,gk));
        g_fock(gi,gk) += g_double_symmetry * (Expr(-0.5) * g_dens(gi,gk));
        g_fock(gk,gi) += g_double_symmetry * (Expr(-0.5) * g_dens(gk,gi));
        g_fock(gk,gk) += g_double_symmetry * (g_dens(gi,gi));

        // symmetry domain C: pairwise symmetry.  i < j, i == k, j == l.  2D iteration space
        RDom g_pairwise_symmetry_dom(0, nbfn, 0, nbfn);
        gi = g_pairwise_symmetry_dom[0]; gj = g_pairwise_symmetry_dom[1];
        g_pairwise_symmetry_dom.where(gi < gj);
        Expr g_pairwise_symmetry = g(gi, gj, gi, gj);
        g_fock(gi,gi) += g_pairwise_symmetry * (Expr(-0.5) * g_dens(gj,gj));
        g_fock(gi,gj) += g_pairwise_symmetry * (g_dens(gi,gj) + g_dens(gj,gi));
        g_fock(gi,gj) += g_pairwise_symmetry * (Expr(-0.5) * g_dens(gj,gi));
        g_fock(gj,gi) += g_pairwise_symmetry * (g_dens(gi,gj) + g_dens(gj,gi));
        g_fock(gj,gi) += g_pairwise_symmetry * (Expr(-0.5) * g_dens(gi,gj));
        g_fock(gj,gj) += g_pairwise_symmetry * (Expr(-0.5) * g_dens(gi,gi));

        // symmetry domain D: single symmetry.  i == j, k < l.  3D iteration space
        RDom g_single_symmetry_dom(0, nbfn, 0, nbfn, 0, nbfn);
        gi = g_single_symmetry_dom[0]; gk = g_single_symmetry_dom[1]; gl = g_single_symmetry_dom[2];
        g_single_symmetry_dom.where(gk < gl);
        Expr g_single_symmetry = g(gi, gi, gk, gl);
        g_fock(gi,gi) += g_single_symmetry * (g_dens(gk,gl) + g_dens(gl,gk));
        g_fock(gi,gk) += g_single_symmetry * (Expr(-0.5) * g_dens(gi,gl));
        g_fock(gi,gl) += g_single_symmetry * (Expr(-0.5) * g_dens(gi,gk));
        g_fock(gk,gi) += g_single_symmetry * (Expr(-0.5) * g_dens(gl,gi));
        g_fock(gk,gl) += g_single_symmetry * (g_dens(gi,gi));
        g_fock(gl,gi) += g_single_symmetry * (Expr(-0.5) * g_dens(gk,gi));
        g_fock(gl,gk) += g_single_symmetry * (g_dens(gi,gi));

        // symmetry domain E: no symmetry.  i < j, i <= k, k < l, j != l.  4D iteration space
        RDom g_no_symmetry_dom(0, nbfn, 0, nbfn, 0, nbfn, 0, nbfn);
        gi = g_no_symmetry_dom[0]; gj = g_no_symmetry_dom[1]; gk = g_no_symmetry_dom[2]; gl = g_no_symmetry_dom[3];
        g_no_symmetry_dom.where(gi < gj && gi <= gk && gk < gl && gi * nbfn + gj < gk * nbfn + gl);
        Expr g_no_symmetry = g(gi, gj, gk, gl);
        g_fock(gi,gj) += g_no_symmetry * (g_dens(gk,gl) + g_dens(gl,gk));
        g_fock(gi,gk) += g_no_symmetry * (Expr(-0.5) * g_dens(gj,gl));
        g_fock(gi,gl) += g_no_symmetry * (Expr(-0.5) * g_dens(gj,gk));
        g_fock(gj,gi) += g_no_symmetry * (g_dens(gk,gl) + g_dens(gl,gk));
        g_fock(gj,gk) += g_no_symmetry * (Expr(-0.5) * g_dens(gi,gl));
        g_fock(gj,gl) += g_no_symmetry * (Expr(-0.5) * g_dens(gi,gk));
        g_fock(gk,gi) += g_no_symmetry * (Expr(-0.5) * g_dens(gl,gj));
        g_fock(gk,gj) += g_no_symmetry * (Expr(-0.5) * g_dens(gl,gi));
        g_fock(gk,gl) += g_no_symmetry * (g_dens(gi,gj) + g_dens(gj,gi));
        g_fock(gl,gi) += g_no_symmetry * (Expr(-0.5) * g_dens(gk,gj));
        g_fock(gl,gj) += g_no_symmetry * (Expr(-0.5) * g_dens(gk,gi));
        g_fock(gl,gk) += g_no_symmetry * (g_dens(gi,gj) + g_dens(gj,gi));

        RDom r(0, nbfn, 0, nbfn);

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
            /* input vectors */
            x_in.dim(0).set_estimate(0, 1024);
            y_in.dim(0).set_estimate(0, 1024);
            z_in.dim(0).set_estimate(0, 1024);
            expnt_in.dim(0).set_estimate(0, 1024);
            rnorm_in.dim(0).set_estimate(0, 1024);
            /* input matrices */
            g_fock_in_in.dim(0).set_estimate(0, 1024);
            g_fock_in_in.dim(1).set_estimate(0, 1024);
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
