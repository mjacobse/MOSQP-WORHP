#include "TestProblems.hpp"
#include "../problem_formulation/MatrixStructure.hpp"
#include "../problem_formulation/MONLP.hpp"
#include <cstddef>


namespace test_problems
{

using mosqp::MatrixStructure;

Jo1::Jo1() : MONLP(
    5,  // number of variables
    3,  // number of constraints
    2,  // number of objectives

    { MatrixStructure({ 1, 2, 3, 4, 5 }),                              // structure of df1
      MatrixStructure({ 1, 2, 3, 4, 5 }) },                            // structure of df2
    MatrixStructure({ 1, 2, 3, 1, 2, 3, 1, 2, 3, 1, 2, 3, 1, 2, 3 },
                    { 1, 1, 1, 2, 2, 2, 3, 3, 3, 4, 4, 4, 5, 5, 5 }),  // structure of dg
    { MatrixStructure({ 1, 2, 3, 4, 5 },{ 1, 2, 3, 4, 5 }),            // structure of hm1
      MatrixStructure({ 1, 2, 3, 4, 5 },{ 1, 2, 3, 4, 5 }) },          // structure of hm2
    true, true, true,                                 // got all user derivatives

    { NEG_INF, NEG_INF, NEG_INF, NEG_INF, NEG_INF },  // lower bounds on x
    { POS_INF, POS_INF, POS_INF, POS_INF, POS_INF },  // upper bounds on x
    { 0.0, 0.0, NEG_INF },                            // lower bounds on g
    { 0.0, 0.0, 0.0 })                                // upper bounds on g
{
}

double Jo1::EvalF_impl(double const *const x, size_t const objective_index) const
{
    switch (objective_index)
    {
    case 0:
        return x[0] * x[0] + x[1] * x[1] + x[2] * x[2] + x[3] * x[3] + x[4] * x[4];
    case 1:
        return 3.0 * x[0] + 2.0 * x[1] - x[2] / 3.0 + (x[3] - x[4]) * (x[3] - x[4]) * (x[3] - x[4]) / 100.0;
    }
    return 0;
}

void Jo1::EvalG_impl(double const *const x, double *const g) const
{
    g[0] = x[0] + 2.0 * x[1] - x[2] - x[3] / 2.0 + x[4] - 2;
    g[1] = 4 * x[0] - 2.0 * x[1] + 0.8 * x[2] + 0.6 * x[3] + 0.5 * x[4] * x[4];
    g[2] = x[0] * x[0] + x[1] * x[1] + x[2] * x[2] + x[3] * x[3] + x[4] * x[4] - 10.0;
}

void Jo1::EvalDF_impl(double const *const x, double *const df, size_t const objective_index) const
{
    switch (objective_index)
    {
    case 0:
        df[0] = 2.0 * x[0];
        df[1] = 2.0 * x[1];
        df[2] = 2.0 * x[2];
        df[3] = 2.0 * x[3];
        df[4] = 2.0 * x[4];
        break;
    case 1:
        df[0] = 3.0;
        df[1] = 2.0;
        df[2] = - 1.0 / 3.0;
        df[3] = + 0.03 * (x[3] - x[4]) * (x[3] - x[4]);
        df[4] = - 0.03 * (x[3] - x[4]) * (x[3] - x[4]);
        break;
    }
}

void Jo1::EvalDG_impl(double const *const x, double *const dg) const
{
    dg[0] = 1.0;
    dg[1] = 4.0;
    dg[2] = 2.0 * x[0];
    dg[3] = 2.0;
    dg[4] = -2.0;
    dg[5] = 2.0 * x[1];
    dg[6] = -1.0;
    dg[7] = 0.8;
    dg[8] = 2.0 * x[2];
    dg[9] = -0.5;
    dg[10] = 0.6;
    dg[11] = 2.0 * x[3];
    dg[12] = 1.0;
    dg[13] = 1.0 * x[4];
    dg[14] = 2.0 * x[4];
}

void Jo1::EvalD2F_impl(double const *const x, double *const d2f, size_t const objective_index) const
{
    switch (objective_index)
    {
    case 0:
        d2f[0] = 2.0;
        d2f[1] = 2.0;
        d2f[2] = 2.0;
        d2f[3] = 2.0;
        d2f[4] = 2.0;
        break;
    case 1:
        d2f[0] = 0.0;
        d2f[1] = 0.0;
        d2f[2] = 0.0;
        d2f[3] = 0.06 * (x[3] - x[4]);
        d2f[4] = 0.06 * (x[3] - x[4]);
        break;
    }
}

void Jo1::EvalD2G_impl(double const *x, double const *mu, double *d2g, size_t const objective_index) const
{
    d2g[0] = 2.0 * mu[2];
    d2g[1] = 2.0 * mu[2];
    d2g[2] = 2.0 * mu[2];
    d2g[3] = 2.0 * mu[2];
    d2g[4] = 1.0 * mu[1] + 2.0 * mu[2];
}

} // namespace test_problems
