#include "TestProblems.hpp"
#include "../problem_formulation/MatrixStructure.hpp"
#include "../problem_formulation/MONLP.hpp"
#include <cstddef>


namespace test_problems
{

using mosqp::MatrixStructure;

Jo2::Jo2() : MONLP(
    3,  // number of variables
    2,  // number of constraints
    2,  // number of objectives

    { MatrixStructure({ 1, 2, 3 }),                   // structure of df1
      MatrixStructure({ 1, 2, 3 }) },                 // structure of df2
    MatrixStructure({ 1, 2, 1, 2, 1, 2 },
                    { 1, 1, 2, 2, 3, 3 }),            // structure of dg
    { MatrixStructure({ 1, 2, 3 },{ 1, 2, 3 }),       // structure of hm1
      MatrixStructure({ 1, 2, 3 },{ 1, 2, 3 }) },     // structure of hm2
    true, true, true,                                 // got all user derivatives

    { NEG_INF, NEG_INF, NEG_INF },  // lower bounds on x
    { POS_INF, POS_INF, POS_INF },  // upper bounds on x
    { 0.0, NEG_INF },               // lower bounds on g
    { 0.0, 0.0 })                   // upper bounds on g
{
}

double Jo2::EvalF_impl(double const *const x, size_t const objective_index) const
{
    switch (objective_index)
    {
    case 0:
        return (x[0] - 1.0) * (x[0] - 1.0) + (x[1] - 2.0) * (x[1] - 2.0) + (x[2] - 3.0) * (x[2] - 3.0);
    case 1:
        return 3.0 * std::sqrt(x[0] * x[0] + 1.0) + 1.0 / 3.0 * std::sqrt(x[1] * x[1] + 1.0) + 0.01 * (x[2] - 1.0) * (x[2] - 1.0);
    }
    return 0;
}

void Jo2::EvalG_impl(double const *const x, double *const g) const
{
    g[0] = x[0] + 2.0 * x[1] - x[2] - 2.0;
    g[1] = x[0] * x[0] + x[1] * x[1] + x[2] * x[2] - 10.0;
}

void Jo2::EvalDF_impl(double const *const x, double *const df, size_t const objective_index) const
{
    switch (objective_index)
    {
    case 0:
        df[0] = 2.0 * x[0] - 2.0;
        df[1] = 2.0 * x[1] - 4.0;
        df[2] = 2.0 * x[2] - 6.0;
        break;
    case 1:
        df[0] = 3.0 / 2.0 * std::pow(x[0] * x[0] + 1.0, -0.5) * 2.0 * x[0];
        df[1] = std::pow(x[1] * x[1] + 1.0, -0.5) / 6.0 * 2.0 * x[1];
        df[2] = 0.02 * (x[2] - 1.0);
        break;
    }
}

void Jo2::EvalDG_impl(double const *const x, double *const dg) const
{
    dg[0] = 1.0;
    dg[1] = 2.0 * x[0];
    dg[2] = 2.0;
    dg[3] = 2.0 * x[1];
    dg[4] = -1.0;
    dg[5] = 2.0 * x[2];
}

void Jo2::EvalD2F_impl(double const *const x, double *const d2f, size_t const objective_index) const
{
    switch (objective_index)
    {
    case 0:
        d2f[0] = 2.0;
        d2f[1] = 2.0;
        d2f[2] = 2.0;
        break;
    case 1:
        d2f[0] = - 0.75 * std::pow(x[0] * x[0] + 1.0, -3.0 / 2.0) * 4.0 * x[0] * x[0]
                 + 3.0 / 2.0 * std::pow(x[0] * x[0] + 1.0, -0.5) * 2.0;
        d2f[1] = - std::pow(x[1] * x[1] + 1.0, -3.0 / 2.0) / 12.0 * 4.0 * x[1] * x[1]
                 + std::pow(x[1] * x[1] + 1.0, -0.5) / 6.0 * 2.0;
        d2f[2] = 0.02;
        break;
    }
}

void Jo2::EvalD2G_impl(double const *x, double const *mu, double *d2g, size_t const objective_index) const
{
    d2g[0] = 2.0 * mu[2];
    d2g[1] = 2.0 * mu[2];
    d2g[2] = 2.0 * mu[2];
}

} // namespace test_problems
