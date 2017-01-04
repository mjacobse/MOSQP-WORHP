#include "TestProblems.hpp"
#include "../problem_formulation/MatrixStructure.hpp"
#include "../problem_formulation/MONLP.hpp"
#include <cstddef>


namespace test_problems
{

using mosqp::MatrixStructure;

/*
x variable
var x{1..2};

minimize f1:
4*x[1]^2+4*x[2]^2;
minimize f2:
(x[1]-5)^2+(x[2]-5)^2;

subject to c1:
(x[1]-5)^2+x[2]^2<=25;
subject to c2:
(x[1]-8)^2+(x[2]+3)^2>=7.7;

subject to bounds1:
0.0 <= x[1] <= 5.0;
subject to bounds2:
0.0 <= x[2] <= 3.0;
*/

BNH::BNH() : MONLP(
    2,  // number of variables
    2,  // number of constraints
    2,  // number of objectives

    { MatrixStructure({ 1, 2 }),                      // structure of df1
      MatrixStructure({ 1, 2 }) },                    // structure of df2
    MatrixStructure({ 1, 2, 1, 2 }, { 1, 1, 2, 2 }),  // structure of dg
    { MatrixStructure({ 1, 2 }, { 1, 2 }),            // structure of hm1
      MatrixStructure({ 1, 2 }, { 1, 2 }) },          // structure of hm2
    true, true, true,  // got all user derivatives

    { 0, 0 },          // lower bounds on x
    { 5, 3 },          // upper bounds on x
    { NEG_INF, 7.7 },  // lower bounds on g
    { 25, POS_INF })   // upper bounds on g
{
}

double BNH::EvalF_impl(double const *const x, size_t const objective_index) const
{
    switch (objective_index)
    {
    case 0:
        return 4 * x[0] * x[0] + 4 * x[1] * x[1];
    case 1:
        return (x[0] - 5)*(x[0] - 5) + (x[1] - 5)*(x[1] - 5);
    }
    return 0;
}

void BNH::EvalG_impl(double const *const x, double *const g) const
{
    g[0] = (x[0] - 5)*(x[0] - 5) + x[1] * x[1];
    g[1] = (x[0] - 8)*(x[0] - 8) + (x[1] + 3)*(x[1] + 3);
}

void BNH::EvalDF_impl(double const *const x, double *const df, size_t const objective_index) const
{
    switch (objective_index)
    {
    case 0:
        df[0] = 8 * x[0];
        df[1] = 8 * x[1];
        break;
    case 1:
        df[0] = 2 * (x[0] - 5);
        df[1] = 2 * (x[1] - 5);
        break;
    }
}

void BNH::EvalDG_impl(double const *const x, double *const dg) const
{
    dg[0] = 2 * (x[0] - 5);
    dg[1] = 2 * (x[0] - 8);
    dg[2] = 2 * x[1];
    dg[3] = 2 * (x[1] + 3);
}

void BNH::EvalD2F_impl(double const *const x, double *const d2f, size_t const objective_index) const
{
    switch (objective_index)
    {
    case 0:
        d2f[0] = 8;
        d2f[1] = 8;
        break;
    case 1:
        d2f[0] = 2;
        d2f[1] = 2;
        break;
    }
}

void BNH::EvalD2G_impl(double const *const x, double const *const mu, double *const d2g, size_t const /*objective_index*/) const
{
    d2g[0] = 2 * mu[0] + 2 * mu[1];
    d2g[1] = 2 * mu[0] + 2 * mu[1];
}

} // namespace test_problems
