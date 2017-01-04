#include "TestProblems.hpp"
#include "../problem_formulation/MatrixStructure.hpp"
#include "../problem_formulation/MONLP.hpp"
#include <cstddef>


namespace test_problems
{

using mosqp::MatrixStructure;

/*
# x variable
var x{1..2} >=-5, <=10;

minimize f1:
x[1]^2+x[2]^2;
minimize f2:
(x[1]-5)^2+(x[2]-5)^2;
*/

BK1::BK1() : MONLP(
    2,  // number of variables
    0,  // number of constraints
    2,  // number of objectives

    { MatrixStructure({ 1, 2 }),             // structure of df1
      MatrixStructure({ 1, 2 }) },           // structure of df2
    MatrixStructure(false),                  // structure of dg
    { MatrixStructure({ 1, 2 },{ 1, 2 }),    // structure of hm1
      MatrixStructure({ 1, 2 },{ 1, 2 }) },  // structure of hm2
    true, true, true,          // got all user derivatives

    { -5.0, -5.0 },            // lower bounds on x
    { 10.0, 10.0 },            // upper bounds on x
    { },         // lower bounds on g
    { })  // upper bounds on g
{
}

double BK1::EvalF_impl(double const *const x, size_t const objective_index) const
{
    switch (objective_index)
    {
    case 0:
        return x[0] * x[0] + x[1] * x[1];
    case 1:
        return (x[0] - 5) * (x[0] - 5) + (x[1] - 5) * (x[1] - 5);
    }
    return 0;
}

void BK1::EvalG_impl(double const *const x, double *const g) const
{
}

void BK1::EvalDF_impl(double const *const x, double *const df, size_t const objective_index) const
{
    switch (objective_index)
    {
    case 0:
        df[0] = 2 * x[0];
        df[1] = 2 * x[1];
        break;
    case 1:
        df[0] = 2 * (x[0] - 5);
        df[1] = 2 * (x[1] - 5);
        break;
    }
}

void BK1::EvalDG_impl(double const *const x, double *const dg) const
{
}

void BK1::EvalD2F_impl(double const *const x, double *const d2f, size_t const objective_index) const
{
    switch (objective_index)
    {
    case 0:
        d2f[0] = 2;
        d2f[1] = 2;
        break;
    case 1:
        d2f[0] = 2;
        d2f[1] = 2;
        break;
    }
}

void BK1::EvalD2G_impl(double const *x, double const *mu, double *d2g, size_t const /*objective_index*/) const
{
    d2g[0] = 0;
    d2g[1] = 0;
}

} // namespace test_problems
