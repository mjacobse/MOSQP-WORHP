#include "TestProblems.hpp"
#include "../problem_formulation/MatrixStructure.hpp"
#include "../problem_formulation/MONLP.hpp"
#include <cstddef>


namespace test_problems
{

using mosqp::MatrixStructure;

/*
... variables
var x1 >= 0, <= 10;		# product 1 in tons / day
var x2 >= 0, <= 10;		# product 2 in tons / day

... objective functions
minimize  f1 : (x1 - 4) ^ 2 + x2 ^ 2 + 5 * x1 + 4 * x2 - x1*x2;	 # ... total cost
minimize  f2 : -x1 - 2 * x2;				                     # ... production

... constraints
subject to

... achieve profit min & limit labor cost
profit : x1*x2 >= 3;
labor:  5 * x1 + 4 * x2 <= 25;

... daily production goal
product : x1 + 2 * x2 >= 5;
*/

ABC_comp::ABC_comp() : MONLP(
    2,  // number of variables
    3,  // number of constraints
    2,  // number of objectives

    { MatrixStructure({ 1, 2 }),                                  // structure of df1
      MatrixStructure({ 1, 2 }) },                                // structure of df2
    MatrixStructure({ 1, 2, 3, 1, 2, 3 }, { 1, 1, 1, 2, 2, 2 }),  // structure of dg
    { MatrixStructure({ 2, 1, 2 }, { 1, 1, 2 }),                  // structure of hm1
      MatrixStructure({ 1, 2 }, { 1, 2 }) },                      // structure of hm2
    true, true, true,          // got all user derivatives

    { 0, 0 },                  // lower bounds on x
    { 10, 10 },                // upper bounds on x
    { 3, NEG_INF, 5 },         // lower bounds on g
    { POS_INF, 25, POS_INF })  // upper bounds on g
{
}

double ABC_comp::EvalF_impl(double const *const x, size_t const objective_index) const
{
    switch (objective_index)
    {
    case 0:
        return ((x[0] - 4) * (x[0] - 4) + x[1] * x[1] + 5 * x[0] + 4 * x[1] - x[0] * x[1]);
    case 1:
        return (-x[0] - 2 * x[1]);
    }
    return 0;
}

void ABC_comp::EvalG_impl(double const *const x, double *const g) const
{
    g[0] = x[0] * x[1];
    g[1] = 5 * x[0] + 4 * x[1];
    g[2] = x[0] + 2 * x[1];
}

void ABC_comp::EvalDF_impl(double const *const x, double *const df, size_t const objective_index) const
{
    switch (objective_index)
    {
    case 0:
        df[0] = (2 * (x[0] - 4) + 5 - x[1]);
        df[1] = (2 * x[1] + 4 - x[0]);
        break;
    case 1:
        df[0] = -1;
        df[1] = -2;
        break;
    }
}

void ABC_comp::EvalDG_impl(double const *const x, double *const dg) const
{
    dg[0] = x[1];
    dg[1] = 5;
    dg[2] = 1;
    dg[3] = x[0];
    dg[4] = 4;
    dg[5] = 2;
}

void ABC_comp::EvalD2F_impl(double const *const x, double *const d2f, size_t const objective_index) const
{
    switch (objective_index)
    {
    case 0:
        d2f[0] = -1;
        d2f[1] = 2;
        d2f[2] = 2;
        break;
    case 1:
        d2f[0] = 0;
        d2f[1] = 0;
        break;
    }
}

void ABC_comp::EvalD2G_impl(double const *x, double const *mu, double *d2g, size_t const objective_index) const
{
    switch (objective_index)
    {
    case 0:
        d2g[0] = 0;
        d2g[1] = 0;
        d2g[2] = 0;
        break;
    case 1:
        d2g[0] = 0;
        d2g[1] = 0;
        break;
    }
}

} // namespace test_problems
