#include "TestProblems.hpp"
#include "../problem_formulation/MatrixStructure.hpp"
#include "../problem_formulation/MONLP.hpp"
#include <cstddef>


namespace test_problems
{

using mosqp::MatrixStructure;

/*
# x variable
var x{1..6};

minimize f1:
-(25*(x[1]-2)^2+(x[2]-2)^2+(x[3]-1)^2+(x[4]-4)^2+(x[5]-1)^2);
minimize f2:
sum {i in 1..6} (x[i]^2);

subject to c1:
x[1]+x[2]-2>=0;
subject to c2:
6-x[1]-x[2]>=0;
subject to c3:
2-x[2]+x[1]>=0;
subject to c4:
2-x[1]+3*x[2]>=0;
subject to c5:
4-(x[3]-3)^2-x[4]>=0;
subject to c6:
(x[5]-3)^2+x[6]-4>=0;

subject to bounds1:
0.0 <= x[1] <= 10;
subject to bounds2:
0.0 <= x[2] <= 10;
subject to bounds3:
1.0 <= x[3] <= 5;
subject to bounds4:
0.0 <= x[4] <= 6;
subject to bounds5:
1.0 <= x[5] <= 5;
subject to bounds6:
0.0 <= x[6] <= 10;
*/
OSY::OSY() : MONLP(
    6,  // number of variables
    6,  // number of constraints
    2,  // number of objectives

    { MatrixStructure({ 1, 2, 3, 4, 5 }),                             // structure of df1
      MatrixStructure({ 1, 2, 3, 4, 5, 6 }) },                        // structure of df2
    MatrixStructure({ 1, 2, 3, 4, 1, 2, 3, 4, 5, 5, 6, 6 },
                    { 1, 1, 1, 1, 2, 2, 2, 2, 3, 4, 5, 6 }),          // structure of dg
    { MatrixStructure({ 1, 2, 3, 4, 5, 6 }, { 1, 2, 3, 4, 5, 6 }),    // structure of hm1
      MatrixStructure({ 1, 2, 3, 4, 5, 6 }, { 1, 2, 3, 4, 5, 6 }) },  // structure of hm2
    true, true, true,          // got all user derivatives

    { 0, 0, 1, 0, 1, 0 },     // lower bounds on x
    { 10, 10, 5, 6, 5, 10 },  // upper bounds on x
    { 0, 0, 0, 0, 0, 0 },     // lower bounds on g
    { POS_INF, POS_INF, POS_INF, POS_INF, POS_INF, POS_INF })  // upper bounds on g
{
}

double OSY::EvalF_impl(double const *const x, size_t const objective_index) const
{
    switch (objective_index)
    {
    case 0:
        return -(25 * (x[0] - 2)*(x[0] - 2) + (x[1] - 2)*(x[1] - 2) +(x[2] - 1)*(x[2] - 1) + (x[3] - 4)*(x[3] - 4) + (x[4] - 1)*(x[4] - 1));
    case 1:
        return x[0] * x[0] + x[1] * x[1] + x[2] * x[2] + x[3] * x[3] + x[4] * x[4] + x[5] * x[5];
    }
    return 0;
}

void OSY::EvalG_impl(double const *const x, double *const g) const
{
    g[0] = x[0] + x[1] - 2;
    g[1] = 6 - x[0] - x[1];
    g[2] = 2 - x[1] + x[0];
    g[3] = 2 - x[0] + 3 * x[1];
    g[4] = 4 - (x[2] - 3)*(x[2] - 3) - x[3];
    g[5] = (x[4] - 3)*(x[4] - 3) + x[5] - 4;
}

void OSY::EvalDF_impl(double const *const x, double *const df, size_t const objective_index) const
{
    switch (objective_index)
    {
    case 0:
        df[0] = -50 * (x[0] - 2);
        df[1] = -2 * (x[1] - 2);
        df[2] = -2 * (x[2] - 1);
        df[3] = -2 * (x[3] - 4);
        df[4] = -2 * (x[4] - 1);
        break;
    case 1:
        df[0] = 2 * x[0];
        df[1] = 2 * x[1];
        df[2] = 2 * x[2];
        df[3] = 2 * x[3];
        df[4] = 2 * x[4];
        df[5] = 2 * x[5];
        break;
    }
}

void OSY::EvalDG_impl(double const *const x, double *const dg) const
{
    dg[0] = 1;
    dg[1] = -1;
    dg[2] = 1;
    dg[3] = -1;
    dg[4] = 1;
    dg[5] = -1;
    dg[6] = -1;
    dg[7] = 3;
    dg[8] = -2 * (x[2] - 3);
    dg[9] = -1;
    dg[10] = 2 * (x[4] - 3);
    dg[11] = 1;
}

void OSY::EvalD2F_impl(double const *const x, double *const d2f, size_t const objective_index) const
{
    switch (objective_index)
    {
    case 0:
        d2f[0] = -50;
        d2f[1] = -2;
        d2f[2] = -2;
        d2f[3] = -2;
        d2f[4] = -2;
        d2f[5] = 0;
        break;
    case 1:
        d2f[0] = 2;
        d2f[1] = 2;
        d2f[2] = 2;
        d2f[3] = 2;
        d2f[4] = 2;
        d2f[5] = 2;
        break;
    }
}

void OSY::EvalD2G_impl(double const *const x, double const *const mu, double *const d2g, size_t const objective_index) const
{
    d2g[0] = 0;
    d2g[1] = 0;
    d2g[2] = -2 * mu[4];
    d2g[3] = 0;
    d2g[4] = 2 * mu[5];
    d2g[5] = 0;
}

} // namespace test_problems
