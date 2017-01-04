#include "TestProblems.hpp"
#include "../problem_formulation/MatrixStructure.hpp"
#include "../problem_formulation/MONLP.hpp"
#include <cstddef>
#include <cmath>


namespace test_problems
{

using mosqp::MatrixStructure;
constexpr double M_PI = 3.14159265358979323846;

/*
# Parameters
param pi := 4*atan(1);

# x variable
var x{1..2} >=0.01, <=pi;


minimize f1:
x[1];
minimize f2:
x[2];

subject to c1:
x[1]^2+x[2]^2-1-0.1*cos(16*atan(x[1]/x[2]))>=0;
subject to c2:
(x[1]-0.5)^2+(x[2]-0.5)^2<=0.5;
*/

GE3::GE3() : MONLP(
    2,  // number of variables
    2,  // number of constraints
    2,  // number of objectives

    { MatrixStructure(std::vector<size_t>({ 1 })),    // structure of df1
      MatrixStructure(std::vector<size_t>({ 2 })) },  // structure of df2
    MatrixStructure({ 1, 2, 1, 2 }, { 1, 1, 2, 2 }),  // structure of dg
    { MatrixStructure({ 2, 1, 2 }, { 1, 1, 2 }),      // structure of hm1
      MatrixStructure({ 2, 1, 2 }, { 1, 1, 2 }) },    // structure of hm2
    true, false, false,  // got all user derivatives

    { 0.01, 0.01 },      // lower bounds on x
    { M_PI, M_PI },      // upper bounds on x
    { 0, NEG_INF },      // lower bounds on g
    { POS_INF, 0.5 })    // upper bounds on g
{
}

double GE3::EvalF_impl(double const *const x, size_t const objective_index) const
{
    switch (objective_index)
    {
    case 0:
        return x[0];
    case 1:
        return x[1];
    }
    return 0;
}

void GE3::EvalG_impl(double const *const x, double *const g) const
{
    g[0] = x[0] * x[0] + x[1] * x[1] - 1 - 0.1*std::cos(16 * std::atan2(x[0], x[1]));
    g[1] = (x[0] - 0.5)*(x[0] - 0.5) + (x[1] - 0.5)*(x[1] - 0.5);
}

void GE3::EvalDF_impl(double const *const x, double *const df, size_t const objective_index) const
{
    switch (objective_index)
    {
    case 0:
        df[0] = 1;
        break;
    case 1:
        df[0] = 1;
        break;
    }
}

} // namespace test_problems
