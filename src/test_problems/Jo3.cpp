#include "TestProblems.hpp"
#include "../problem_formulation/MatrixStructure.hpp"
#include "../problem_formulation/MONLP.hpp"
#include <cstddef>


namespace test_problems
{

using mosqp::MatrixStructure;

Jo3::Jo3() : MONLP(
    4,  // number of variables
    0,  // number of constraints
    2,  // number of objectives

    { MatrixStructure({ 1, 2, 3, 4 }),                              // structure of df1
      MatrixStructure({ 1, 2, 3, 4 }) },                            // structure of df2
    MatrixStructure(),           // structure of dg
    { MatrixStructure(false),    // structure of hm1
      MatrixStructure(false) },  // structure of hm2
    true, false, false,

    { -4.0, -4.0, -4.0, -4.0 },  // lower bounds on x
    {  4.0,  4.0,  4.0,  4.0 },  // upper bounds on x
    {  },                        // lower bounds on g
    {  })                        // upper bounds on g
{
}

double Jo3::EvalF_impl(double const *const x, size_t const objective_index) const
{
    switch (objective_index)
    {
    case 0:
        return 1.0 - std::exp(- (x[0] + 0.5) * (x[0] + 0.5) - (x[1] + 0.5) * (x[1] + 0.5) - (x[2] + 0.5) * (x[2] + 0.5) - (x[3] + 0.5) * (x[3] + 0.5));
    case 1:
        return 1.0 - std::exp(- (x[0] - 0.5) * (x[0] - 0.5) - (x[1] - 0.5) * (x[1] - 0.5) - (x[2] - 0.5) * (x[2] - 0.5) - (x[3] - 0.5) * (x[3] - 0.5));
    }
    return 0;
}

void Jo3::EvalG_impl(double const *const x, double *const g) const
{
}

void Jo3::EvalDF_impl(double const *const x, double *const df, size_t const objective_index) const
{
    switch (objective_index)
    {
    case 0:
        df[0] = std::exp(-(x[0] + 0.5) * (x[0] + 0.5) - (x[1] + 0.5) * (x[1] + 0.5) - (x[2] + 0.5) * (x[2] + 0.5) - (x[3] + 0.5) * (x[3] + 0.5)) * 2.0 * (x[0] + 0.5);
        df[1] = std::exp(-(x[0] + 0.5) * (x[0] + 0.5) - (x[1] + 0.5) * (x[1] + 0.5) - (x[2] + 0.5) * (x[2] + 0.5) - (x[3] + 0.5) * (x[3] + 0.5)) * 2.0 * (x[1] + 0.5);
        df[2] = std::exp(-(x[0] + 0.5) * (x[0] + 0.5) - (x[1] + 0.5) * (x[1] + 0.5) - (x[2] + 0.5) * (x[2] + 0.5) - (x[3] + 0.5) * (x[3] + 0.5)) * 2.0 * (x[2] + 0.5);
        df[3] = std::exp(-(x[0] + 0.5) * (x[0] + 0.5) - (x[1] + 0.5) * (x[1] + 0.5) - (x[2] + 0.5) * (x[2] + 0.5) - (x[3] + 0.5) * (x[3] + 0.5)) * 2.0 * (x[3] + 0.5);
        break;
    case 1:
        df[0] = std::exp(-(x[0] - 0.5) * (x[0] - 0.5) - (x[1] - 0.5) * (x[1] - 0.5) - (x[2] - 0.5) * (x[2] - 0.5) - (x[3] - 0.5) * (x[3] - 0.5)) * 2.0 * (x[0] - 0.5);
        df[1] = std::exp(-(x[0] - 0.5) * (x[0] - 0.5) - (x[1] - 0.5) * (x[1] - 0.5) - (x[2] - 0.5) * (x[2] - 0.5) - (x[3] - 0.5) * (x[3] - 0.5)) * 2.0 * (x[1] - 0.5);
        df[2] = std::exp(-(x[0] - 0.5) * (x[0] - 0.5) - (x[1] - 0.5) * (x[1] - 0.5) - (x[2] - 0.5) * (x[2] - 0.5) - (x[3] - 0.5) * (x[3] - 0.5)) * 2.0 * (x[2] - 0.5);
        df[3] = std::exp(-(x[0] - 0.5) * (x[0] - 0.5) - (x[1] - 0.5) * (x[1] - 0.5) - (x[2] - 0.5) * (x[2] - 0.5) - (x[3] - 0.5) * (x[3] - 0.5)) * 2.0 * (x[3] - 0.5);
        break;
    }
}

void Jo3::EvalDG_impl(double const *const x, double *const dg) const
{
}

void Jo3::EvalD2F_impl(double const *const x, double *const d2f, size_t const objective_index) const
{
}

void Jo3::EvalD2G_impl(double const *x, double const *mu, double *d2g, size_t const objective_index) const
{
}

} // namespace test_problems
