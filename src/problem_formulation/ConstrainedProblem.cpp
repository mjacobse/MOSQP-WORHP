#include "ConstrainedProblem.hpp"
#include "MatrixStructure.hpp"
#include <cstddef>
#include <vector>


namespace mosqp
{

ConstrainedProblem::ConstrainedProblem(size_t num_variables, size_t num_constraints,
    MatrixStructure structure_dg, bool user_dg,
    std::vector<double> x_l, std::vector<double> x_u,
    std::vector<double> g_l, std::vector<double> g_u)
    : numVariables(num_variables), numConstraints(num_constraints), structureDG(structure_dg),
      userDG(user_dg), xl(x_l), xu(x_u), gl(g_l), gu(g_u)
{
}

size_t ConstrainedProblem::GetNumVariables() const
{
    return numVariables;
}

size_t ConstrainedProblem::GetNumConstraints() const
{
    return numConstraints;
}

MatrixStructure const & ConstrainedProblem::GetStructureDG() const
{
    return structureDG;
}

bool ConstrainedProblem::UserDG() const
{
    return userDG;
}

std::vector<double> const & ConstrainedProblem::GetXL() const
{
    return xl;
}

std::vector<double> const & ConstrainedProblem::GetXU() const
{
    return xu;
}

std::vector<double> const & ConstrainedProblem::GetGL() const
{
    return gl;
}

std::vector<double> const & ConstrainedProblem::GetGU() const
{
    return gu;
}

void ConstrainedProblem::EvalG(double const *const x, double *const g) const
{
    numEvalG += 1;
    EvalG_impl(x, g);
}

size_t ConstrainedProblem::GetNumEvalG() const
{
    return numEvalG;
}

size_t ConstrainedProblem::GetNumEvalDG() const
{
    return numEvalDG;
}

void ConstrainedProblem::EvalDG(double const *const x, double *const dg) const
{
    numEvalDG += 1;
    EvalDG_impl(x, dg);
}

void ConstrainedProblem::EvalDG_impl(double const *const x, double *const dg) const
{
}

} // namespace mosqp
