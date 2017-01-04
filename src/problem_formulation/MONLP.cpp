#include "MONLP.hpp"
#include "ConstrainedProblem.hpp"
#include "MatrixStructure.hpp"
#include <cstddef>
#include <vector>


namespace mosqp
{

MONLP::MONLP(size_t num_variables, size_t num_constraints, size_t num_objectives,
             std::vector<MatrixStructure> structure_df,
             MatrixStructure structure_dg,
             std::vector<MatrixStructure> structure_hm,
             bool user_df, bool user_dg, bool user_hm,
    std::vector<double> x_l, std::vector<double> x_u,
    std::vector<double> g_l, std::vector<double> g_u)
    : ConstrainedProblem(num_variables, num_constraints, structure_dg, user_dg, x_l, x_u, g_l, g_u),
      numObjectives(num_objectives), structureDF(structure_df), structureHM(structure_hm),
      userDF(user_df), userHM(user_hm)
{
}

size_t MONLP::GetNumObjectives() const
{
    return numObjectives;
}

MatrixStructure const & MONLP::GetStructureDF(size_t const objective_index) const
{
    return structureDF[objective_index];
}

MatrixStructure const & MONLP::GetStructureHM(size_t const objective_index) const
{
    return structureHM[objective_index];
}

bool MONLP::UserDF() const
{
    return userDF;
}

bool MONLP::UserHM() const
{
    return userHM;
}

void MONLP::EvalF(double const *const x, double *const f) const
{
    for (size_t i = 0; i < numObjectives; i += 1)
    {
        numEvalF += 1;
        f[i] = EvalF_impl(x, i);
    }
}

double MONLP::EvalF(double const *const x, size_t const objective_index) const
{
    numEvalF += 1;
    return EvalF_impl(x, objective_index);
}

void MONLP::EvalDF(double const *const x, double *const df, size_t const objective_index) const
{
    numEvalDF += 1;
    EvalDF_impl(x, df, objective_index);
}

void MONLP::EvalD2F(double const *const x, double *const d2f, size_t const objective_index) const
{
    numEvalD2F += 1;
    EvalD2F_impl(x, d2f, objective_index);
}

void MONLP::EvalD2G(double const *const x, double const *mu, double *const d2g, size_t const objective_index) const
{
    numEvalD2G += 1;
    EvalD2G_impl(x, mu, d2g, objective_index);
}

size_t MONLP::GetNumEvalF() const
{
    return numEvalF;
}

size_t MONLP::GetNumEvalDF() const
{
    return numEvalDF;
}

size_t MONLP::GetNumEvalD2F() const
{
    return numEvalD2F;
}

size_t MONLP::GetNumEvalD2G() const
{
    return numEvalD2G;
}

void MONLP::EvalDF_impl(double const *const x, double *const df, size_t const objective_index) const
{
}

void MONLP::EvalD2F_impl(double const *const x, double *const d2f, size_t const objective_index) const
{
}

void MONLP::EvalD2G_impl(double const *x, double const *mu, double *d2g, size_t const objective_index) const
{
}

} // namespace mosqp
