#include "NLP.hpp"
#include "ConstrainedProblem.hpp"
#include "MatrixStructure.hpp"
#include <cstddef>
#include <vector>


namespace mosqp
{

NLP::NLP(size_t num_variables, size_t num_constraints,
         MatrixStructure structure_df,
         MatrixStructure structure_dg,
         MatrixStructure structure_hm,
         bool user_df, bool user_dg, bool user_hm,
         std::vector<double> x_l, std::vector<double> x_u,
         std::vector<double> g_l, std::vector<double> g_u)
    : ConstrainedProblem(num_variables, num_constraints, structure_dg, user_dg, x_l, x_u, g_l, g_u),
      structureDF(structure_df), structureHM(structure_hm), userDF(user_df), userHM(user_hm)
{
}

MatrixStructure const & NLP::GetStructureDF() const
{
    return structureDF;
}

MatrixStructure const & NLP::GetStructureHM() const
{
    return structureHM;
}

bool NLP::UserDF() const
{
    return userDF;
}

bool NLP::UserHM() const
{
    return userHM;
}

double NLP::EvalF(double const *const x) const
{
    numEvalF += 1;
    return EvalF_impl(x);
}

void NLP::EvalDF(double const *const x, double *const df) const
{
    numEvalDF += 1;
    EvalDF_impl(x, df);
}

void NLP::EvalHM(double const *const x, double const *mu, double const scale_obj, double *const hm) const
{
    numEvalHM += 1;
    EvalHM_impl(x, mu, scale_obj, hm);
}

void NLP::SetParameters(std::vector<double> parameters)
{
    this->parameters = parameters;
}

void NLP::EvalDF_impl(double const *const x, double *const df) const
{
}

void NLP::EvalHM_impl(double const *const x, double const *mu, double const scale_obj, double *const hm) const
{
}

} // namespace mosqp
