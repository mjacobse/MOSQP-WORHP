#include "SingleMONLP.hpp"
#include "MONLP.hpp"
#include <cstddef>


namespace mosqp
{

SingleMONLP::SingleMONLP(MONLP const &monlp, size_t const objective_index)
    : NLP(monlp.GetNumVariables(), monlp.GetNumConstraints(),
          monlp.GetStructureDF(objective_index),
          monlp.GetStructureDG(),
          monlp.GetStructureHM(objective_index),
          monlp.UserDF(), monlp.UserDG(), monlp.UserHM(),
          monlp.GetXL(), monlp.GetXU(),
          monlp.GetGL(), monlp.GetGU()),
    monlp(monlp), objectiveIndex(objective_index)
{
}

double SingleMONLP::EvalF_impl(double const *const x) const
{
    return monlp.EvalF(x, objectiveIndex);
}

double SingleMONLP::GetF(Point const &point) const
{
    return point.GetObjectiveValue(objectiveIndex);
}

void SingleMONLP::EvalG_impl(double const *const x, double *const g) const
{
    monlp.EvalG(x, g);
}

void SingleMONLP::GetG(Point const &point, double *g) const
{
    std::vector<double> const &original_g = point.GetConstraints();
    std::copy(original_g.cbegin(), original_g.cend(), g);
}

void SingleMONLP::EvalDF_impl(double const *const x, double *const df) const
{
    return monlp.EvalDF(x, df, objectiveIndex);
}

void SingleMONLP::EvalDG_impl(double const *const x, double *const dg) const
{
    return monlp.EvalDG(x, dg);
}

void SingleMONLP::EvalHM_impl(double const *const x, double const *mu, double const scale_obj, double *const hm) const
{
    size_t const hm_nnz = monlp.GetStructureHM(objectiveIndex).GetNumNonZeros();
    std::vector<double> d2g(hm_nnz, 0.0);

    monlp.EvalD2F(x, hm, objectiveIndex);
    monlp.EvalD2G(x, mu, d2g.data(), objectiveIndex);

    for (size_t i = 0; i < hm_nnz; i += 1)
    {
        hm[i] = scale_obj * hm[i] + d2g[i];
    }
}

}
