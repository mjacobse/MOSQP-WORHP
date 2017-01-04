#pragma once

#include "NLP.hpp"
#include "MONLP.hpp"
#include <cstddef>


namespace mosqp
{

class SingleMONLP : public NLP
{
public:
    SingleMONLP(MONLP const &monlp, size_t objective_index);

    double GetF(Point const &point) const override;
    void GetG(Point const &point, double *g) const override;

protected:
    size_t const objectiveIndex;
    MONLP const &monlp;

    double EvalF_impl(double const *x) const override;
    void EvalG_impl(double const *x, double *g) const override;
    void EvalDF_impl(double const *x, double *df) const override;
    void EvalDG_impl(double const *x, double *dg) const override;
    void EvalHM_impl(double const *x, double const *mu, double scale_obj, double *hm) const override;
};

} // namespace mosqp
