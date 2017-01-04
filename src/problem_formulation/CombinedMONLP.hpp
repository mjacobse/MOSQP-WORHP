#pragma once

#include "MatrixStructure.hpp"
#include "NLP.hpp"
#include "MONLP.hpp"
#include "../algorithm/Point.hpp"
#include <vector>


namespace mosqp
{

// TODO: Give user DF and HM and implement those functions.
// TODO: Scale different objective functions.
class CombinedMONLP : public NLP
{
public:
    CombinedMONLP(MONLP const &monlp, std::vector<double> const &scalings);

    virtual std::vector<double> const & GetGL() const override;
    virtual std::vector<double> const & GetGU() const override;

    double GetF(Point const &point) const override;
    void GetG(Point const &point, double *g) const override;

    void SetParameters(std::vector<double> parameters) override;

private:
    MONLP const &monlp;
    std::vector<double> combined_gl;
    std::vector<double> combined_gu;
    std::vector<double> scalings;

    static MatrixStructure GetCombinedStructureDF(MONLP const &monlp);
    static MatrixStructure GetCombinedStructureDG(MONLP const &monlp);
    static MatrixStructure GetCombinedStructureHM(MONLP const &monlp);

    double EvalF_impl(double const *x) const override;
    void EvalG_impl(double const *x, double *g) const override;
    void EvalDF_impl(double const *x, double *df) const override;
    //void EvalDG_impl(double const *x, double *dg) const override;
    void EvalHM_impl(double const *x, double const *mu, double scale_obj, double *hm) const override;
};

} // namespace mosqp
