#pragma once

#include "ConstrainedProblem.hpp"
#include "MatrixStructure.hpp"
#include "../algorithm/Point.hpp"
#include <cstddef>
#include <vector>


namespace mosqp
{

class NLP : public ConstrainedProblem
{
public:
    MatrixStructure const & GetStructureDF() const;
    MatrixStructure const & GetStructureHM() const;
    bool UserDF() const;
    bool UserHM() const;

    double EvalF(double const *x) const;
    virtual double GetF(Point const &point) const = 0;
    virtual void GetG(Point const &point, double *g) const = 0;
    void EvalDF(double const *x, double *df) const;
    void EvalHM(double const *x, double const *mu, double scale_obj, double *hm) const;

    virtual void SetParameters(std::vector<double> parameters);

protected:
    NLP(size_t num_variables, size_t num_constraints,
        MatrixStructure structure_df, MatrixStructure structure_dg, MatrixStructure structure_hm,
        bool user_df, bool user_dg, bool user_hm,
        std::vector<double> x_l, std::vector<double> x_u,
        std::vector<double> g_l, std::vector<double> g_u);
    MatrixStructure const structureDF;
    MatrixStructure const structureHM;
    bool const userDF;
    bool const userHM;

    std::vector<double> parameters;

    // bookkeeping variables
    mutable size_t numEvalF = 0;
    mutable size_t numEvalDF = 0;
    mutable size_t numEvalHM = 0;

    virtual double EvalF_impl(double const *x) const = 0;
    virtual void EvalDF_impl(double const *x, double *df) const;
    virtual void EvalHM_impl(double const *x, double const *mu, double scale_obj, double *hm) const;
};

} // namespace mosqp
