#pragma once

#include "ConstrainedProblem.hpp"
#include "MatrixStructure.hpp"
#include <cstddef>
#include <vector>


namespace mosqp
{

// TODO: support only some of DF, HM being given by user
//       => vector of userDF and userHM bool
class MONLP : public ConstrainedProblem
{
public:
    size_t GetNumObjectives() const;

    MatrixStructure const & GetStructureDF(size_t objective_index) const;
    MatrixStructure const & GetStructureHM(size_t objective_index) const;
    bool UserDF() const;
    bool UserHM() const;

    void EvalF(double const *x, double *f) const;
    double EvalF(double const *x, size_t objective_index) const;
    void EvalDF(double const *x, double *df, size_t objective_index) const;
    void EvalD2F(double const *x, double *d2f, size_t objective_index) const;
    void EvalD2G(double const *x, double const *mu, double *d2g, size_t objective_index) const;

    size_t GetNumEvalF() const;
    size_t GetNumEvalDF() const;
    size_t GetNumEvalD2F() const;
    size_t GetNumEvalD2G() const;

protected:
    MONLP(size_t num_variables, size_t num_constraints, size_t num_objectives,
          std::vector<MatrixStructure> structure_df,
          MatrixStructure structure_dg,
          std::vector<MatrixStructure> structure_hm,
          bool user_df, bool user_dg, bool user_hm,
          std::vector<double> x_l, std::vector<double> x_u,
          std::vector<double> g_l, std::vector<double> g_u);

    size_t const numObjectives;
    std::vector<MatrixStructure> const structureDF;
    std::vector<MatrixStructure> const structureHM;
    bool const userDF;
    bool const userHM;

    mutable size_t numEvalF = 0;
    mutable size_t numEvalDF = 0;
    mutable size_t numEvalD2F = 0;
    mutable size_t numEvalD2G = 0;

    virtual double EvalF_impl(double const *x, size_t objective_index) const = 0;
    virtual void EvalDF_impl(double const *x, double *df, size_t objective_index) const;
    virtual void EvalD2F_impl(double const *x, double *d2f, size_t objective_index) const;
    virtual void EvalD2G_impl(double const *x, double const *mu, double *d2g, size_t objective_index) const;
};

} // namespace mosqp
