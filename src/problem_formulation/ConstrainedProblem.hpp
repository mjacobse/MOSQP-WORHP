#pragma once

#include "MatrixStructure.hpp"
#include <cstddef>
#include <vector>


namespace mosqp
{

// TODO: remove pointers, use iterators or vectors or something
class ConstrainedProblem
{
public:
    size_t GetNumVariables() const;
    size_t GetNumConstraints() const;

    virtual MatrixStructure const & GetStructureDG() const;
    bool UserDG() const;

    std::vector<double> const & GetXL() const;
    std::vector<double> const & GetXU() const;
    virtual std::vector<double> const & GetGL() const;
    virtual std::vector<double> const & GetGU() const;

    void EvalG(double const *x, double *g) const;
    void EvalDG(double const *x, double *dg) const;

    size_t GetNumEvalG() const;
    size_t GetNumEvalDG() const;

protected:
    ConstrainedProblem(size_t num_variables, size_t num_constraints,
                       MatrixStructure structure_dg, bool user_dg,
                       std::vector<double> x_l, std::vector<double> x_u,
                       std::vector<double> g_l, std::vector<double> g_u);
    size_t const numVariables;
    size_t const numConstraints;
    MatrixStructure const structureDG;
    bool const userDG;
    std::vector<double> const xl;
    std::vector<double> const xu;
    std::vector<double> const gl;
    std::vector<double> const gu;

    mutable size_t numEvalG = 0;
    mutable size_t numEvalDG = 0;

    virtual void EvalG_impl(double const *x, double *g) const = 0;
    virtual void EvalDG_impl(double const *x, double *dg) const;
};

} // namespace mosqp
