#pragma once

#include "../problem_formulation/MONLP.hpp"
#include <limits>


namespace test_problems
{

double constexpr POS_INF = std::numeric_limits<double>::infinity();
double constexpr NEG_INF = -POS_INF;


class ABC_comp : public mosqp::MONLP
{
public:
    ABC_comp();
private:
    double EvalF_impl(double const *x, size_t objective_index) const override;
    void EvalG_impl(double const *x, double *g) const override;
    void EvalDF_impl(double const *x, double *df, size_t objective_index) const override;
    void EvalDG_impl(double const *x, double *dg) const override;
    void EvalD2F_impl(double const *x, double *d2f, size_t objective_index) const override;
    void EvalD2G_impl(double const *x, double const *mu, double *d2g, size_t objective_index) const override;
};


class BK1 : public mosqp::MONLP
{
public:
    BK1();
private:
    double EvalF_impl(double const *x, size_t objective_index) const override;
    void EvalG_impl(double const *x, double *g) const override;
    void EvalDF_impl(double const *x, double *df, size_t objective_index) const override;
    void EvalDG_impl(double const *x, double *dg) const override;
    void EvalD2F_impl(double const *x, double *d2f, size_t objective_index) const override;
    void EvalD2G_impl(double const *x, double const *mu, double *d2g, size_t objective_index) const override;
};


class BNH : public mosqp::MONLP
{
public:
    BNH();
private:
    double EvalF_impl(double const *x, size_t objective_index) const override;
    void EvalG_impl(double const *x, double *g) const override;
    void EvalDF_impl(double const *x, double *df, size_t objective_index) const override;
    void EvalDG_impl(double const *x, double *dg) const override;
    void EvalD2F_impl(double const *x, double *d2f, size_t objective_index) const override;
    void EvalD2G_impl(double const *x, double const *mu, double *d2g, size_t objective_index) const override;
};


class GE3 : public mosqp::MONLP
{
public:
    GE3();
private:
    double EvalF_impl(double const *x, size_t objective_index) const override;
    void EvalG_impl(double const *x, double *g) const override;
    void EvalDF_impl(double const *x, double *df, size_t objective_index) const override;
};


class OSY : public mosqp::MONLP
{
public:
    OSY();
private:
    double EvalF_impl(double const *x, size_t objective_index) const override;
    void EvalG_impl(double const *x, double *g) const override;
    void EvalDF_impl(double const *x, double *df, size_t objective_index) const override;
    void EvalDG_impl(double const *x, double *dg) const override;
    void EvalD2F_impl(double const *x, double *d2f, size_t objective_index) const override;
    void EvalD2G_impl(double const *x, double const *mu, double *d2g, size_t objective_index) const override;
};

} // namespace test_problems
