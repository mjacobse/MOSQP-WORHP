#pragma once

#include "../problem_formulation/MONLP.hpp"
#include <cstddef>
#include <vector>


namespace mosqp
{

class Point
{
public:
    static double TOL_FEAS;
    static double TOL_DOMINATION;

    Point(std::vector<double> const &x, MONLP const &monlp);
    Point(std::vector<double> const &x, std::vector<double> const &lambda, std::vector<double> const &mu,
          MONLP const &monlp);
    // Creates a random point within the given bounds.
    Point(std::vector<double> const &lower_bounds, std::vector<double> const &upper_bounds, MONLP const &monlp);

    // Marks this point as stopped.
    void SetStopped(bool stopped) const;

    double GetObjectiveValue(size_t objective_index) const;
    std::vector<double> const & GetObjectiveValues() const;
    std::vector<double> const & GetX() const;
    std::vector<double> const & GetLambda() const;
    std::vector<double> const & GetMu() const;
    std::vector<double> const & GetConstraints() const;
    double GetDistance(double *other_x) const;
    bool HasMultipliers() const;
    bool IsStopped() const;
    bool IsFeasible() const;
    bool IsDominated(Point const &point) const;
    bool IsSmaller(Point const &point, size_t objective_index) const;

private:
    std::vector<double> x;
    std::vector<double> lambda;
    std::vector<double> mu;
    std::vector<double> f;
    std::vector<double> g;
    std::vector<double> cv;
    mutable bool stopped;

    void UpdateFunctionValues(MONLP const &monlp);
};

} // namespace mosqp
