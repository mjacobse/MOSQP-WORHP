#include "Point.hpp"
#include "../problem_formulation/MONLP.hpp"
#include <algorithm>
#include <cstddef>
#include <random>
#include <vector>


namespace mosqp
{

double Point::TOL_FEAS;
double Point::TOL_DOMINATION;

Point::Point(std::vector<double> const &x, MONLP const &monlp)
    : x(x), lambda(), mu(), stopped(false),
      f(monlp.GetNumObjectives()), g(monlp.GetNumConstraints()),
      cv(monlp.GetNumVariables() * 2 + monlp.GetNumConstraints() * 2)
{
    UpdateFunctionValues(monlp);
}

Point::Point(std::vector<double> const &x, std::vector<double> const &lambda,
             std::vector<double> const &mu, MONLP const &monlp)
    : x(x), lambda(lambda), mu(mu), stopped(false),
      f(monlp.GetNumObjectives()), g(monlp.GetNumConstraints()),
      cv(monlp.GetNumVariables() * 2 + monlp.GetNumConstraints() * 2)
{
    UpdateFunctionValues(monlp);
}

static std::default_random_engine re;
Point::Point(std::vector<double> const &lower_bounds, std::vector<double> const &upper_bounds,
             MONLP const &monlp)
    : x(monlp.GetNumVariables()), lambda(), mu(),
      f(monlp.GetNumObjectives()), g(monlp.GetNumConstraints()),
      cv(monlp.GetNumVariables() * 2 + monlp.GetNumConstraints() * 2)
{
    std::vector<double> new_x;
    std::uniform_real_distribution<double> unif;
    for (size_t i = 0; i < lower_bounds.size(); i += 1)
    {
        unif = std::uniform_real_distribution<double>(lower_bounds[i], upper_bounds[i]);
        new_x.emplace_back(unif(re));
    }

    x = new_x;
    UpdateFunctionValues(monlp);
}

void Point::SetStopped(bool const stopped) const
{
    this->stopped = stopped;
}

void Point::UpdateFunctionValues(MONLP const &monlp)
{
    monlp.EvalF(x.data(), f.data());

    size_t const num_constraints = monlp.GetNumConstraints();
    size_t const num_variables = monlp.GetNumVariables();
    g = std::vector<double>(num_constraints);
    std::vector<double> const &g_l = monlp.GetGL();
    std::vector<double> const &g_u = monlp.GetGU();
    std::vector<double> const &x_l = monlp.GetXL();
    std::vector<double> const &x_u = monlp.GetXU();
    monlp.EvalG(x.data(), g.data());

    size_t j = 0;
    for (size_t i = 0; i < num_variables; i += 1)
    {
        cv[j++] = std::max(x_l[i] - x[i], 0.0);
        cv[j++] = std::max(x[i] - x_u[i], 0.0);
    }

    for (size_t i = 0; i < num_constraints; i += 1)
    {
        cv[j++] = std::max(g_l[i] - g[i], 0.0);
        cv[j++] = std::max(g[i] - g_u[i], 0.0);
    }
}

double Point::GetObjectiveValue(size_t const objective_index) const
{
    return f[objective_index];
}

std::vector<double> const & Point::GetObjectiveValues() const
{
    return f;
}

std::vector<double> const & Point::GetX() const
{
    return x;
}

std::vector<double> const & Point::GetLambda() const
{
    return lambda;
}

std::vector<double> const & Point::GetMu() const
{
    return mu;
}

std::vector<double> const & Point::GetConstraints() const
{
    return g;
}

double Point::GetDistance(double *const other_x) const
{
    double distance = 0;
    for (int i = 0; i < x.size(); i += 1)
    {
        distance += (x[i] - other_x[i]) * (x[i] - other_x[i]);
    }

    distance = std::sqrt(distance);
    return distance;
}

bool Point::HasMultipliers() const
{
    return lambda.size() != 0;
}

bool Point::IsStopped() const
{
    return stopped;
}

bool Point::IsFeasible() const
{
    for (double violation : cv)
    {
        if (violation > TOL_FEAS)
        {
            return false;
        }
    }
    return true;
}

bool Point::IsSmaller(Point const &point, size_t const objective_index) const
{
    return f[objective_index] < point.f[objective_index];
}

bool Point::IsDominated(Point const &point) const
{
    for (size_t i = 0; i < f.size(); i += 1)
    {
        if (f[i] < point.f[i])
        {
            return false;
        }
    }

    if (*std::max_element(cv.begin(), cv.end()) < *std::max_element(point.cv.begin(), point.cv.end()))
    {
        return false;
    }

    return true;
}

} // namespace mosqp
