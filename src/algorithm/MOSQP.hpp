#pragma once

#include "ParetoFront.hpp"
#include "Point.hpp"
#include "../problem_formulation/MONLP.hpp"
#include <fstream>
#include <vector>


namespace mosqp
{

class Parameters
{
public:
    Parameters();
    int maxPoints;
    int numCompletionTries;
    double TOL_FEAS;
    double TOL_DOMINATION;
    int SPREAD_MAX_STEPS;
    double SPREAD_ARMIJO_MIN_ALPHA;
    double SPREAD_ARMIJO_BETA;
    double SPREAD_MIN_SEARCH_LENGTH;
    int REFINE_MAX_STEPS;
    double REFINE_ARMIJO_MIN_ALPHA;
    double REFINE_ARMIJO_BETA;
    double REFINE_MIN_SEARCH_LENGTH;
};


class MOSQP
{
public:
    MOSQP(MONLP &monlp,
          std::vector<Point> initial_points = std::vector<Point>(),
          Parameters parameters = Parameters());
    ~MOSQP();

    ParetoFront Solve();

private:
    // The instance of the (derived) multiobjective problem.
    MONLP const &monlp;
    // Initial points provided by the user (empty if none).
    std::vector<Point> initialPoints;
    // Parameters currently used by MOSQP.
    Parameters const parameters;
    // Current iteration of the Pareto front. During solving this will be updated
    // and eventually returned once finished.
    ParetoFront paretoFront;
    // scalings that are determined when finding extreme pareto points
    // and are used in the combined objective function to scale each objective to
    // similar values
    std::vector<double> scalings;

    std::ofstream log;

    // First solving stage as proposed in the MOSQP paper.
    // Takes the given initial points by the user and adds some more, depending
    // on the chosen parameters.
    void CompleteInitialPoints();
    // Second solving stage.
    // Tries to spread the Pareto front as much as possible.
    void SpreadParetoFront();
    // Proposed step after the second solving stage.
    // Adds the extreme points, i.e. the (ideally) global minima of each objective function.
    void AddExtremeParetoPoints();
    // Third solving stage.
    // Drives the spread front to Pareto optimality.
    void RefineParetoFront();
};

} // namespace mosqp
