#include "MOSQP.hpp"
#include "ParetoFront.hpp"
#include "Point.hpp"
#include "../../include/worhp/worhp.h"
#include "../nlp_solver/WorhpSolver.hpp"
#include "../problem_formulation/MONLP.hpp"
#include "../problem_formulation/CombinedMONLP.hpp"
#include "../problem_formulation/SingleMONLP.hpp"
#include <algorithm>
#include <cassert>
#include <cstddef>
#include <iostream>
#include <memory>


namespace mosqp
{

// set default parameters
Parameters::Parameters()
    : maxPoints(100),
      numCompletionTries(200),
      TOL_FEAS(1e-3),
      TOL_DOMINATION(1e-5),
      SPREAD_MAX_STEPS(20),
      SPREAD_ARMIJO_MIN_ALPHA(1e-5),
      SPREAD_ARMIJO_BETA(0.5),
      SPREAD_MIN_SEARCH_LENGTH(1e-5),
      REFINE_MAX_STEPS(200),
      REFINE_ARMIJO_MIN_ALPHA(1e-5),
      REFINE_ARMIJO_BETA(0.5),
      REFINE_MIN_SEARCH_LENGTH(5e-5)
{
}

MOSQP::MOSQP(MONLP &monlp, std::vector<Point> initial_points, Parameters parameters)
    : monlp(monlp), initialPoints(initial_points), parameters(parameters),
      paretoFront(parameters.maxPoints, monlp.GetNumObjectives(), initial_points),
      log("log.txt")
{
    Point::TOL_FEAS = parameters.TOL_FEAS;
    Point::TOL_DOMINATION = parameters.TOL_DOMINATION;

    log << "Stage 0: Initialisation - " << monlp.GetName() << std::endl;
    CompleteInitialPoints();
    std::cout << "================= Complete Initial Points =================" << std::endl;
    paretoFront.WriteF(log);
}

MOSQP::~MOSQP()
{
    log.close();
}

ParetoFront MOSQP::Solve()
{
    log << "Stage 1: Spread - " << monlp.GetName() << std::endl;
    std::cout << "==================== SpreadParetoFront ====================" << std::endl;
    SpreadParetoFront();

    /*log << "Stage 1.5: Extreme Points - " << monlp.GetName() << std::endl;
    std::cout << "================= AddExtremeParetoPoints ==================" << std::endl;
    AddExtremeParetoPoints();*/
    scalings = { 1.0, 1.0 };

    log << "Stage 2: Refine - " << monlp.GetName() << std::endl;
    std::cout << "==================== RefineParetoFront ====================" << std::endl;
    RefineParetoFront();

    assert(paretoFront.AllFeasible());
    assert(paretoFront.AllNonDominated());

    return paretoFront;
}

void MOSQP::CompleteInitialPoints()
{
    std::vector<double> const &x_l = monlp.GetXL();
    std::vector<double> const &x_u = monlp.GetXU();

    int tries;
    for (tries = 0; tries < parameters.numCompletionTries; tries += 1)
    {
        paretoFront.AddPoint(Point(x_l, x_u, monlp));
        if (paretoFront.IsFull() && paretoFront.AllFeasible() && paretoFront.AllNonDominated())
        {
            break;
        }
    }

    std::cout << "CompleteInitialPoints: Paretofront has " << paretoFront.NumPoints() << " points"
              << " after " << tries << " tries." << std::endl;
}

void MOSQP::SpreadParetoFront()
{
    paretoFront.UnstopAll();

    size_t const num_objectives = monlp.GetNumObjectives();
    // 1 solver for each objective
    std::vector<WorhpSolver> worhp;
    worhp.reserve(num_objectives);
    // the solvers store a reference to the problem so we gotta make sure they
    // actually survive throughout so make a vector of em (idk man, this feels so bad,
    // maybe using unique_ptr or something would be better here)
    std::vector<SingleMONLP> problems;
    problems.reserve(num_objectives);
    for (size_t objective_index = 0; objective_index < num_objectives; objective_index += 1)
    {
        problems.emplace_back(monlp, objective_index);
        worhp.emplace_back(problems.back());
        worhp[objective_index].par.ArmijoMinAlpha = parameters.SPREAD_ARMIJO_MIN_ALPHA;
        worhp[objective_index].par.ArmijoBeta = parameters.SPREAD_ARMIJO_BETA;
        worhp[objective_index].par.ArmijoBetaAres = parameters.SPREAD_ARMIJO_BETA;
        worhp[objective_index].par.KeepAcceptableSol = false;
        worhp[objective_index].par.LowPassFilter = false;
        worhp[objective_index].par.MaxIter = std::numeric_limits<int>::max();
        worhp[objective_index].par.NLPmethod = 1;  // use merit function instead of filter
        worhp[objective_index].par.TolFeas = 1e-20;
        worhp[objective_index].par.TolOpti = 1e-20;
    }

    std::vector<Point> new_points;
    std::vector<double> x;
    std::vector<double> lambda;
    std::vector<double> mu;
    std::vector<double> penalties;
    double step_length;

    for (int step = 0; step < parameters.SPREAD_MAX_STEPS; step += 1)
    {
        new_points.clear();
        for (auto it_point = paretoFront.begin(); it_point != paretoFront.end(); )
        {
            if (!it_point->IsStopped())
            {
                for (size_t i = 0; i < num_objectives; i += 1)
                {
                    worhp[i].SetInitialGuess(*it_point);
                    worhp[i].DoMajorIter();

                    if (worhp[i].cnt.status <= TerminateError)
                    {
                        // couldn't find step in any way so not gonna add this
                        std::cout << "SpreadParetoFront: WORHP terminated with status '"
                                  << worhp[i].cnt.status << "'!" << std::endl;
                    }
                    else
                    {
                        step_length = it_point->GetDistance(worhp[i].opt.X);
                        if (step_length < parameters.SPREAD_MIN_SEARCH_LENGTH)
                        {
                            // search length too small, TODO: go into feasibility restoration
                            std::cout << "SpreadParetoFront: Search length too small!" << std::endl;
                        }
                        else
                        {
                            // everything fine
                            x.assign(worhp[i].opt.X, worhp[i].opt.X + monlp.GetNumVariables());
                            lambda.assign(worhp[i].opt.Lambda, worhp[i].opt.Lambda + monlp.GetNumVariables());
                            mu.assign(worhp[i].opt.Mu, worhp[i].opt.Mu + monlp.GetNumConstraints());
                            double *penalty = RWS_PTR((&worhp[i].wsp), worhp[i].wsp.penalty);
                            penalties.assign(penalty, penalty + monlp.GetNumConstraints());
                            new_points.emplace_back(x, lambda, mu, penalties, worhp[i].wsp.MeritNewValue, monlp);
                        }
                    }
                }

                it_point->SetStopped(true);
                // they do this for some reason, but why remove points just because they are
                // infeasible when they will be thrown out by cleanup method anyways. These
                // might still be decent starting points for the refinement stage.
                if (!it_point->IsFeasible())
                {
                    it_point = paretoFront.RemovePoint(it_point);
                    continue;
                }
            }

            it_point += 1;
        }

        int num_added = paretoFront.AddPoints(new_points);
        std::cout << "SpreadParetoFront: Added " << num_added << " points" << std::endl;
        paretoFront.WriteF(log);
        if (paretoFront.AllStopped())
        {
            std::cout << "SpreadParetoFront: All points stopped!" << std::endl;
            break;
        }
    }
}

void MOSQP::AddExtremeParetoPoints()
{
    double constexpr POS_INF = std::numeric_limits<double>::infinity();
    double constexpr NEG_INF = -POS_INF;

    // initial guess
    size_t const num_variables = monlp.GetNumVariables();
    std::vector<double> x0(num_variables);
    std::vector<double> x_l = monlp.GetXL();
    std::vector<double> x_u = monlp.GetXU();
    for (size_t i = 0; i < num_variables; i += 1)
    {
        if (x_l[i] > NEG_INF && x_u[i] < POS_INF)
        {
            x0[i] = (x_l[i] + x_u[i]) / 2;
        }
        else if (x_l[i] > NEG_INF)
        {
            x0[i] = x_l[i];
        }
        else if (x_u[i] < POS_INF)
        {
            x0[i] = x_u[i];
        }
        else
        {
            x0[i] = 0;
        }
    }

    // find extreme points
    std::vector<double> new_coordinates;
    size_t const num_objectives = monlp.GetNumObjectives();
    for (size_t objective_index = 0; objective_index < num_objectives; objective_index += 1)
    {
        SingleMONLP problem(monlp, objective_index);
        std::unique_ptr<WorhpSolver> worhp(std::make_unique<WorhpSolver>(problem));
        worhp->SetInitialGuess(x0);
        worhp->Solve();

        if (worhp->cnt.status >= TerminateSuccess)
        {
            new_coordinates.assign(worhp->opt.X, worhp->opt.X + monlp.GetNumVariables());
            Point extreme_point(new_coordinates, monlp);
            scalings.push_back(1.0 + std::abs(extreme_point.GetObjectiveValue(objective_index)));
            paretoFront.AddPoint(extreme_point);
        }
        else
        {
            scalings.push_back(1.0);
        }
    }

    paretoFront.WriteF(log);
}

void MOSQP::RefineParetoFront()
{
    paretoFront.UnstopAll();

    CombinedMONLP combinedProblem(monlp, scalings);
    std::unique_ptr<WorhpSolver> worhp(std::make_unique<WorhpSolver>(combinedProblem));
    worhp->par.ArmijoMinAlpha = parameters.REFINE_ARMIJO_MIN_ALPHA;
    worhp->par.ArmijoBeta = parameters.REFINE_ARMIJO_BETA;
    worhp->par.ArmijoBetaAres = parameters.REFINE_ARMIJO_BETA;
    worhp->par.KeepAcceptableSol = false;
    worhp->par.LowPassFilter = false;
    worhp->par.MaxIter = std::numeric_limits<int>::max();
    worhp->par.NLPmethod = 1;  // use merit function instead of filter
    worhp->par.TolFeas = 1e-20;
    worhp->par.TolOpti = 1e-20;

    std::vector<Point> new_points;
    std::vector<double> x;
    std::vector<double> lambda;
    std::vector<double> mu;
    std::vector<double> penalties;
    double step_length;

    for (int step = 0; step < parameters.REFINE_MAX_STEPS; step += 1)
    {
        new_points.clear();
        for (auto it_point = paretoFront.begin(); it_point != paretoFront.end();)
        {
            if (!it_point->IsStopped())
            {
                // The feasibility restoration will be a bit different from what the paper does, so
                // maybe assert that direction is descent direction for all f.
                combinedProblem.SetParameters(it_point->GetObjectiveValues());
                worhp->SetInitialGuess(*it_point);
                worhp->DoMajorIter();

                if (worhp->cnt.status <= TerminateError)
                {
                    std::cout << "RefineParetoFront: WORHP terminated with status '"
                              << worhp->cnt.status << "'!" << std::endl;
                }
                else
                {
                    // only interested in point if worhp didnt terminate with an error
                    x.assign(worhp->opt.X, worhp->opt.X + combinedProblem.GetNumVariables());
                    lambda.assign(worhp->opt.Lambda, worhp->opt.Lambda + combinedProblem.GetNumVariables());
                    mu.assign(worhp->opt.Mu, worhp->opt.Mu + combinedProblem.GetNumConstraints());
                    double *penalty = RWS_PTR((&worhp->wsp), worhp->wsp.penalty);
                    penalties.assign(penalty, penalty + monlp.GetNumConstraints());
                    new_points.emplace_back(x, lambda, mu, penalties, worhp->wsp.MeritNewValue, monlp);
                    step_length = it_point->GetDistance(worhp->opt.X);

                    if (worhp->cnt.status >= TerminateSuccess)
                    {
                        new_points.back().SetStopped(true);
                        std::cout << "RefineParetoFront: Optimal point found, status '"
                                  << worhp->cnt.status << "'!" << std::endl;
                    }
                    else if (step_length < parameters.REFINE_MIN_SEARCH_LENGTH)
                    {
                        if (new_points.back().IsFeasible())
                        {
                            new_points.back().SetStopped(true);
                            std::cout << "RefineParetoFront: Optimal point found, search length small!" << std::endl;
                        }
                        else
                        {
                            new_points.erase(new_points.end() - 1);
                        }
                    }
                }

                if (it_point->IsFeasible())
                {
                    // point was used for finding descent direction regarding combined function
                    // so stop it, as we dont want to get the same descent direction again
                    it_point->SetStopped(true);
                }
                else
                {
                    // point was used for finding descent direction regarding combined function,
                    // so it will not be useful to find further points. Also it's not feasible and
                    // as we only want feasible points in the paretoFront throw it out.
                    it_point = paretoFront.RemovePoint(it_point);
                    continue;
                }
            }

            ++it_point;
        }

        int num_added = paretoFront.AddPoints(new_points);
        std::cout << "RefineParetoFront: Added " << num_added << " points" << std::endl;
        paretoFront.WriteF(log);
        if (paretoFront.AllStopped())
        {
            break;
        }
    }

    // remove infeasible points
    for (auto it_point = paretoFront.begin(); it_point != paretoFront.end();)
    {
        if (!it_point->IsFeasible())
        {
            it_point = paretoFront.RemovePoint(it_point);
        }
        else
        {
            ++it_point;
        }
    }

    paretoFront.WriteF(log);
}

} // namespace mosqp
