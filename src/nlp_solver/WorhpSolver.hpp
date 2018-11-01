#pragma once

#include "worhp/worhp.h"
#include "../algorithm/Point.hpp"
#include "../problem_formulation/NLP.hpp"
#include <vector>


namespace mosqp
{

class WorhpSolver
{
public:
    OptVar opt;
    Workspace wsp;
    Params par;
    Control cnt;

    WorhpSolver(NLP const &nlp);
    ~WorhpSolver();

    void DoMajorIter();
    void Solve();
    void SetInitialGuess(std::vector<double> const &x);
    void SetInitialGuess(Point const &point);
private:
    NLP const &nlp;
    bool Loop();
    void Init();
};

} // namespace mosqp
