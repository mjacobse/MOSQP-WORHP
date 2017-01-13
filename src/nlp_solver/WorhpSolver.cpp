#include "WorhpSolver.hpp"
#include "../../include/worhp/worhp.h"
#include "../algorithm/Point.hpp"
#include "../problem_formulation/MatrixStructure.hpp"
#include "../problem_formulation/NLP.hpp"
#include <algorithm>
#include <cassert>
#include <cstddef>
#include <vector>


namespace mosqp
{

void silentPrint(int mode, const char *s)
{
}

WorhpSolver::WorhpSolver(NLP const &nlp)
    : nlp(nlp)
{
    //SetWorhpPrint(silentPrint);
    Init();
}

WorhpSolver::~WorhpSolver()
{
    StatusMsg(&opt, &wsp, &par, &cnt);
    WorhpFree(&opt, &wsp, &par, &cnt);
}

void WorhpSolver::DoMajorIter()
{
    int major_iter = wsp.MajorIter;
    bool terminated = false;
    while (!terminated && (wsp.MajorIter - major_iter) < 20 &&
           (GetPreviousStage(&cnt, -1) != Finalise || wsp.CurrentFeasMode != 0))
    {
        terminated = !Loop();
    }
}

void WorhpSolver::Solve()
{
    while (Loop())
    {
    }
}

void WorhpSolver::SetInitialGuess(std::vector<double> const &x)
{
    std::copy(x.begin(), x.end(), opt.X);
    if (!(GetCurrentStage(&cnt) == Init_Data))
    {
        SetNextStage(&cnt, Pre_KKT);
        cnt.status = Iterating;
    }
}

void WorhpSolver::SetInitialGuess(Point const &point)
{
    SetInitialGuess(point.GetX());
    opt.F = nlp.GetF(point);
    nlp.GetG(point, opt.G);
    if (point.HasMultipliers())
    {
        std::vector<double> const &lambda = point.GetLambda();
        std::vector<double> const &mu = point.GetMu();
        std::vector<double> const &penalties = point.GetPenalties();
        std::copy(lambda.cbegin(), lambda.cend(), opt.Lambda);
        std::copy(mu.cbegin(), mu.cend(), opt.Mu);
        std::copy(penalties.cbegin(), penalties.cend(), RWS_PTR((&wsp), wsp.penalty));
        wsp.MeritOldValue = point.GetMeritValue();
    }
}

bool WorhpSolver::Loop()
{
    if (GetUserAction(&cnt, callWorhp))
    {
        Worhp(&opt, &wsp, &par, &cnt);
    }

    if (GetUserAction(&cnt, iterOutput))
    {
        IterationOutput(&opt, &wsp, &par, &cnt);
        DoneUserAction(&cnt, iterOutput);
    }

    if (GetUserAction(&cnt, evalF))
    {
        opt.F = wsp.ScaleObj * nlp.EvalF(opt.X);
        DoneUserAction(&cnt, evalF);
    }

    if (GetUserAction(&cnt, evalDF))
    {
        nlp.EvalDF(opt.X, wsp.DF.val);
        for (size_t i = 0; i < wsp.DF.nnz; i += 1)
        {
            wsp.DF.val[i] *= wsp.ScaleObj;
        }
        DoneUserAction(&cnt, evalDF);
    }

    if (GetUserAction(&cnt, evalG))
    {
        nlp.EvalG(opt.X, opt.G);
        DoneUserAction(&cnt, evalG);
    }

    if (GetUserAction(&cnt, evalDG))
    {
        nlp.EvalDG(opt.X, wsp.DG.val);
        DoneUserAction(&cnt, evalDG);
    }

    if (GetUserAction(&cnt, evalHM))
    {
        nlp.EvalHM(opt.X, opt.Mu, wsp.ScaleObj, wsp.HM.val);
        DoneUserAction(&cnt, evalHM);
    }

    if (GetUserAction(&cnt, fidif))
    {
        WorhpFidif(&opt, &wsp, &par, &cnt);
    }

    return (cnt.status < TerminateSuccess && cnt.status > TerminateError);
}

void WorhpSolver::Init()
{
    WorhpPreInit(&opt, &wsp, &par, &cnt);

    int status = 0;
    char param_file[] = "../../../src/worhp.xml";
    ReadParams(&status, param_file, &par);

    opt.n = static_cast<int>(nlp.GetNumVariables());
    opt.m = static_cast<int>(nlp.GetNumConstraints());

    MatrixStructure const &structure_df = nlp.GetStructureDF();
    MatrixStructure const &structure_dg = nlp.GetStructureDG();
    MatrixStructure const &structure_hm = nlp.GetStructureHM();

    wsp.DF.nnz = structure_df.IsDefined() ? static_cast<mat_int>(structure_df.GetNumNonZeros()) : WorhpMatrix_Init_Dense;
    wsp.DG.nnz = structure_dg.IsDefined() ? static_cast<mat_int>(structure_dg.GetNumNonZeros()) : WorhpMatrix_Init_Dense;
    wsp.HM.nnz = structure_hm.IsDefined() ? static_cast<mat_int>(structure_hm.GetNumNonZeros()) : WorhpMatrix_Init_Dense;

    par.UserDF = nlp.UserDF();
    par.UserDG = nlp.UserDG();
    par.UserHM = nlp.UserHM();

    WorhpInit(&opt, &wsp, &par, &cnt);
    assert(cnt.status == FirstCall);

    std::vector<double> const &xl = nlp.GetXL();
    std::vector<double> const &xu = nlp.GetXU();
    std::vector<double> const &gl = nlp.GetGL();
    std::vector<double> const &gu = nlp.GetGU();
    std::copy(xl.begin(), xl.end(), opt.XL);
    std::copy(xu.begin(), xu.end(), opt.XU);
    std::copy(gl.begin(), gl.end(), opt.GL);
    std::copy(gu.begin(), gu.end(), opt.GU);

    if (wsp.DF.NeedStructure)
    {
        std::copy(structure_df.GetRowBegin(), structure_df.GetRowEnd(), wsp.DF.row);
    }

    if (wsp.DG.NeedStructure)
    {
        std::copy(structure_dg.GetRowBegin(), structure_dg.GetRowEnd(), wsp.DG.row);
        std::copy(structure_dg.GetColBegin(), structure_dg.GetColEnd(), wsp.DG.col);
    }

    if (wsp.HM.NeedStructure)
    {
        std::copy(structure_hm.GetRowBegin(), structure_hm.GetRowEnd(), wsp.HM.row);
        std::copy(structure_hm.GetColBegin(), structure_hm.GetColEnd(), wsp.HM.col);
    }
}

} // namespace mosqp
