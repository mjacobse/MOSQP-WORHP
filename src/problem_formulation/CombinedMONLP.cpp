#include "CombinedMONLP.hpp"
#include "MatrixStructure.hpp"
#include "MONLP.hpp"
#include <algorithm>
#include <cstddef>
#include <limits>
#include <numeric>
#include <vector>


namespace mosqp
{

CombinedMONLP::CombinedMONLP(MONLP const &monlp, std::vector<double> const &scalings)
    : NLP(monlp.GetNumVariables(), monlp.GetNumConstraints() + monlp.GetNumObjectives(),
          GetCombinedStructureDF(monlp), GetCombinedStructureDG(monlp), GetCombinedStructureHM(monlp),
          monlp.UserDF(), false, monlp.UserHM(),
          monlp.GetXL(), monlp.GetXU(), monlp.GetGL(), monlp.GetGU()),
      monlp(monlp), scalings(scalings)
{
    combined_gl = monlp.GetGL();
    combined_gu = monlp.GetGU();
    combined_gl.insert(combined_gl.end(), monlp.GetNumObjectives(), -std::numeric_limits<double>::infinity());
    combined_gu.insert(combined_gu.end(), monlp.GetNumObjectives(), 0);
}

void CombinedMONLP::SetParameters(std::vector<double> parameters)
{
    NLP::SetParameters(parameters);
}

MatrixStructure CombinedMONLP::GetCombinedStructureDF(MONLP const &monlp)
{
    MatrixStructure combined_structure;
    for (size_t objective_index = 0; objective_index < monlp.GetNumObjectives(); objective_index += 1)
    {
        MatrixStructure const &structure = monlp.GetStructureDF(objective_index);
        if (!structure.IsDefined())
        {
            return MatrixStructure(false);
        }

        for (auto row = structure.GetRowBegin(); row != structure.GetRowEnd(); ++row)
        {
            combined_structure.SetEntryNNZ(*row, 1);
        }
    }

    return combined_structure;
}

MatrixStructure CombinedMONLP::GetCombinedStructureDG(MONLP const &monlp)
{
    if (!monlp.GetStructureDG().IsDefined())
    {
        return MatrixStructure(false);
    }

    MatrixStructure combined_structure = monlp.GetStructureDG();
    size_t const num_objectives = monlp.GetNumObjectives();
    size_t const num_constraints = monlp.GetNumConstraints();
    for (size_t objective_index = 0; objective_index < num_objectives; objective_index += 1)
    {
        MatrixStructure const &df_structure = monlp.GetStructureDF(objective_index);
        if (!df_structure.IsDefined())
        {
            return MatrixStructure(false);
        }
        
        for (auto row = df_structure.GetRowBegin(); row != df_structure.GetRowEnd(); ++row)
        {
            combined_structure.SetEntryNNZ(num_constraints + objective_index + 1, *row);
        }
    }

    return combined_structure;
}

MatrixStructure CombinedMONLP::GetCombinedStructureHM(MONLP const &monlp)
{
    size_t const num_variables = monlp.GetNumVariables();
    size_t const num_objectives = monlp.GetNumObjectives();

    std::vector<size_t> diag(num_variables);
    std::iota(diag.begin(), diag.end(), 1);
    MatrixStructure combined_structure(diag, diag);

    for (size_t objective_index = 0; objective_index < num_objectives; objective_index += 1)
    {
        MatrixStructure const &structure = monlp.GetStructureHM(objective_index);
        if (!structure.IsDefined())
        {
            return MatrixStructure(false);
        }

        std::vector<size_t>::const_iterator row = structure.GetRowBegin();
        std::vector<size_t>::const_iterator col = structure.GetColBegin();
        while (row != structure.GetRowEnd() && col != structure.GetColEnd())
        {
            if (*row != *col)
            {
                combined_structure.SetEntryNNZ(*row, *col, num_variables);
            }

            ++row;
            ++col;
        }
    }

    return combined_structure;
}

std::vector<double> const & CombinedMONLP::GetGL() const
{
    return combined_gl;
}

std::vector<double> const & CombinedMONLP::GetGU() const
{
    return combined_gu;
}

double CombinedMONLP::EvalF_impl(double const *const x) const
{
    double f = 0;
    for (size_t i = 0; i < monlp.GetNumObjectives(); i += 1)
    {
        f += monlp.EvalF(x, i) / scalings[i];
    }

    return f;
}

double CombinedMONLP::GetF(Point const &point) const
{
    double f = 0;
    for (size_t i = 0; i < monlp.GetNumObjectives(); i += 1)
    {
        f += point.GetObjectiveValue(i) / scalings[i];
    }

    return f;
}

void CombinedMONLP::EvalG_impl(double const *const x, double *const g) const
{
    monlp.EvalG(x, g);

    size_t const num_constraints = monlp.GetNumConstraints();
    size_t const num_objectives = monlp.GetNumObjectives();
    for (size_t i = 0; i < num_objectives; i += 1)
    {
        g[i + num_constraints] = monlp.EvalF(x, i) - parameters[i];
    }
}

void CombinedMONLP::GetG(Point const &point, double *g) const
{
    std::vector<double> const &original_g = point.GetConstraints();
    std::copy(original_g.cbegin(), original_g.cend(), g);

    size_t const num_constraints = monlp.GetNumConstraints();
    size_t const num_objectives = monlp.GetNumObjectives();
    for (size_t i = 0; i < num_objectives; i += 1)
    {
        g[i + num_constraints] = point.GetObjectiveValue(i) - parameters[i];
    }
}

void CombinedMONLP::EvalDF_impl(double const *const x, double *const df) const
{
    size_t const num_objectives = monlp.GetNumObjectives();
    MatrixStructure structure_df = GetStructureDF();
    std::vector<size_t>::const_iterator df_row = structure_df.GetRowBegin();
    std::fill(df, df + structure_df.GetNumNonZeros(), 0.0);

    size_t i;
    for (size_t objective_index = 0; objective_index < num_objectives; objective_index += 1)
    {
        i = 0;

        MatrixStructure orig_structure_df = monlp.GetStructureDF(objective_index);
        std::vector<double> df_orig_val_vector(orig_structure_df.GetNumNonZeros());
        monlp.EvalDF(x, df_orig_val_vector.data(), objective_index);

        std::vector<size_t>::const_iterator df_orig_row = orig_structure_df.GetRowBegin();
        std::vector<size_t>::const_iterator df_orig_row_end = orig_structure_df.GetRowEnd();
        std::vector<double>::const_iterator df_orig_val = df_orig_val_vector.cbegin();

        while (df_orig_row != df_orig_row_end)
        {
            if (*df_orig_row == df_row[i])
            {
                df[i] += *df_orig_val / scalings[objective_index];
                ++df_orig_row;
                ++df_orig_val;
            }

            i += 1;
        }
    }
}

/*void CombinedMONLP::EvalDG_impl(double const *const x, double *const dg) const
{
    monlp.EvalDG(x, dg);
    size_t const num_objectives = monlp.GetNumObjectives();
    size_t current_index = monlp.GetStructureDG().col.size();

    for (size_t column = 0; column < num_variables; column += 1)
    {
        for (size_t objective_index = 0; objective_index < num_objectives; objective_index += 1)
        {
            // TODO: fix this so that its in columnwise order
            monlp.EvalDF(x, dg + current_index, objective_index);
            current_index += monlp.GetStructureDF(objective_index).row.size();
        }
    }
}*/

void CombinedMONLP::EvalHM_impl(double const *x, double const *mu, double scale_obj, double *hm) const
{
    size_t const num_objectives = monlp.GetNumObjectives();
    size_t const num_constraints = monlp.GetNumConstraints();
    MatrixStructure structure_hm = GetStructureHM();
    std::vector<size_t>::const_iterator hm_row = structure_hm.GetRowBegin();
    std::vector<size_t>::const_iterator hm_col = structure_hm.GetColBegin();
    std::fill(hm, hm + structure_hm.GetNumNonZeros(), 0.0);

    size_t i;
    for (size_t objective_index = 0; objective_index < num_objectives; objective_index += 1)
    {
        i = 0;

        MatrixStructure orig_structure_hm = monlp.GetStructureHM(objective_index);
        std::vector<double> df2_orig_val_vector(orig_structure_hm.GetNumNonZeros());
        monlp.EvalD2F(x, df2_orig_val_vector.data(), objective_index);

        std::vector<size_t>::const_iterator d2f_orig_row = orig_structure_hm.GetRowBegin();
        std::vector<size_t>::const_iterator d2f_orig_col = orig_structure_hm.GetColBegin();
        std::vector<double>::const_iterator d2f_orig_val = df2_orig_val_vector.cbegin();
        std::vector<double>::const_iterator d2f_orig_val_end = df2_orig_val_vector.cend();

        while (d2f_orig_val != d2f_orig_val_end)
        {
            if (*d2f_orig_row == hm_row[i] && *d2f_orig_col == hm_col[i])
            {
                hm[i] += *d2f_orig_val * (1.0 / scalings[objective_index] + mu[objective_index + num_constraints]);
                ++d2f_orig_row;
                ++d2f_orig_col;
                ++d2f_orig_val;
            }

            i += 1;
        }
    }

    MatrixStructure orig_structure_d2g = monlp.GetStructureHM(1);
    std::vector<double> dg_val_vector(orig_structure_d2g.GetNumNonZeros());
    monlp.EvalD2G(x, mu, dg_val_vector.data(), 1);

    std::vector<size_t>::const_iterator d2g_orig_row = orig_structure_d2g.GetRowBegin();
    std::vector<size_t>::const_iterator d2g_orig_col = orig_structure_d2g.GetColBegin();
    std::vector<double>::const_iterator d2g_orig_val = dg_val_vector.cbegin();
    std::vector<double>::const_iterator d2g_orig_val_end = dg_val_vector.cend();

    i = 0;
    while (d2g_orig_val != d2g_orig_val_end)
    {
        if (*d2g_orig_row == hm_row[i] && *d2g_orig_col == hm_col[i])
        {
            hm[i] += *d2g_orig_val;
            ++d2g_orig_row;
            ++d2g_orig_col;
            ++d2g_orig_val;
        }

        i += 1;
    }
}

} // namespace mosqp
