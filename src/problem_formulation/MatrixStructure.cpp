#include "MatrixStructure.hpp"
#include <cassert>
#include <cstddef>
#include <vector>


namespace mosqp
{

MatrixStructure::MatrixStructure(bool const defined)
    : defined(defined),  nnz(0), row(), col()
{
}

MatrixStructure::MatrixStructure(std::vector<size_t> const & row)
    : defined(true), nnz(row.size()), row(row), col(nnz, 1)
{
}

MatrixStructure::MatrixStructure(std::vector<size_t> const & row, std::vector<size_t> const & col)
    : defined(true), nnz(row.size()), row(row), col(col)
{
    assert(row.size() == col.size());
}

void MatrixStructure::SetEntryNNZ(size_t const new_row, size_t const new_col, size_t const num_special_diag)
{
    size_t index = 0;
    while (index < nnz - num_special_diag && new_col != col[index])
    {
        index += 1;
    }

    while (index < nnz - num_special_diag && new_col == col[index] && new_row > row[index])
    {
        index += 1;
    }

    if (index == nnz - num_special_diag || row[index] != new_row || col[index] != new_col)
    {
        row.insert(row.begin() + index, new_row);
        col.insert(col.begin() + index, new_col);
        nnz += 1;
    }

    assert(row.size() == nnz);
    assert(col.size() == nnz);
    assert(IsSortingCorrect(num_special_diag));
}

bool MatrixStructure::IsDefined() const
{
    return defined;
}

size_t MatrixStructure::GetNumNonZeros() const
{
    return nnz;
}

std::vector<size_t>::const_iterator MatrixStructure::GetRowBegin() const
{
    return row.cbegin();
}

std::vector<size_t>::const_iterator MatrixStructure::GetRowEnd() const
{
    return row.cend();
}

std::vector<size_t>::const_iterator MatrixStructure::GetColBegin() const
{
    return col.cbegin();
}

std::vector<size_t>::const_iterator MatrixStructure::GetColEnd() const
{
    return col.cend();
}

bool MatrixStructure::IsSortingCorrect(size_t num_special_diag)
{
    for (size_t i = 0; i < nnz - num_special_diag - 1; i += 1)
    {
        if (col[i] > col[i + 1])
        {
            return false;
        }
        else if (col[i] == col[i + 1])
        {
            if (row[i] >= row[i + 1])
            {
                return false;
            }
        }
    }

    return true;
}

} // namespace MOSQP
