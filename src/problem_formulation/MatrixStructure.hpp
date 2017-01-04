#pragma once

#include <cstddef>
#include <vector>


namespace mosqp
{

class MatrixStructure
{
public:
    MatrixStructure(bool defined = true);
    MatrixStructure(std::vector<size_t> const &row);
    MatrixStructure(std::vector<size_t> const &row, std::vector<size_t> const &col);

    void SetEntryNNZ(size_t new_row, size_t new_col, size_t num_special_diag = 0);

    bool IsDefined() const;
    size_t GetNumNonZeros() const;
    std::vector<size_t>::const_iterator GetRowBegin() const;
    std::vector<size_t>::const_iterator GetRowEnd() const;
    std::vector<size_t>::const_iterator GetColBegin() const;
    std::vector<size_t>::const_iterator GetColEnd() const;
private:
    bool const defined;
    size_t nnz;
    std::vector<size_t> row;
    std::vector<size_t> col;

    bool IsSortingCorrect(size_t num_special_diag = 0);
};

} // namespace mosqp
