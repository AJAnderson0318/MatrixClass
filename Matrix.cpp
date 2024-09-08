#include <algorithm>
#include <limits> // std::numeric_limits

#include "Matrix.hpp"
//your code here...

//You'll implement this method in part two of the project
Matrix Matrix::invMod(Scalar m) const
{
  //For part two of the project you will use the elementary row operations to peform
  //Gauss-Jordan Elimination to find the inverse of a matrix...you may start early, if you wish
  return Matrix();
}

//Instructor provided methods and functions below; do not modify
Matrix::Matrix(const Matrix &A): n(A.n) 
{
  elements = A.elements;
  col_ind = A.col_ind;
  row_ind = A.row_ind;
}

Matrix& Matrix::operator=(Matrix A)
{
  std::swap(elements, A.elements);
  std::swap(col_ind, A.col_ind);
  std::swap(row_ind, A.row_ind);

  return *this;
}
  
std::ostream& operator<<(std::ostream& os, const Matrix &A)
{
  std::vector<Elem> A_d; //dense representation of the Matrix A
  std::size_t r_ind; //index in col_ind/elements where row i begins
  std::size_t r_cnt; //number of entries in col_ind/elements for the i-th row
  std::size_t i, j; //for looping
  
  A_d.resize((A.row_ind.size() - 1) * A.n); //matrix of zeros with size m * n 

  //add entries from sparse matrix to dense matrix
  for(i = 0; i < A.row_ind.size() - 1; i++) //iterate over rows
    {
      r_ind = A.row_ind[i]; //index in elements/col_ind for first entry in row i
      r_cnt = A.row_ind[i+1] - A.row_ind[i]; //number of elements in row i

      for(j = r_ind; j < r_ind + r_cnt; j++) //iterate over elements in row
	A_d[ (i * A.n) + A.col_ind[j] ] = A.elements[j]; //convert (i,j) to row-major order index
    }
 
  //output matrix
  for(i = 0; i < A_d.size(); i++)
    {
      //beginning of a row (column zero), insert a line break, unless we're in the first row
      if(i != 0 && i % A.n == 0)
	os << std::endl;
      
      os << A_d[i] << " ";
    }

  return os;
}
