#include <algorithm>
#include <limits> // std::numeric_limits

#include "Matrix.hpp"
//your code here...


  /**
   * Default constructor. It should create an empty 2-by-2 matrix.
   */ 
Matrix::Matrix():  n(2) 
{
    //row_ind need to be size of row plus 1
    row_ind = {0, 0, 0};
}

  /**
   * Parameterized constructor.  It should create an empty matrix of user-supplied dimensions.
   * @param m - number of row for the new matrix.
   * @param n - number of columns for the new matrix.
   */ 
Matrix::Matrix(std::size_t m, std::size_t n): n(n)
{
    for (int i = 0; i <= m; i++)
    {
        row_ind.push_back(0);
    }
}
  
  /**
   * Parameterized constructor.  Use the parameters to set the matrix elements; if parameters are inconsistent then create a 0-by-0 (empty) matrix.
   * @param A - values for matrix elements, specified using a dense, row-major order vector
   * @param n - number of columns for the new matrix.
   */ 
Matrix::Matrix(const std::vector<Elem> &A, std::size_t n): n(n)
{
    if ( (A.size()%n) == 0) //Checks that this combination of matrix size and columns can work
    { 
        
        int row_count = A.size() / n;
        
        for (int i = 0; i <= row_count+1; i++)
        {
            row_ind.push_back(0);
        }
        for (int i = 0; i < A.size(); i++)
        {
            if (A[i] != 0)
            {
                elements.push_back(A[i]);
                col_ind.push_back(i%n);
                row_ind[(i/n)+1] = row_ind[(i/n)+1] +1;
            }
        }


    }
}

  /**
   * Another parameterized constructor.  Use the parameters to set the matrix element; if parameters are inconsistent then create a 0-by-0 (empty) matrix.
   * @param ptr_A - pointer to a dense, row-major order array of matrix elements
   * @param m - number of rows for the new matrix.
   * @param n - number of columns for the new matrix.
   */ 
Matrix::Matrix(const Elem *ptr_A, std::size_t m, std::size_t n): n(n)
{

        for (int i = 0; i < m+1; i++)
        {
            row_ind.push_back(0);
        }
            
        int colIndex = 0;
        for (int i = 0; i < m*n; i++)
        {
            if (colIndex == n)
            {
                colIndex = 0;
            }
            if (ptr_A[i] != 0)
            {
             elements.push_back(ptr_A[i]);
             row_ind[i+1] = row_ind[i+1] + 1;
             col_ind.push_back(colIndex);
            }
            colIndex++;
        }
}

  
  /**
   * Returns the element at specified row, column index.
   * @param i - row index of object.
   * @param j - column index of object.
   * @return element at specified row, column index or smallest possible value for Elem if index is invalid.
   */ 
  Elem Matrix::e(std::size_t i, std::size_t j) const
{
    if ( (i > -1) && (j > -1)&& (i < row_ind.size()) && (j < n)) // check for valid indexes
    {
        for (int x = row_ind[i]; x < row_ind[i+1]; x++)
        {
            if (col_ind[x] == j)
            {
                return elements[x];
            }
        }
    }
        return 0;
}
  
  /**
   * Sets the element at specified row, column index to given value; if either index is invalid matrix should not be modified.
   * @param i - row index of object to set.
   * @param j - column index of object to set.
   * @param aij - value for element at index i, j
   * @return true if set is successful, false otherwise.
   */
  bool Matrix::e(std::size_t i, std::size_t j, Elem aij)
{
    //Check that the index is valid
    if ((i < row_ind.size()) && (i > -1) && ( j > n) && ( j > -1))
    {
        //if element exists then replace existing value in elements
        for (int x = row_ind[i]; x < row_ind[i+1]; i++)
        {
            if (col_ind[x] == j)
            {
                elements[x] = aij;
                return true;
            }
        }
        col_ind.insert(col_ind.begin()+row_ind[i+1], j);
        elements.insert(elements.begin()+row_ind[i+1], aij);
        row_ind[i+1] = row_ind[i+1] + 1;
        return true;

    }
    return false;
}

  /**
   * Returns the size of the matrix along a given dimension (i.e., number of row(s) or column(s))
   * @param dim - 1 for row, 2 for column
   * @return the length of the dimension specified, if dimension is not valid return 0
   */ 
  std::size_t Matrix::size(std::size_t dim) const
{
    int dimension = 0;
    if (dim == 1)
    {
        dimension = row_ind.size() - 1;
    }
    else if (dim == 2)
    {
        dimension = n;
    }
    return dimension;
}
  
 /**
   * Returns true if the matrices this and rhs are the same, false otherwise.
   * @param rhs - the Matrix object to compare to this object.
   * @return true if all the elements in both objects are the same, false otherwise.
   */ 
  bool Matrix::equal( const Matrix& rhs ) const
{
    bool ans = false;
    
    if ( (elements == rhs.elements) && (row_ind == rhs.row_ind) && (col_ind == rhs.col_ind) && (n == rhs.n))
    {
        ans = true;
    }

    return ans;
}

  /**
   * Creates and returns a new Matrix object representing the matrix addition of two Matrix objects.
   * @return a new Matrix object that contains the appropriate summed elements, a 0-by-0 matrix if matrices can't be added.
   * @param rhs - the Matrix object to add to this object.
   */
  Matrix Matrix::add( const Matrix &rhs ) const
{
    
    Matrix result = Matrix(0,0);
  if ( (row_ind.size() == rhs.row_ind.size()) && (n == rhs.n))
  {
    std::vector<Elem> newValues;

    for (int i = 0; i < row_ind.size(); i++)
    {
        for ( int j = 0; i < n; i++)
        {
            newValues.push_back( e(i,j) + rhs.e(i,j) );
        }
    }
    result = Matrix(newValues, n); 
  }
  return result;
}

  /**
   * Creates and returns a new Matrix object representing the matrix subtraction of two Matrix objects.
   * @return a new Matrix object that contains the appropriate difference elements, a 0-by-0 matrix if matrices can't be subtracted.
   * @param rhs - the Matrix object to subtract from this object.
   */
  Matrix Matrix::sub( const Matrix &rhs ) const
{
    
    Matrix result = Matrix(0,0);
  if ( (row_ind.size() == rhs.row_ind.size()) && (n == rhs.n))
  {
    std::vector<Elem> newValues;

    for (int i = 0; i < row_ind.size(); i++)
    {
        for ( int j = 0; i < n; i++)
        {
            newValues.push_back( e(i,j) - rhs.e(i,j) );
        }
    }
    result = Matrix(newValues, n); 
  }
  return result;
}

  /**
   * Creates and returns a new Matrix object that is the multiplication of this and the given Matrix object.
   * @return a new Matrix object that contains the multiplication of this and the given Matrix object, a 0-by-0 matrix if matrices can't be multiplied.
   * @param rhs - the Matrix object to multiply with this object.
   */
  Matrix Matrix::mult( const Matrix &rhs ) const
{
    Matrix ans = Matrix(0,0);

    if ( n == rhs.row_ind.size()-1) //if the matrices can be multiplied
    {
        Matrix mult = Matrix(row_ind.size() -1, rhs.n);


    }


    return ans;

}

  /**
   * Creates and returns a new Matrix object that is the multiplication of this and the given scalar.
   * @return a new Matrix object that contains the multiplication of this and the given scalar.
   * @param rhs - the scalar value to multiply with this object.
   */
  Matrix Matrix::mult( Scalar k ) const
{

    Matrix test;
    return test;
}

  /**
   * Creates and returns a new Matrix object that is the power of this.
   * @return a new Matrix object that raises this and to the given power.
   * @param k - the power to which this object should be raised, a 0-by-0 matrix if matrix can't be raised to power.
   */
  Matrix Matrix::pow( Scalar k ) const
{

    Matrix test;
    return test;
}

  /**
   * Creates and returns a new Matrix object that is the transpose of this.
   * @return a new Matrix object that is the transpose of this object.
   */
  Matrix Matrix::trans() const
{

    Matrix test;
    return test;
}

  /**
   * Creates and returns a new Matrix object that is the concatenation of this and the given Matrix object.
   * @return a new Matrix object that is the vertical or horizontal concatenation of two matrices, a 0-by-0 matrix if matrices can't be concatenated.
   * @param dim - 1 for vertical cat (RHS below this), 2 for horizontal cat (RHS to the right of this)
   */
  Matrix Matrix::cat(const Matrix &rhs, std::size_t dim) const
{

    Matrix test;
    return test;
}

  
  /**
   * Switch (swap) rows within the matrix (in-place operation)
   * Ri <=> Rj
   * @return true if rows i and j were switched
   * @param i - row number
   * @param j - row number
   */
  bool Matrix::rowSwitch(std::size_t i, std::size_t j)
{
    return false;
}

  /**
   * Multiply a row within the matrix by a scalar (in-place operation)
   * Ri = k*Ri
   * @return true if rows i was multiplied by scalar k
   * @param i - row number
   * @param k - scalar to multiply row i by
   */
  bool Matrix::rowMult(std::size_t i, Scalar k)
{
    return false;
}

  /**
   * Replace a row with the sum of itself and the scalar multiple of another row (in-place operation)
   * Ri = Ri + k*Rj
   * @return true if row i Ri
   * @param i - row number
   * @param k - scalar to multiply row i by
   */
  bool Matrix::rowAdd(std::size_t i, std::size_t j, Scalar k)
{
    return false;
}


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
