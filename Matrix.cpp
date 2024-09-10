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
    elements.resize(4);
    col_ind.resize(4);
    row_ind.resize(3);
    for (int x=0; x< row_ind.size(); x++)
    {
        row_ind[x] = 0;
    }

}

  /**
   * Parameterized constructor.  It should create an empty matrix of user-supplied dimensions.
   * @param m - number of row for the new matrix.
   * @param n - number of columns for the new matrix.
   */ 
Matrix::Matrix(std::size_t m, std::size_t n): n((n==0 || m==0)?0:n) 
{
    //add check for parameters = 0
   //inialize  
   if ( (m == 0) || (n == 0))
   {
       row_ind.resize(1);
       elements.resize(0);
       col_ind.resize(0);
   }
   else
   {
       row_ind.resize(m+1);
       col_ind.resize(n*m);
       elements.resize(n*m);

        for (int x=0; x< row_ind.size(); x++)
        {
            row_ind[x] = 0;
        }
   }
//if m or n is zero then row resize should be 1
}
  
  /**
   * Parameterized constructor.  Use the parameters to set the matrix elements; if parameters are inconsistent then create a 0-by-0 (empty) matrix.
   * @param A - values for matrix elements, specified using a dense, row-major order vector
   * @param n - number of columns for the new matrix.
   */ 
Matrix::Matrix(const std::vector<Elem> &A, std::size_t n): n((n==0)?0:n)
{
    int vecLength = A.size();
    
//if or n is zero then row resize should be 1
    if ((n != 0) )
    {
        int rowSize = A.size() / n + 1 ;    // number of rows plus one
        
         //initialize the private vectors
             row_ind.resize(rowSize);
             col_ind.resize(vecLength);
             elements.resize(vecLength);

        for (int x=0; x< row_ind.size(); x++)
        {
            row_ind[x] = 0;
        }

        for (int i = 0; i < vecLength; i++)
        {
                e(i/n, i%n, A[i]);
        }
    }
    else
    {

     //initialize the private vectors
         row_ind.resize(1);
         col_ind.resize(0);
         elements.resize(0);

    }
    //iterate from zero to size of A

}

  /**
   * Another parameterized constructor.  Use the parameters to set the matrix element; if parameters are inconsistent then create a 0-by-0 (empty) matrix.
   * @param ptr_A - pointer to a dense, row-major order array of matrix elements
   * @param m - number of rows for the new matrix.
   * @param n - number of columns for the new matrix.
   */ 
Matrix::Matrix(const Elem *ptr_A, std::size_t m, std::size_t n): n((n==0 || m==0)?0:n)
{
    //m*n 

        
       if ( (m == 0) || (n == 0))
       {
           row_ind.resize(1);
           elements.resize(0);
           col_ind.resize(0);
       }
       else
       {
           row_ind.resize(m+1);
           elements.resize(m*n);
           col_ind.resize(m*n);

            for (int x=0; x< row_ind.size(); x++)
            {
                row_ind[x] = 0;
            }


            for (int i = 0; i < sizeof(ptr_A); i++)
            {
                    e(i/n, i%n,ptr_A[i]);
            }

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
    Elem test = 0;

    //Check that the index is valid
    if ( (i<0) || (j<0)|| (i >= row_ind.size()-1) || (j >= n))
    {
        return test;
    }
    // Find the row index of the value
    

    //Iterate over col_ind from index start to end
    for (int x = row_ind[i]; x < row_ind[i+1]; x++)
    {
        if (col_ind[x] == j)
        {
            return elements[x];
        }
    }
        return test;
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
    Elem test = false;

    //Check that the index is valid
    if ( (i<0) || (j<0)|| (i >= row_ind.size()-1) || (j > n))
    {
        return test;
    }
    if ( e(i,j) > 0)
    {
        for (int x = row_ind[i]; x < row_ind[i+1]; x++)
        {
            if (col_ind[x] == j)
            {
                elements[x] = aij;
                return true;
            }
        }

    } 
    else
    {

        for (int x = i; x < row_ind.size(); x++)
        {
            row_ind[i+1] = row_ind[i+1] + 1;
        }
        col_ind.insert(col_ind.begin()+row_ind[i], j);
        col_ind.pop_back();
        elements.insert(elements.begin()+row_ind[i], aij);
        elements.pop_back();
        
        return true;
        
    }

    
        return test;
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
    Matrix result;
  return result;
}

  /**
   * Creates and returns a new Matrix object representing the matrix subtraction of two Matrix objects.
   * @return a new Matrix object that contains the appropriate difference elements, a 0-by-0 matrix if matrices can't be subtracted.
   * @param rhs - the Matrix object to subtract from this object.
   */
  Matrix Matrix::sub( const Matrix &rhs ) const
{
    
    Matrix result;
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
  
  for(std::size_t i = 0; i < A.elements.size(); i++)
    {
      //beginning of a row (column zero): insert a line break unless we're in the first row
      if(i != 0 && i % A.n == 0)
	os << std::endl;
      
      os << A.elements[i] << " ";
    }

  return os;
}
