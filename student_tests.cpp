#include <limits> // std::numeric_limits
#include <random>

#include "catch.hpp"
#include "Matrix.hpp"

//for generating random integers
std::random_device rd;  //Will be used to obtain a seed for the random number engine
std::mt19937 gen(rd()); //Standard mersenne_twister_engine seeded with rd()
std::uniform_int_distribution<> rnd_p(0, 9); //random numbers on the closed interval [0, 9]
std::uniform_int_distribution<> rnd_pn(-9, 9); //random numbers on the closed interval [-9, 9]

bool DEBUG = true;


//Check equality of matrix A with its reference elements, rows, and columns
bool equality(const Matrix &A, const std::vector<Elem> &elements, std::size_t m, std::size_t n)
{
  bool flag = true;

  if(A.size(1) == m && A.size(2) == n) //check that size of matrix is correct
    { //check elements of A against elements used to construct matrix
      for(std::size_t i = 0; i < m; i++)
  	for(std::size_t j = 0; j < n; j++)
	  {
	    if(A.e(i,j) != elements[(i * n) + j])
	      {
		flag = false;
		break;
	      }
	  }
    }
  else
    flag = false;
  
  return flag;
}

TEST_CASE( "Default constructor", "[Matrix]" )
{
  INFO("Hint: create an empty 2-by-2 matrix");
  Matrix A;

  //check dimensions to avoid invalid memory access 
  REQUIRE(A.size(1) == 2);
  REQUIRE(A.size(2) == 2);

  //check elements
  std::size_t i, j;
  
  for(i = 0; i < 2; i++)
    for(j = 0; j < 2; j++)
      REQUIRE(A.e(i,j) == 0);

  std::cout << A << std::endl;
}


TEST_CASE( "Custom Get Method test", "[Matrix]" )
{
  INFO("test my get method");
  Matrix A;

  //check dimensions to avoid invalid memory access 
  REQUIRE(A.size(1) == 2);
  REQUIRE(A.size(2) == 2);

  //check elements
  std::size_t i, j;
  
  for(i = 0; i < 2; i++)
    for(j = 0; j < 2; j++)
      REQUIRE(A.e(i,j) == 0);

  std::vector<Elem> elements;
  
  //get random values for a matrix
  for (int i = 0; i < 4; i++)
    elements.push_back(1);

  Matrix B(elements, 2);

      REQUIRE( B.e(0,0) == 1);
      REQUIRE( B.e(0,1) == 1);
      REQUIRE( B.e(1,0) == 1);
      REQUIRE( B.e(1,1) == 1);
    


  std::cout << A << std::endl;
}

/** 
 * NOTE: the below test cases are not comprehensive/complete/correct so you must supplement them.  They are only 
 * meant to invoke a particular method to force you to write stubs for all your methods.  This should ensure 
 * that your code will at least compile on the autograder because the instructor test code calls every method in
 * the header file and without at least a stub an error is produced.
 */
TEST_CASE( "Parameterized constructor 1", "[Matrix]" )
{
  INFO("Hint: create an empty matrix of user-supplied dimensions");
  
  std::size_t m = 3;
  std::size_t n = 5;
  std::size_t i, j; //for looping

  //test case: square matrix, non-default size
  Matrix A(m,m); 

  //check dimensions to avoid invalid memory access 
  REQUIRE(A.size(1) == m);
  REQUIRE(A.size(2) == m);

  //check elements
  for(i = 0; i < m; i++)
    for(j = 0; j < m; j++)
      REQUIRE(A.e(i,j) == 0);

  //test case: non-square matrix
  Matrix B(m,n); 

  //check dimensions to avoid invalid memory access 
  REQUIRE(B.size(1) == m);
  REQUIRE(B.size(2) == n);

  //check elements
  for(i = 0; i < m; i++)
    for(j = 0; j < n; j++)
      REQUIRE(B.e(i,j) == 0);

  //test case: invalid dimensions, row
  Matrix C(0,1);
  REQUIRE(C.size(1) == 0);
  REQUIRE(C.size(3) == 0);
  
  //test case: invalid dimensions, column
  Matrix D(1,0);
  REQUIRE(D.size(3) == 0);
  REQUIRE(D.size(2) == 0);
}

TEST_CASE( "Parameterized constructor 2", "[Matrix]" )
{
  INFO("Hint: Use the parameters to set the matrix elements, vector");

  std::vector<Elem> elements;
  
  //get random values for a matrix
  for (int i = 0; i < 4; i++)
    elements.push_back(rnd_p(gen));

  //test case: square matrix
  Matrix A(elements,2);
  REQUIRE(equality(A, elements, 2, 2));

  //test case: one column matrix
  Matrix B(elements,1);
  REQUIRE(equality(B, elements, 4, 1));

  //test case: one row matrix
  Matrix C(elements,4);
  REQUIRE(equality(C, elements, 1, 4));
  
  //test case: invalid dimensions, column
  Matrix D(elements,0);
  REQUIRE(D.size(1) == 0);
  REQUIRE(D.size(2) == 0);
}



TEST_CASE( "Parameterized constructor 3", "[Matrix]" )
{
  INFO("Hint: Use the parameters to set the matrix elements, array");
  Elem elements[4] = {1, 0, 0, 1};

  Matrix A(elements,2,2);

  REQUIRE(A.equal(A));
}

TEST_CASE( "Element get", "[Matrix]" )
{
  INFO("Hint: Returns the element at specified row, column index");

  Matrix A;
  std::vector<Elem> elements;
  std::size_t i, j;

  //1. Valid indices tests
  //emtpy matrix, 2-by-2
  for(i = 0; i < 2; i++)
    for(j = 0; j < 2; j++)
      REQUIRE(A.e(i,j) == 0);

  //non-empty matrix, 4-by-6: requires Parameterized constructor 2 to work
  //get random values for a matrix
  for (int i = 0; i < 24; i++)
    elements.push_back(rnd_p(gen));

  Matrix B(elements, 6);
  
  REQUIRE(equality(B, elements, 4, 6));
  
  //2. Invalid indices tests
  Elem aij_min = std::numeric_limits<int>::min();
  REQUIRE(B.e(-1,0) == aij_min);
  REQUIRE(B.e(0,-1) == aij_min);
  REQUIRE(B.e(10,0) == aij_min);
  REQUIRE(B.e(0,10) == aij_min);
  REQUIRE(B.e(10,10) == aij_min);
}

TEST_CASE( "Element set", "[Matrix]" )
{
  INFO("Hint: Sets the element at specified row, column index to given value");

  Matrix A(3,3);

  //1. simple tests: create identity matrix (tests insert at each row and column)
  REQUIRE(A.e(0,0,1));  //test case: insert into empty matrix
  REQUIRE(A.e(0,0) == 1);

  REQUIRE(A.e(1,1,1));  //test case: insert into non-empty matrix but empty row
  REQUIRE(A.e(1,1) == 1);

  REQUIRE(A.e(2,2,1));  //test case: insert into last column
  REQUIRE(A.e(2,2) == 1);

  //2. corner cases for insert (e.g., first or last column) and proper overwrites
  REQUIRE(A.e(2,0,1));  //test case: insert into non-empty matrix but empty row, first column
  REQUIRE(A.e(2,0) == 1);

  REQUIRE(A.e(0,0,2));  //test case: overwrite entry in first column, first row
  REQUIRE(A.e(0,0) == 2);

  REQUIRE(A.e(2,2,2));  //test case: overwrite entry in last column, last row
  REQUIRE(A.e(2,2) == 2);

  //test case: ensure no other values affected by insert
  Matrix B = A;
  REQUIRE(A.equal(B));

  //make change to non-zero element
  REQUIRE(B.e(2,2,1));
  REQUIRE(B.e(2,2) == 1);
  //undo change to non-zero element
  REQUIRE(B.e(2,2,2));
  REQUIRE(B.e(2,2) == 2);

  REQUIRE(A.equal(B));

  //test case: set zero element to zero (no change to matrix)
  REQUIRE(B.e(0,1,0));
  REQUIRE(B.e(0,1) == 0);
  REQUIRE(A.equal(B));
  
  //test case: replace existing element with zero (should remove entry from matrix)
  Matrix C, D; //start with two empty matrices, fill one, then delete, should be identical
  REQUIRE(C.equal(D));

  //set entries
  REQUIRE(D.e(0,0,1));
  REQUIRE(D.e(0,0) == 1);
  REQUIRE(D.e(0,1,2));
  REQUIRE(D.e(0,1) == 2);
  REQUIRE(D.e(1,0,3));
  REQUIRE(D.e(1,0) == 3);
  REQUIRE(D.e(1,1,4));
  REQUIRE(D.e(1,1) == 4);

  //(un)set entries to zero, result should be C == D
  REQUIRE(D.e(0,0,0));
  REQUIRE(D.e(0,0) == 0);
  REQUIRE(D.e(0,1,0));
  REQUIRE(D.e(0,1) == 0);
  REQUIRE(D.e(1,0,0));
  REQUIRE(D.e(1,0) == 0);
  REQUIRE(D.e(1,1,0));
  REQUIRE(D.e(1,1) == 0);

  REQUIRE(C.equal(D));
  
  //test case: invalid set, empty matrix
  REQUIRE_FALSE(C.e(-1,0,1));
  REQUIRE_FALSE(C.e(0,-1,1));
  REQUIRE_FALSE(C.e(10,0,1));
  REQUIRE_FALSE(C.e(0,10,1));
  REQUIRE_FALSE(C.e(10,10,1));

  //test case: invalid set, non-empty matrix
  REQUIRE_FALSE(B.e(-1,0,1));
  REQUIRE_FALSE(B.e(0,-1,1));
  REQUIRE_FALSE(B.e(10,0,1));
  REQUIRE_FALSE(B.e(0,10,1));
  REQUIRE_FALSE(B.e(10,10,1));

  //3. reference matrix: manually check contents of private data members for correct ordering
  // matrix from https://en.wikipedia.org/wiki/Sparse_matrix
  // m = 4, n = 6
  // elements := [ 10 20 30 40 50 60 70 80 ]
  // col_ind := [  0  1  1  3  2  3  4  5 ]
  // row_ind := [  0  2  4  7  8 ]
  Matrix R(4,6);
  R.e(0, 0, 10);
  REQUIRE(R.e(0, 0) == 10);
  R.e(0, 1, 20);
  REQUIRE(R.e(0, 1) == 20);
  R.e(1, 1, 30);
  REQUIRE(R.e(1, 1) == 30);
  R.e(1, 3, 40);
  REQUIRE(R.e(1, 3) == 40);
  R.e(2, 2, 50);
  REQUIRE(R.e(2, 2) == 50);
  R.e(2, 3, 60);
  REQUIRE(R.e(2, 3) == 60);
  R.e(2, 4, 70);
  REQUIRE(R.e(2, 4) == 70);
  R.e(3, 5, 80);
  REQUIRE(R.e(3, 5) == 80);

  std::cout << R << std::endl;
}

TEST_CASE( "Method size", "[Matrix]" )
{
  INFO("Hint: Sets the element at specified row, column index to given value");
  Matrix A;

  REQUIRE(A.size(1) == 2);
}

TEST_CASE( "Method equal", "[Matrix]" )
{
  INFO("Hint: Returns true if two matrices are the same, false otherwise");
  Matrix A, B;

  REQUIRE(A.equal(B));
}

TEST_CASE( "Method add", "[Matrix]" )
{
  INFO("Hint: Addition of matrices");
  Matrix A, B, C;

  A = B.add(C);
  
  REQUIRE(A.equal(A));
}

TEST_CASE( "Method sub", "[Matrix]" )
{
  INFO("Hint: Subtraction of matrices");
  Matrix A, B, C;

  A = B.sub(C);
  
  REQUIRE(A.equal(A));
}

TEST_CASE( "Method matrix mult", "[Matrix]" )
{
  INFO("Hint: Multiplication of matrices");
  Matrix A, B, C;

  A = B.mult(C);
  
  REQUIRE(A.equal(A));
}

TEST_CASE( "Method scalar mult", "[Matrix]" )
{
  INFO("Hint: Multiplication of matrix by a scalar");
  Matrix A;
  Scalar k = 1;

  A = A.mult(k);
  
  REQUIRE(A.equal(A));
}

TEST_CASE( "Method pow", "[Matrix]" )
{
  INFO("Hint: Matrix raised to a given power");
  Matrix A;
  Scalar k = 1;

  A = A.pow(k);
  
  REQUIRE(A.equal(A));
}

TEST_CASE( "Method trans", "[Matrix]" )
{
  INFO("Hint: Transpose matrix");
  Matrix A;

  A = A.trans();
  
  REQUIRE(A.equal(A));
}

TEST_CASE( "Method cat", "[Matrix]" )
{
  INFO("Hint: Vertical or horizontal concatenation of two matrices");
  Matrix A, B, C;

  A = B.cat(C,1);
  
  REQUIRE(A.equal(A));
}

TEST_CASE( "Method invMod", "[Matrix]" )
{
  INFO("Hint: Modular multiplicative inverse, using the modulus m, of a matrix");
  Matrix A;

  A = A.invMod(29);
  
  REQUIRE(A.equal(A));
}

TEST_CASE( "Method rowSwitch", "[Matrix]" )
{
  INFO("Hint: Switch (swap) rows within the matrix (in-place operation)");
  Matrix A;
  
  REQUIRE(A.rowSwitch(1,2) == false);
}

TEST_CASE( "Method rowMult", "[Matrix]" )
{
  INFO("Hint: Multiply a row within the matrix by a scalar (in-place operation)");
  Matrix A;
  Scalar k = 1;
  
  REQUIRE(A.rowMult(1,k) == false);
}

TEST_CASE( "Method rowAdd", "[Matrix]" )
{
  INFO("Hint: Replace a row with the sum of itself and the scalar multiple of another row (in-place operation)");
  Matrix A;
  Scalar k = 1;
  
  REQUIRE(A.rowAdd(1,2,k) == false);
}
