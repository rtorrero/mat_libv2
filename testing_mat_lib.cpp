// testing_mat_lib.hpp
// author: Antonio C. Domínguez Brito <antonio.dominguez@ulpgc.es>
// creation date: october 3rd 2020
// Description: This is a problem for testing mat_lib library

#include <iostream>
using namespace std;

#include "matrix.hpp"
#include "vector.hpp"

using matrix_t=mat_lib::matrix<double>;
using vector_t=mat_lib::vector<double>;
using complex_matrix_t=mat_lib::matrix<complex<double>>;

int main(int arc, char* argv[])
try
{
  cout<<"This is testing_mat_lib!"<<endl;

  matrix_t a{
    {1,  2,  3,  4},
    {5,  6,  7,  8},
    {9, 10, 11, 12}
  };
  cout<<"a="<<a<<endl;

  matrix_t b{a};
  cout<<"b="<<b<<endl;

  a+=b;
  cout<<"a="<<a<<endl;

  matrix_t  c=3.5*(a-b);
  cout<<"c="<<c<<endl;

  matrix_t d=a*3-b;
  cout<<"d="<<d<<endl;

  matrix_t e=-a;
  cout<<"e="<<e<<endl;

  matrix_t f=a.make_transpose();
  cout<<"f="<<f<<endl;

  a/=3;
  cout<<"a="<<a<<endl;

  matrix_t a2=a*a.make_transpose();
  cout<<"a2="<<a2<<endl;

  matrix_t l{1, 2, 3, 4, 5, 6};
  cout<<"l="<<l<<endl;

  cout<<"l*lt="<<l*~l<<endl;

  cout<<"lt*l="<<~l*l<<endl;

  matrix_t at=~a;
  cout<<"at="<<at<<endl;
  cout<<"at*a="<<~a*a<<endl;
  cout<<"a*at="<<a*~a<<endl;

  a.save_as("a.matrix");
  at.save_as("at.matrix");
  (~l*l).save_as("ltl.matrix");

  matrix_t g{"ltl.matrix"};
  cout<<"g="<<g<<endl;

  cout<<~matrix_t{"at.matrix"}<<endl;

  // Complex matrix tests
  cout << "\n=== Complex Matrix Tests ===\n";

  complex_matrix_t cm{
    {{1.0, 2.0}, {2.0, -1.0}, {3.0, 1.0}},
    {{0.0, 1.0}, {4.0, 0.0}, {2.0, -2.0}}
  };
  cout << "Complex matrix cm=" << cm << endl;

  // Save complex matrix
  cm.save_as("complex.matrix");

  // Load complex matrix
  complex_matrix_t cm_loaded{"complex.matrix"};
  cout << "Loaded complex matrix=" << cm_loaded << endl;

  // Test equality
  cout << "Complex matrix: Original equals loaded: " << boolalpha << (cm == cm_loaded) << endl;

  cout << "Double matrix: Original equals loaded: " << boolalpha << (at==matrix_t{"at.matrix"}) << endl;

  vector_t v{1, 2, 3, 4};
  cout << "v=" << v << endl;

  // Scalar multiplication
  vector_t v_scaled = v * 2.5;
  cout << "v * 2.5=" << v_scaled << endl;

  v_scaled /= 2.5;
  cout << "v_scaled / 2.5=" << v_scaled << endl;

  // Negation operator
  vector_t v_neg = -v;
  cout << "-v=" << v_neg << endl;

  // Element access
  cout << "v[2]=" << v.at(2) << endl;

  try {
    cout << "v[10] (out of range)=" << v.at(10) << endl;
  } catch (const exception &e) {
    cout << "[EXPECTED EXCEPTION] " << e.what() << endl;
  }

  // Vector addition and subtraction
  vector_t v_add = v + vector_t{4, 3, 2, 1};
  cout << "v + {4, 3, 2, 1}=" << v_add << endl;

  vector_t v_sub = v - vector_t{1, 1, 1, 1};
  cout << "v - {1, 1, 1, 1}=" << v_sub << endl;

  // Dot product
  vector_t v2{5, 6, 7, 8};
  double dot_product = dot(v, v2);
  cout << "dot(v, v2)=" << dot_product << endl;

  // Save and load a vector from file
  v.save_as("v.vector");
  vector_t v_loaded{"v.vector"};
  cout << "v_loaded=" << v_loaded << endl;

  // --- Matrix-vector multiplication tests ---
  // Matrix * Vector (column vector multiplication)
  vector_t mv = a * v;
  cout << "a * v=" << mv << endl;

  // Vector * Matrix (row vector multiplication)
  vector_t vm = vector_t{1, 2, 3} * a;
  cout << "vector{1, 2, 3} * a=" << vm << endl;

  // --- Advanced matrix tests ---
  // Test get_row and get_column
  cout << "Row 1 of a: " << a.get_row(1) << endl;
  cout << "Column 2 of a: " << a.get_column(2) << endl;

  // Check out-of-range access for get_row and get_column
  try {
    cout << "Invalid row of a: " << a.get_row(10) << endl;
  } catch (const exception &e) {
    cout << "[EXPECTED EXCEPTION] " << e.what() << endl;
  }

  try {
    cout << "Invalid column of a: " << a.get_column(5) << endl;
  } catch (const exception &e) {
    cout << "[EXPECTED EXCEPTION] " << e.what() << endl;
  }
  return 0;
}
catch(exception& e)
{ cerr<<"[EXCEPTION] "<<e.what()<<endl; }
