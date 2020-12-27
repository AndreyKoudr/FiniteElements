/*
BSD 2-Clause License

Copyright (c) 2020, Andrey Kudryavtsev (andrewkoudr@hotmail.com)
All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice, this
   list of conditions and the following disclaimer.

2. Redistributions in binary form must reproduce the above copyright notice,
   this list of conditions and the following disclaimer in the documentation
   and/or other materials provided with the distribution.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*/

#include "Vector.h"
#include "Conforming3D.h"
#include "Matrix.h"
#include "FEMVirtMatrix.h"

#include <cstdlib>
#include <algorithm>
#include <functional>
#include <numeric>
#include <iostream>
#include <limits>
#include <chrono>

/**
  Full finite-element C++ code with SIMD-accelerated equation solver based on conforming linear 
finite elements with tests and example solution to Laplace equation. From scratch.
*/

/**
  Conforming finite elements
  --------------------------
  This is an example code (VS 2019 without Windows specifics) on how to write approximation and 
solution to a partial differential equation with conforming linear finite elements. Everything from
scratch, no third-party code.
  The element with all associated code is here in <I>Conforming3D.h</I>. "Conforming" means that you 
do not have to use a conformal mesh but build approximation on such constructions like balanced 
octrees. The element provides continuous approximation across element faces, but its shape functions 
are not trivial.

  Conforming linear finite element
  --------------------------------
  <I>K.-W. Cheng, T.-P. Fries. XFEM with Hanging Nodes for Two-phase 
Incompressible Flow</I>.

         (7)----------18-----------(6)
         /|                        /|    
        / |                       / |    
      19  |         25           17 |     
      /   |                     /   |          
     /   15           23       /    14         
   (4)-----------16----------(5)    |    
    |     |                   |     |     
    |  24 |                   |  22 |    
    |     |                   |     |
    |    (3)----------10------|----(2)
   12    /       21           13   / 
    |   /                     |   / 
    |  11          20         |  9
  W ^ /V                      | /
    |/                        |/
   (0)->---------8-----------(1)
       U

  Any node from 8 to 25 (18 nodes in all) may not exist. That is, 8 standard and
18 hanging. The element provides C0-continuous approximation across faces of 
neighbour elements. 
  Parameters U,V,W change from 0 to 1 from the node 0.

  T - class for real numbers
  Tvector - class for 3D vector with up to 4 components

  Derivatives
  -----------
    Derivatives are not continuous across parametric value of 0.5 due to e.g. 
  T c0a = T1 - std::abs(c0); in shapeFunctionsConforming() - abs value 
  function expressions, derivatives at and near hanging nodes with a parameter 0.5
  are not well defined as abs() is discontinuous, therefore the two functions 
  shapeDerivativesConformingSimple() based on finite differences and shapeDerivativesExact() - 
  exact expression MAY GIVE VERY DIFFERENT results. Because the derivative does not exist at x = 0 
  for function y = |x|.
    In general, first derivatives for a linear approximation must be constant, so the use
  of shapeDerivativesConformingFromSubOctant() which calculates derivatives at centres 
  of sub-octants looks most consistent.

  Matrix
  ------
  This class is for operations with not very large matrices but mainly for manipulations with
finite element shape functions. The class does not contain any special features for speedup,
it is a regular code for matrix ariphmetics and inversion for 2 x 2, 3 x 3 and 4 x 4 matrices.
These formulae are 100% reliable.

  Solution to a system of equations
  ---------------------------------
  A FEM system is banded and the code stores matrix by rows taking into account its banded nature with the 
solution by Gauss elemination. The built-in SIMD makes the code sufficiently faster. But multithreading does not
make any good. I do not know why, did not try to understand why. The maximum reasonable size of the 
system is something up to 70..100 thousands with 1/10 band width; the bigger solution will be destroyed by 
roundup errors, storage probems and time of solution. A better parallelisation and a more sophisticated 
algorithm is necessary for a good commercial code.

  Tests
  -----
  A group of tests, define/undefine TEST_MATRIX, TEST_ELEMENT, TEST_SOLUTION below.
  You can select various types for templates, all variants are possible :

  #define T float     - float variables are 4-byte floats
  #define T double    - floats are 8-bytes long

  #define Tvector Vector4           - SIMD-accelerated 4-component vector based on 4-byte floats
  #define Tvector TVector<float>    - plain 4-component vector on floats
  #define Tvector TVector<double>   - plain 4-component vector on doubles

  When testing solution to the system of equations, the matrix is not quite well positive definite,
so it would be much better to use doubles, like this :
  #define T double
  #define T TVector<double>

  If T and Tvector have different basic types, the code produces some type-conversion warnings, I left them
as they are.

  List of tests : disabled/enable by #define TEST_... :
  - Tests of Matrix class : ariphmetics, multiplies, inversion etc.
  - Tests of finite element shape functions : sum must be 1, for derivatives 0 etc. Some tests may fail
due to random nature of data when inconsistent with tolerance value currently set.
  - Tests of parallel solution to a banded system of equations. The matrix is not well positive definite,
so residuals maybe not very good with 4-byte float types. Another thing is that the increase in the 
number of threads does not help or even makes the code slower.
  - ... and finally finite element test. Mesh is not generated here; it is a matter of further publications.
It is as follows :


         (3)-----------------------(12)---------(21)--------(30) Z=1.0
         /|                        /|           /|          /|
        / |                       / |          / |         / |
       /  |                     (11)---------(20)--------(29)|
      /   |                     /|  |        /|  |       /|  |
     /    |                    / | (9)------/-|-(18)--- /-| (27)
   (2)-----------------------(10)---------(19)--------(28)| /|
    |     |                   |  |/ |      |  |/ |     |  |/ |
    |     |                   | (8)------- |-(17)----- | (26)|
    |     |                   | /|  |      | /|  |     | /|  |
    |    (1)------------------|/-|-(6)-----|/---(15)---|/-|-(24) Y=1.0
    |    /                   (7)----------(16)--------(25)| /
    |   /                     |  |/        |  |/       |  |/
    |  /                      | (5)--------|-(14)------|-(23)
  Z ^ /Y                      | /          | /         | /
    |/                        |/           |/          |/
   (0)->---------------------(4)----------(13)--------(22)
       X                     X=1.0                    X=2.0

- one big conforming finite element with five hanging nodes (5,7,8,9,11) and eight 
non-conforming small finite elements each with 8 nodes.
  The "non-conforming" plane is X = 0.5 to achieve continuous approximation
across it with our conforming finite element Confroming3D.

  The problem is formulated as follows : solve Laplace equation with essential boundary
conditions at plane X = 0 (potential is 0.0) and at plane X = 2.0 (specify potential 1.0).
On all other planes natural boundary conditions will be automatically satisfied. Read 
<I>J.J.Connor, C.A.Brebbia Finite Element Techniques for Fluid Flow</I>. Very good about weak 
formulations, stiffness and mass matrices etc.

  So we expect that potential value will grow linearly from 0.0 at X = 0 to 1.0 at X = 2 with
the derivative of potential (velocity) along X equal everywhere to 1.0 / 2.0 = 0.5 and zero velocity
components in Y and Z. 

  A test printout
  --------------
  Seed 1609068822

  ===== Tests of matrix =====

  Transpose

  Multiply (AB)T = BT AT

  Ariphmetic

  Inversion

  ===== Tests of 3D conforming finite element =====

  Testing incorrectly defined elements (arbitrary hanging nodes), shape function value at a node must be 1, at all other nodes 0

  Testing correctly defined elements (whole faces with hanging nodes), shape function value at a node must be 1, at all other nodes 0

  Calculation of shape function derivatives by four methods must yield close results

  Calculation of values at some specific cases (correct hanging nodes on two faces)

  Derivatives on real coordinates, element 2 x 0.5 x 0.75

  Derivatives on parameters and real coordinates, random tests :
(1) element derivatives on X,Y,Z must be units (identity matrix)
(2) element derivatives on parameters must be equal to box sizes
(3) element volume must be equal to Jacobi matrix determinant

  ===== Solution to banded system with SIMD and multiple threads  =====

Allocator stored file : FEMVirtMatrix.bin, size 40000000
Solving system of size 5000 by 500
1. solving with simple C solver...
4900 of 5000
2. solving with parallel solver...
4900
Time : simple solver 2.75633 sec, SIMD parallel solver 0.374472, acceleration 7.36059
Tesing results
Comparing solutions...
Max difference b/w solutions is 9.19168e-13, max solution value 0.109975, min value -0.110092, residual 8.99503e-13, parallel residual 1.33116e-12
Allocator stored file deleted : FEMVirtMatrix.bin

  ===== Solving Laplace equation with conforming finite elements =====

Allocator stored file : FEMVirtMatrix.bin, size 6944
0
Max residual 1.4988e-15
        X       Y       Z       potential       gradient

        0.5     0.5     0.5     0.25            0.5     5.55112e-17     -5.55112e-17
        1.25    0.25    0.25    0.625           0.5     -5.55112e-17    5.55112e-17
        1.25    0.75    0.25    0.625           0.5     -5.55112e-17    -1.66533e-16
        1.25    0.25    0.75    0.625           0.5     -1.11022e-16    0
        1.25    0.75    0.75    0.625           0.5     -2.77556e-17    -2.77556e-17
        1.75    0.25    0.25    0.875           0.5     0       0
        1.75    0.75    0.25    0.875           0.5     -1.11022e-16    -1.11022e-16
        1.75    0.25    0.75    0.875           0.5     -1.11022e-16    -1.11022e-16
        1.75    0.75    0.75    0.875           0.5     -1.11022e-16    -1.11022e-16
Allocator stored file deleted : FEMVirtMatrix.bin

  What is good and what is bad
  ----------------------------
  - Ready 100% reliable code for conforming finite elemnet. Remember that derivatives of shape functions are
discontinuous within the element.
  - Matrix class for small matrices with 100% reliable inversion formulae. If you run into a problem, 
you might have forgotten to transpose a matrix somewhere. Matrix operations are slow, no special attention 
to code speed.
  - Vector class for 4-component vectors Vector4 is accelerated by SIMD, but it cannot help much without
reorganisation of the whole code.
  - FEMVirtMatrix is SIMD accelerated with normal speedup 4 times or more, but multithreading does not help 
very much, maybe vecause of my computer is with just a single hardware set of XMM registers.
  - The last test - solution to a Laplace equation is very very simple but produces correct results. We need
a mesh generation code for proper testing (to come).

  References
  ----------
  - J.J.Connor, C.A.Brebbia. Finite Element Techniques for Fluid Flow (about finite elements).
  - K.-W. Cheng, T.-P. Fries. XFEM with Hanging Nodes for Two-phase Incompressible Flow (shape 
functions for conforming finite element). 
*/

using namespace std;

// Random integer
static int random(const int max)
{
  return (rand() - 1) * max / RAND_MAX;
}

// Random real
template <typename T> T random()
{
  return static_cast<T>(rand()) / static_cast<T>(RAND_MAX);
}

// Compute difference between two arrays of vectors
template <typename T, typename Tvector, size_t N> T diff(const std::array<Tvector,N> &a0, 
  const std::array<Tvector,N> &a1, const bool print = false)
{
#if 0
  T max = -std::numeric_limits<T>::max();
  int index = -1;
  Tvector sum0(0,0,0), sum1(0,0,0);
  for (int i = 0; i < N; i++)
  {
    T l = !(a0[i] - a1[i]);
    if (l > max)
    {
      max = l;
      index = i;
    }

    sum0 += a0[i];
    sum1 += a1[i];
  }
  if (print)
   cout << "max " << max << " index " << index << 
      " sum0 = " << sum0.X  << " " << sum0.Y << " " << sum0.Z << 
      " sum1 = " << sum1.X  << " " << sum1.Y << " " << sum1.Z << endl;
    return max;
#else
  std::array<T,N> l;
  std::transform(a0.begin(),a0.end(),a1.begin(),l.begin(),[](auto v0, auto v1) { return !(v0 - v1); });

  T sum = std::accumulate(l.begin(),l.end(),static_cast<T>(0.0),std::plus<T>());
  sum /= static_cast<T>(N);
  return sum;
#endif
}

// Create random matrix N * M
template <typename T, size_t N, size_t M> Matrix<T,N,M> randomMatrix(const T max)
{
  Matrix<T,N,M> m; 

  for (int i = 0; i < N; i++)
  {
    for (int j = 0; j < M; j++)
    {
#if 1
      m[i][j] = T(random(int(max)));
#else
      m[i][j] = random<T>() * max;
#endif
    }
  }

  return m;
}

template <typename T, typename Tvector> void addElementCoordinates(const Tvector& min, const Tvector& max,
  std::vector<Tvector> &totalCoordinates)
{
  // fill all coordinates including all hanging
  std::array<Tvector,Conforming3D<T,Tvector>::NUM_NODES> basicCoordinates = {
    Tvector(min.X,min.Y,min.Z),
    Tvector(max.X,min.Y,min.Z),
    Tvector(max.X,max.Y,min.Z),
    Tvector(min.X,max.Y,min.Z),

    Tvector(min.X,min.Y,max.Z),
    Tvector(max.X,min.Y,max.Z),
    Tvector(max.X,max.Y,max.Z),
    Tvector(min.X,max.Y,max.Z),

    Tvector(T(0.0),T(0.0),T(0.0)),
    Tvector(T(0.0),T(0.0),T(0.0)),
    Tvector(T(0.0),T(0.0),T(0.0)),
    Tvector(T(0.0),T(0.0),T(0.0)),
    Tvector(T(0.0),T(0.0),T(0.0)),
    Tvector(T(0.0),T(0.0),T(0.0)),
    Tvector(T(0.0),T(0.0),T(0.0)),
    Tvector(T(0.0),T(0.0),T(0.0)),
    Tvector(T(0.0),T(0.0),T(0.0)),
    Tvector(T(0.0),T(0.0),T(0.0)),
    Tvector(T(0.0),T(0.0),T(0.0)),
    Tvector(T(0.0),T(0.0),T(0.0)),
    Tvector(T(0.0),T(0.0),T(0.0)),
    Tvector(T(0.0),T(0.0),T(0.0)),
    Tvector(T(0.0),T(0.0),T(0.0)),
    Tvector(T(0.0),T(0.0),T(0.0)),
    Tvector(T(0.0),T(0.0),T(0.0)),
    Tvector(T(0.0),T(0.0),T(0.0))
  };

  for (int i = 0; i < Conforming3D<T,Tvector>::NUM_NODES; i++)
  {
    std::array<T,Conforming3D<T,Tvector>::NUM_NODES> func;
    Conforming3D<T,Tvector>::shapeFunctions(Conforming3D<T,Tvector>::nodeParms()[i],func);
    Tvector c = shapeFuncValue(func,basicCoordinates,Tvector(T(0.0),T(0.0),T(0.0)));
    totalCoordinates.push_back(c);
  }
}


#define TEST_MATRIX
#define TEST_ELEMENT
#define TEST_SOLUTION

int main()
{
  // can change it here into float/Vector4
  #define T double
  #define Tvector TVector<double>

  typedef Conforming3D<T,Tvector> Element;

  // start random generator
  auto seed = static_cast<unsigned int>(time(NULL));
  srand(seed);
  cout << "  Seed " << seed << std::endl << std::endl;

#ifdef TEST_MATRIX

  cout << "  ===== Tests of matrix =====" << endl << endl;

  cout << "  Transpose" << endl << endl;
  {
    T max = static_cast<T>(8.0);

    for (int i = 0; i < 100; i++)
    {
      Matrix<T,4,2> m = randomMatrix<T,4,2>(max);
      
      Matrix<T,2,4> mT = m.transpose();

      Matrix<T,4,2> mold = mT.transpose();

      assert(m == mold);
    }

    for (int i = 0; i < 100; i++)
    {
      Matrix<T,4,4> m = randomMatrix<T,4,4>(max);
      Matrix<T,4,4> mold = m;
      
      m.transposeSelf();

      assert(m != mold);
      
      m.transposeSelf();

      assert(m == mold);
    }

int gsggsgs = 0;
  }

  cout << "  Multiply (AB)T = BT AT" << endl << endl;
  {
    T max = static_cast<T>(18.0);

    for (int i = 0; i < 100; i++)
    {
      Matrix<T,13,7> A = randomMatrix<T,13,7>(max);
      Matrix<T,7,13> B = randomMatrix<T,7,13>(max);

      Matrix<T,13,13> AB = A * B;
      AB.transposeSelf();

      Matrix<T,13,7> BT = B.transpose();
      Matrix<T,7,13> AT = A.transpose();
      Matrix<T,13,13> BTAT = BT * AT;

      assert(AB == BTAT);
    }
  }

  cout << "  Ariphmetic" << endl << endl;
  {
    T max = static_cast<T>(3.0);

    for (int i = 0; i < 100; i++)
    {
      Matrix<T,13,7> A = randomMatrix<T,13,7>(max);
      Matrix<T,13,7> B = randomMatrix<T,13,7>(max);

      Matrix<T,13,7> AB = A + B * max;

      Matrix<T,13,7> A1 = AB - B * max;

      assert(A == A1);
    }

    for (int i = 0; i < 100; i++)
    {
      Matrix<T,4,4> A = randomMatrix<T,4,4>(max);
      Matrix<T,4,4> B = randomMatrix<T,4,4>(max);

      Matrix<T,4,4> AB = A + B * max;

      Matrix<T,4,4> A1 = AB - B * max;

      assert(A == A1);
    }

    for (int i = 0; i < 100; i++)
    {
      Matrix<T,3,3> A = randomMatrix<T,3,3>(max);
      Matrix<T,3,3> B = randomMatrix<T,3,3>(max);

      Matrix<T,3,3> AB = A + B * max;

      Matrix<T,3,3> A1 = AB - B * max;

      assert(A == A1);
    }

    for (int i = 0; i < 100; i++)
    {
      Matrix<T,2,2> A = randomMatrix<T,2,2>(max);
      Matrix<T,2,2> B = randomMatrix<T,2,2>(max);

      Matrix<T,2,2> AB = A + B * max;

      Matrix<T,2,2> A1 = AB - B * max;

      assert(A == A1);
    }
  }

  cout << "  Inversion" << endl << endl;
  {
    T max = static_cast<T>(13.0);

    Matrix<T,4,4>::setTolerance(std::numeric_limits<T>::epsilon() * static_cast<T>(100.0) * max * max);
    Matrix<T,3,3>::setTolerance(std::numeric_limits<T>::epsilon() * static_cast<T>(100.0) * max * max);

    for (int i = 0; i < 100; i++)
    {
      Matrix<T,4,4> A = randomMatrix<T,4,4>(max);
      A.setDiag(max);

      Matrix<T,4,4> A1 = +A;

      Matrix<T,4,4> I = A * A1;

      assert(I.isIdentity());
    }

    for (int i = 0; i < 100; i++)
    {
      Matrix<T,3,3> A = randomMatrix<T,3,3>(max);
      A.setDiag(max);

      Matrix<T,3,3> A1 = +A;

      Matrix<T,3,3> I = A * A1;

      assert(I.isIdentity());
    }

    for (int i = 0; i < 100; i++)
    {
      Matrix<T,2,2> A = randomMatrix<T,2,2>(max);
      A.setDiag(max);

      Matrix<T,2,2> A1 = +A;

      Matrix<T,2,2> I = A * A1;

      assert(I.isIdentity());
    }

    Matrix<T,4,4>::setDefaultTolerance();

  }

#endif

#ifdef TEST_ELEMENT

  cout << "  ===== Tests of 3D conforming finite element =====" << endl << endl;

  cout << "  Testing incorrectly defined elements (arbitrary hanging nodes), shape function value at a node must be 1, at all other nodes 0" << endl << endl;
  {
    T tolerance = static_cast<T>(0.000001);

    for (int i = 0; i < 100; i++)
    {
      std::array<LINT,Element::NUM_BASICNODES> nodes = {0,1,2,3,4,5,6,7};

      std::array<LINT,Element::NUM_HANGINGNODES> hangingNodes = {
        -1, // 8  0
        -1, // 9  1
        -1, // 10 2
        -1, // 11 3
        -1, // 12 4
        -1, // 13 5
        -1, // 14 6
        -1, // 15 7
        -1, // 16 8
        -1, // 17 9
        -1, // 18 10
        -1, // 19 11
        -1, // 20 12
        -1, // 21 13
        -1, // 22 14 
        -1, // 23 15
        -1, // 24 16
        -1  // 25 17
      };         

      // random conforming nodes
      int n = random(Element::NUM_HANGINGNODES);
      for (int j = 0; j < n; j++)
      {
        hangingNodes[random(Element::NUM_HANGINGNODES - 1)] = 1;
      }

      // create element
      Element e(nodes,hangingNodes);

      // Standard shape functions : check 1.0 is only value among all shape functions
      for (int j = 0; j < Element::NUM_BASICNODES; j++)
      {
        Tvector coords = e.nodeParms()[j];

        std::array<T,Element::NUM_NODES> func;
        e.shapeFunctions(coords,func);

        for (int k = 0; k < Element::NUM_BASICNODES; k++)
        {
          if (k == j)
          {
            if (std::abs(func[k] - T1) > tolerance)
            {
              assert(false);
            }
          } else
          {
            if (std::abs(func[k]) > tolerance)
            {
              assert(false);
            }
          }
        }
      }

      // Sum of all shape functions must be 1.0
      for (int j = 0; j < Element::NUM_NODES; j++)
      {
        Tvector coords = e.nodeParms()[j];

        std::array<T,Element::NUM_NODES> func;
        e.shapeFunctionsConforming(coords,func);

        T sum = std::accumulate(func.begin(),func.end(),T0,std::plus<T>());

        assert(std::abs(sum - T1) < tolerance);
      }

      // Conforming functions : check 1.0 
      for (int j = 0; j < Element::NUM_NODES; j++)
      {
        Tvector coords = e.nodeParms()[j];

        std::array<T,Element::NUM_NODES> func;
        e.shapeFunctionsConforming(coords,func);

        if (j < Element::NUM_BASICNODES)
        {
          for (int k = 0; k < Element::NUM_BASICNODES; k++)
          {
            if (k == j)
            {
              if (std::abs(func[k] - T1) > tolerance)
              {
                assert(false);
              }
            } else
            {
              if (std::abs(func[k]) > tolerance)
              {
                assert(false);
              }
            }
          }
          for (int k = Element::NUM_BASICNODES; k < Element::NUM_NODES; k++)
          {
            if (std::abs(func[k]) > tolerance)
            {
              assert(false);
            }
          }
        } else
        {
          if (!e.nodeExists(j))
            continue;

          // This test is incorrect : the element maybe wrong with arbitrary
          // selection of hanging nodes; i.e. all nodes on a face must be hanging
          // to produce zeroes at [0..7] and single unit in [8..25]
          //for (int k = 0; k < Element::NUM_BASICNODES; k++)
          //{
          //  if (std::abs(func[k]) > tolerance)
          //  {
          //    assert(false);
          //  }
          //}
          for (int k = Element::NUM_BASICNODES; k < Element::NUM_NODES; k++)
          {
            if (k == j)
            {
              if (std::abs(func[k] - T1) > tolerance)
              {
                assert(false);
              }
            } else
            {
              // see comment above
              //if (std::abs(func[k]) > tolerance)
              //{
              //  assert(false);
              //}
            }
          }
        }
      }
    }
  }

  cout << "  Testing correctly defined elements (whole faces with hanging nodes), shape function value at a node must be 1, at all other nodes 0" << endl << endl;
  {
    T tolerance = static_cast<T>(0.000001);

    for (int i = 0; i < 100; i++)
    {
      std::array<LINT,Element::NUM_BASICNODES> nodes = {0,1,2,3,4,5,6,7};

      std::array<LINT,Element::NUM_HANGINGNODES> hangingNodes = {
        -1, // 8  0
        -1, // 9  1
        -1, // 10 2
        -1, // 11 3
        -1, // 12 4
        -1, // 13 5
        -1, // 14 6
        -1, // 15 7
        -1, // 16 8
        -1, // 17 9
        -1, // 18 10
        -1, // 19 11
        -1, // 20 12
        -1, // 21 13
        -1, // 22 14 
        -1, // 23 15
        -1, // 24 16
        -1  // 25 17
      };         

      // random conforming nodes
      int n = random(6);
      for (int j = 0; j < n; j++)
      {
        int face = random(6);
        for (int k = 4; k < 9; k++)
          hangingNodes[Element::faces()[face][k] - 8] = 1;
      }

      // create element
      Element e(nodes,hangingNodes);

      // Standard shape functions : check 1.0 is only value among all shape functions
      for (int j = 0; j < Element::NUM_BASICNODES; j++)
      {
        Tvector coords = e.nodeParms()[j];

        std::array<T,Element::NUM_NODES> func;
        e.shapeFunctions(coords,func);

        for (int k = 0; k < Element::NUM_BASICNODES; k++)
        {
          if (k == j)
          {
            if (std::abs(func[k] - T1) > tolerance)
            {
              assert(false);
            }
          } else
          {
            if (std::abs(func[k]) > tolerance)
            {
              assert(false);
            }
          }
        }
      }

      // Sum of all shape functions must be 1.0
      for (int j = 0; j < Element::NUM_NODES; j++)
      {
        Tvector coords = e.nodeParms()[j];

        std::array<T,Element::NUM_NODES> func;
        e.shapeFunctionsConforming(coords,func);

        T sum = std::accumulate(func.begin(),func.end(),T0,std::plus<T>());

        assert(std::abs(sum - T1) < tolerance);
      }

      // Conforming functions : check 1.0 
      for (int j = 0; j < Element::NUM_NODES; j++)
      {
        Tvector coords = e.nodeParms()[j];

        std::array<T,Element::NUM_NODES> func;
        e.shapeFunctionsConforming(coords,func);

        if (j < Element::NUM_BASICNODES)
        {
          for (int k = 0; k < Element::NUM_BASICNODES; k++)
          {
            if (k == j)
            {
              if (std::abs(func[k] - T1) > tolerance)
              {
                assert(false);
              }
            } else
            {
              if (std::abs(func[k]) > tolerance)
              {
                assert(false);
              }
            }
          }
          for (int k = Element::NUM_BASICNODES; k < Element::NUM_NODES; k++)
          {
            if (std::abs(func[k]) > tolerance)
            {
              assert(false);
            }
          }
        } else
        {
          if (!e.nodeExists(j))
            continue;

          for (int k = 0; k < Element::NUM_BASICNODES; k++)
          {
            if (std::abs(func[k]) > tolerance)
            {
              assert(false);
            }
          }
          for (int k = Element::NUM_BASICNODES; k < Element::NUM_NODES; k++)
          {
            if (k == j)
            {
              if (std::abs(func[k] - T1) > tolerance)
              {
                assert(false);
              }
            } else
            {
              // see comment above
              if (std::abs(func[k]) > tolerance)
              {
                assert(false);
              }
            }
          }
        }
      }
    }
  }

  cout << "  Calculation of shape function derivatives by four methods must yield close results" << endl << endl;
  {
    T tolerance = static_cast<T>(0.0002);

    for (int i = 0; i < 10000; i++)  //!!!!!!
    {
      std::array<LINT,Element::NUM_BASICNODES> nodes = {0,1,2,3,4,5,6,7};

      std::array<LINT,Element::NUM_HANGINGNODES> hangingNodes = {
        -1, // 8  0
        -1, // 9  1
        -1, // 10 2
        -1, // 11 3
        -1, // 12 4
        -1, // 13 5
        -1, // 14 6
        -1, // 15 7
        -1, // 16 8
        -1, // 17 9
        -1, // 18 10
        -1, // 19 11
        -1, // 20 12
        -1, // 21 13
        -1, // 22 14 
        -1, // 23 15
        -1, // 24 16
        -1  // 25 17
      };         

      // random conforming nodes
      int n = random(6);
      for (int j = 0; j < n; j++)
      {
        int face = random(6);
        for (int k = 4; k < 9; k++)
          hangingNodes[Element::faces()[face][k] - 8] = 1;
      }

      // create element
      Element e(nodes,hangingNodes);

      std::array<Tvector,Element::NUM_NODES> dfunc0,dfunc1,dfunc2,dfunc3;

      // coordinates must not be close to 0.0,0.5,1.0; otherwise
      // difference between analytic and finite difference calculations
      // may be large (derivatives are discontinuous at hanging nodes due to
      // absolute value expression in shape functions (see explanation in
      // function decsriptions)  
      Tvector coords(random<T>(),random<T>(),random<T>());
      for (int k = 0; k < 3; k++)
      {
        // must be far from 0.0,0.5,1.0 at a distance larger than
        // finite difference step (0.001)
        T d = coords[k] - T(0.5);
        if (std::abs(d < T(0.01)))
        {
          if (coords[k] > T(0.5))
          {
            coords[k] = T(0.51);
          } else
          {
            coords[k] = T(0.49);
          }
        }
      }

      e.shapeDerivativesConformingExact(coords,dfunc0);
      e.shapeDerivativesConformingSimple(coords,dfunc1);
      e.shapeDerivativesNonConforming(coords,dfunc2);
      e.shapeDerivativesConformingFromSubOctant(coords,dfunc3);

      // check differences
      T diff01 = diff<T,Tvector,Element::NUM_NODES>(dfunc0,dfunc1,true);
      T diff02 = diff<T,Tvector,Element::NUM_NODES>(dfunc0,dfunc2);
      T diff03 = diff<T,Tvector,Element::NUM_NODES>(dfunc0,dfunc3);

      // compare only two close methods : exact and finite difference approximation;
      // they must be close
      if (diff01 >= 0.02)
      {
        cout << diff01 << endl;
        assert(false);
      }

      // check sum
      Tvector sum0, sum1, sum2, sum3;
      sum0 = sum1 = sum2 = sum3 = Tvector(T0,T0,T0);
      for (int i = 0; i < Element::NUM_NODES; i++)
      {
        sum0 += dfunc0[i];
        sum1 += dfunc1[i];
        sum2 += dfunc2[i];
        sum3 += dfunc3[i];
      }
      T len0 = !sum0, len1 = !sum1, len2 = !sum2, len3 = !sum3;
      assert(len0 < tolerance);
      assert(len1 < tolerance);
      assert(len2 < tolerance);
      assert(len3 < tolerance);
    }
  }

  cout << "  Calculation of values at some specific cases (correct hanging nodes on two faces)" << endl << endl;
  {
    T tolerance = static_cast<T>(0.000001);

    Tvector coords(T0,T05,T05);

    std::array<LINT,Element::NUM_BASICNODES> nodes = {0,1,2,3,4,5,6,7};

    std::array<LINT,Element::NUM_HANGINGNODES> hangingNodes = {
        -1, // 8  0
        -1, // 9  1
         1, // 10 2
         1, // 11 3
         1, // 12 4
        -1, // 13 5
         1, // 14 6
         1, // 15 7
        -1, // 16 8
        -1, // 17 9
         1, // 18 10
         1, // 19 11
        -1, // 20 12
        -1, // 21 13
        -1, // 22 14 
         1, // 23 15
         1, // 24 16
        -1  // 25 17
    };         

    // create element
    Element e(nodes,hangingNodes);

    std::array<T,Element::NUM_NODES> simplefunc,func;

    e.shapeFunctions(coords,simplefunc);
    e.shapeFunctionsConforming(coords,func);

    for (int i = 0; i < 100; i++)
    {
      std::array<T,26> values =
      {5,5,5,5,5, 5,5,5,5,5, 5,5,5,5,5, 5,5,5,5,5, 5,5,5,5,50, 5};

      values[random(26)] = random<T>() * T(10.0);

      values[24] = 50;

      T v1 = shapeFuncValue<T,T,26>(func,values,T0);
      assert(std::abs(v1 - 50) < tolerance);
    }

  }

  cout << "  Derivatives on real coordinates, element 2 x 0.5 x 0.75" << endl << endl;
  {
    Tvector min(T(0.0),T(0.0),T(0.0)), max(T(2.0),T(0.5),T(0.75)), d = max - min;

    std::vector<Tvector> totalCoordinates;
    addElementCoordinates<T,Tvector>(min,max,totalCoordinates);

    Element e({0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25});

    Tvector parms(0.75,0.75,0.75);
    std::array<Tvector,Element::NUM_NODES> derivatives;
    e.shapeDerivativesConformingExact(parms,derivatives);

    Matrix<T,3,Element::NUM_NODES> B;
    T jacobian;
    bool ok = e.Bmatrix(totalCoordinates,derivatives,B,jacobian);

    // coordinates of this element nodes
    Matrix<T,Element::NUM_NODES,3> coords = e.coordMatrix(totalCoordinates);

    Matrix<T,3,Element::NUM_NODES> mderivatives = e.shapeFunctionDerivativesMatrix(derivatives).transpose();

    // derivatives on parameters
    Matrix<T,3,3> derivOnParm = mderivatives * coords;
    // derivatives on real coordinates
    Matrix<T,3,3> derivOnXYZ = B * coords;

    // convert to 3-component vector
    Tvector dd;
    for (int j = 0; j < 3; j++)
      dd.XYZ[j] = derivOnParm[j][j];

    assert(derivOnXYZ.isIdentity());
    assert(dd == d);
    T volume = dd.XYZ[0] * dd.XYZ[1] * dd.XYZ[2];
    assert(std::abs(jacobian - volume) < dd.tolerance());
  }

  cout << "  Derivatives on parameters and real coordinates, random tests :" << endl;
  cout << "(1) element derivatives on X,Y,Z must be units (identity matrix)" << endl;
  cout << "(2) element derivatives on parameters must be equal to box sizes" << endl;
  cout << "(3) element volume must be equal to Jacobi matrix determinant" << endl << endl;
  {
    Matrix<T,3,3>::setTolerance(T(0.0001));
//    Matrix<T,3,3>::setTolerance(std::numeric_limits<T>::epsilon() * static_cast<T>(10000.0));

    for (int i = 0; i < 100; i++)
    {
      Tvector min(random<T>(),random<T>(),random<T>()), 
        max(min.X + T(0.01) + random<T>(),min.Y + T(0.01) + random<T>(),min.Z + T(0.01) + random<T>()), 
        d = max - min;

      std::vector<Tvector> totalCoordinates;
      addElementCoordinates<T,Tvector>(min,max,totalCoordinates);

      std::array<LINT,Element::NUM_NODES> nodes = 
        {0,1,2,3,4,5,6,7,-8,-9,-10,-11,-12,-13,-14,-15,-16,-17,-18,-19,-20,-21,-22,-23,-24,-25};

      // random conforming nodes
      int n = random(6);
      for (int j = 0; j < n; j++)
      {
        int face = random(6);
        for (int k = 4; k < 9; k++)
          nodes[Element::faces()[face][k]] = std::abs(nodes[Element::faces()[face][k]]);
      }

      // create element
      Element e(nodes);

      // random parameter
      Tvector parms(random<T>(),random<T>(),random<T>());

      std::array<Tvector,Element::NUM_NODES> derivatives;
      e.shapeDerivativesConformingExact(parms,derivatives);

      Matrix<T,3,Element::NUM_NODES> B;
      T jacobian;
      bool ok = e.Bmatrix(totalCoordinates,derivatives,B,jacobian);

      if (ok)
      {
        // coordinates of this element nodes
        Matrix<T,Element::NUM_NODES,3> coords = e.coordMatrix(totalCoordinates);

        Matrix<T,3,Element::NUM_NODES> mderivatives = e.shapeFunctionDerivativesMatrix(derivatives).transpose();

        // derivatives on parameters
        Matrix<T,3,3> derivOnParm = mderivatives * coords;
        // derivatives on real coordinates
        Matrix<T,3,3> derivOnXYZ = B * coords;

        // convert to 3-component vector
        Tvector dd;
        for (int j = 0; j < 3; j++)
          dd.XYZ[j] = derivOnParm[j][j];

        assert(derivOnXYZ.isIdentity());
        assert(dd == d);
        T volume = dd.XYZ[0] * dd.XYZ[1] * dd.XYZ[2];
        assert(std::abs(jacobian - volume) < dd.tolerance());
      }
    }
  }

#endif

#ifdef TEST_SOLUTION

  cout << "  ===== Solution to banded system with SIMD and multiple threads  =====" << endl << endl;

  {
                                  // this is system order
    #define N 5000
                                  // this is half band width including diagonal
    #define MB 500
                                  // number of threads
    #define NUM_THREADS 1

    using namespace std::chrono;
    high_resolution_clock::time_point t1,t2,t3,t4;

    progressprint = true;
															    // create matrices
    Allocator allocator;
    FEMVirtMatrix<T> matrix(&allocator,N,MB);
    int NN = ((N + 4) / 4) * 4;

	  std::vector<T> bvirt(NN,T0);
	  std::vector<T> b(NN,T0);
	  std::vector<T> rvirt(NN,T0);
	  std::vector<T> bsimple(NN,T0);
															// fill matrices
	  for (int i = 0; i < N; i++)
	  {
		  for (int j = i - (MB - 1); j <= i + (MB - 1); j++)
		  {
			  if ((j >= 0) && (j < N))
			  {
				  if (i == j)
				  {
					  *matrix.getElement(i,j) = 4.0;
				  } else
				  {
					  *matrix.getElement(i,j) = (j % 2) ? T(0.0) : T(1.0);
				  }
			  }
		  }

		  bvirt[i] = bsimple[i] = b[i] = T1;
	  }

                              // store matrix
    matrix.storeMatrix();
                              // we solve the system with a usual C solver, which, 
                              // however, takes into account banded structure 
                              // of the system matrix
    cout << "Solving system of size " << N << " by " << MB << endl;
    cout << "1. solving with simple C solver..." << endl;

                              // solve system
	  t1 = high_resolution_clock::now();
	  bool ok1 = matrix.solveSystemSimple(&bsimple[0],T(0.0000001));
	  t2 = high_resolution_clock::now();
                              // and now solve the system with the parallel solver
    printf("2. solving with parallel solver...\n");
                              // restore matrix
    matrix.restoreMatrix(false);
                              // solve
	  t3 = high_resolution_clock::now();
	  bool ok2 = matrix.solveSystem(&bvirt[0],T(0.0000001),NUM_THREADS);
	  t4 = high_resolution_clock::now();
                              // report time
    duration<double> time_span1 = duration_cast<duration<double>>(t2 - t1);
    duration<double> time_span2 = duration_cast<duration<double>>(t4 - t3);

    cout << "Time : simple solver " << time_span1.count() << " sec, SIMD parallel solver " << 
      time_span2.count() << ", acceleration " << time_span1.count() / time_span2.count() << endl;

    cout << "Tesing results" << endl;
                              // compute resuduals
    std::vector<T> resvirt(NN,T0);
    std::vector<T> ressimple(NN,T0);

    matrix.restoreMatrix(false);
    matrix.multiply(&bvirt[0],&resvirt[0],NUM_THREADS);
    vectorsSubtract<T>(&resvirt[0],&b[0],NN);

    matrix.restoreMatrix(false);
    matrix.multiply(&bsimple[0],&ressimple[0],NUM_THREADS);
    vectorsSubtract<T>(&ressimple[0],&b[0],NN);

                              // compare RHS vectors
    printf("Comparing solutions...\n");
    T maxdiff = T0;
    T maxv = T0;
    T minv = T0;
    T diff;
    T maxresvirt = T0;
    T maxressimple = T0;
    for (size_t i = 0; i < N; i++)
    {
      diff = std::abs(bvirt[i] - bsimple[i]);
//!!! cout << i << " " << diff << " " << bvirt[i] << " " << bsimple[i] << endl;
      if (diff > maxdiff) maxdiff = diff;
      if (i == 0 || bvirt[i] > maxv) maxv = bvirt[i];
      if (i == 0 || bvirt[i] < minv) minv = bvirt[i];

      if (std::abs(resvirt[i]) > maxresvirt)
        maxresvirt = std::abs(resvirt[i]);
      if (std::abs(ressimple[i]) > maxressimple)
        maxressimple = std::abs(ressimple[i]);
    };
    cout << "Max difference b/w solutions is " << maxdiff << ", max solution value " << maxv <<
      ", min value " << minv << ", residual " << maxressimple << ", parallel residual " << maxresvirt << endl;

    #undef N
    #undef MB
    #undef NUM_THREADS
  }
  cout << endl;

#endif

  cout << "  ===== Solving Laplace equation with conforming finite elements =====" << endl << endl;

  {
    // See explanation at the top

    #define NUM_THREADS 4

/*
         (3)-----------------------(12)---------(21)--------(30) Z=1.0
         /|                        /|           /|          /|
        / |                       / |          / |         / |
       /  |                     (11)---------(20)--------(29)|
      /   |                     /|  |        /|  |       /|  |
     /    |                    / | (9)------/-|-(18)--- /-| (27)
   (2)-----------------------(10)---------(19)--------(28)| /|
    |     |                   |  |/ |      |  |/ |     |  |/ |
    |     |                   | (8)------- |-(17)----- | (26)|
    |     |                   | /|  |      | /|  |     | /|  |
    |    (1)------------------|/-|-(6)-----|/---(15)---|/-|-(24) Y=1.0
    |    /                   (7)----------(16)--------(25)| /
    |   /                     |  |/        |  |/       |  |/
    |  /                      | (5)--------|-(14)------|-(23)
  Z ^ /Y                      | /          | /         | /
    |/                        |/           |/          |/
   (0)->---------------------(4)----------(13)--------(22)
       X                     X=1.0                    X=2.0
*/

    std::vector<Tvector> totalCoordinates = {
      Tvector(T(0.0),T(0.0),T(0.0)),
      Tvector(T(0.0),T(1.0),T(0.0)),
      Tvector(T(0.0),T(0.0),T(1.0)),
      Tvector(T(0.0),T(1.0),T(1.0)),  // 3

      Tvector(T(1.0),T(0.0),T(0.0)),
      Tvector(T(1.0),T(0.5),T(0.0)),
      Tvector(T(1.0),T(1.0),T(0.0)),
      Tvector(T(1.0),T(0.0),T(0.5)),
      Tvector(T(1.0),T(0.5),T(0.5)),
      Tvector(T(1.0),T(1.0),T(0.5)),
      Tvector(T(1.0),T(0.0),T(1.0)),
      Tvector(T(1.0),T(0.5),T(1.0)),
      Tvector(T(1.0),T(1.0),T(1.0)), // 12

      Tvector(T(1.5),T(0.0),T(0.0)),
      Tvector(T(1.5),T(0.5),T(0.0)),
      Tvector(T(1.5),T(1.0),T(0.0)),
      Tvector(T(1.5),T(0.0),T(0.5)),
      Tvector(T(1.5),T(0.5),T(0.5)),
      Tvector(T(1.5),T(1.0),T(0.5)),
      Tvector(T(1.5),T(0.0),T(1.0)),
      Tvector(T(1.5),T(0.5),T(1.0)),
      Tvector(T(1.5),T(1.0),T(1.0)), // 21

      Tvector(T(2.0),T(0.0),T(0.0)),
      Tvector(T(2.0),T(0.5),T(0.0)),
      Tvector(T(2.0),T(1.0),T(0.0)),
      Tvector(T(2.0),T(0.0),T(0.5)),
      Tvector(T(2.0),T(0.5),T(0.5)),
      Tvector(T(2.0),T(1.0),T(0.5)),
      Tvector(T(2.0),T(0.0),T(1.0)),
      Tvector(T(2.0),T(0.5),T(1.0)),
      Tvector(T(2.0),T(1.0),T(1.0))  // 30
    };

    std::vector<Conforming3D<T,Tvector>> elements = {
      Conforming3D<T,Tvector>({0,4,6,1,2,10,12,3, -1,5,-1,-1, -1,7,9,-1, -1,11,-1,-1, -1,-1,8,-1,-1,-1}),

      Conforming3D<T,Tvector>({4,13,14,5,7,16,17,8,     -1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1}),
      Conforming3D<T,Tvector>({5,14,15,6,8,17,18,9,     -1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1}),
      Conforming3D<T,Tvector>({7,16,17,8,10,19,20,11,   -1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1}),
      Conforming3D<T,Tvector>({8,17,18,9,11,20,21,12,   -1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1}),

      Conforming3D<T,Tvector>({13,22,23,14,16,25,26,17, -1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1}),
      Conforming3D<T,Tvector>({14,23,24,15,17,26,27,18, -1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1}),
      Conforming3D<T,Tvector>({16,25,26,17,19,28,29,20, -1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1}),
      Conforming3D<T,Tvector>({17,26,27,18,20,29,30,21, -1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1})
    };

    // calculate bandwidth
    size_t bandwidth = 0;
    for (auto e : elements)
    {
      LINT min,max;
      e.nodeMinMax(min,max);
      LINT diff = max - min;
      assert(diff > 0);

      bandwidth = std::max<LINT>(diff,bandwidth);
    }
    bandwidth++;

    // create FEM matrix
    Allocator allocator;
    FEMVirtMatrix<T> matrix(&allocator,totalCoordinates.size(),bandwidth);

    // fill matrix
    for (auto e : elements)
    {
      Matrix<T,Element::NUM_NODES,Element::NUM_NODES> ematrix;
      e.stiffnessMatrix(totalCoordinates,ematrix,GAUSSINT_2);

      for (int i = 0; i < Element::NUM_NODES; i++)
      {
        if (e.nodes()[i] < 0)
          continue;
        for (int j = 0; j < Element::NUM_NODES; j++)
        {
          if (e.nodes()[j] < 0)
            continue;
          
          matrix.addToElement(e.nodes()[i],e.nodes()[j],ematrix[i][j]);
        }
      }
    }

    // array size must have length of 4 or 2 floats to svvoid buffer override with SIMD
    size_t N = ((totalCoordinates.size() + 4) / 4) * 4;
    assert(N >= totalCoordinates.size());
  
    // right-hand side of the system
    std::vector<T> b(N,T0); // added for SIMD to avoid buffer override

    // apply boundary conditions at X = 0 - zero potential; b[0]..1,2,3 is already zero
    matrix.degenerateEquation(0);
    matrix.degenerateEquation(1);
    matrix.degenerateEquation(2);
    matrix.degenerateEquation(3);

    // apply boundary conditions at X = 1 - unit potential;
    for (size_t i = 22; i <= 30; i++)
    {
      matrix.degenerateEquation(i);
      b[i] = T1;
    }

    // store matrix
    matrix.storeMatrix();

    // to test simple solution
    auto rhs = b;

/*
         (3)-----------------------(12)---------(21)--------(30) Z=1.0
         /|                        /|           /|          /|
        / |                       / |          / |         / |
       /  |                     (11)---------(20)--------(29)|
      /   |                     /|  |        /|  |       /|  |
     /    |                    / | (9)------/-|-(18)--- /-| (27)
   (2)-----------------------(10)---------(19)--------(28)| /|
    |     |                   |  |/ |      |  |/ |     |  |/ |
    |     |                   | (8)------- |-(17)----- | (26)|
    |     |                   | /|  |      | /|  |     | /|  |
    |    (1)------------------|/-|-(6)-----|/---(15)---|/-|-(24) Y=1.0
    |    /                   (7)----------(16)--------(25)| /
    |   /                     |  |/        |  |/       |  |/
    |  /                      | (5)--------|-(14)------|-(23)
  Z ^ /Y                      | /          | /         | /
    |/                        |/           |/          |/
   (0)->---------------------(4)----------(13)--------(22)
       X                     X=1.0                    X=2.0
*/

    //bool ok1 = matrix.solveSystemSimple(&bsimple[0],T(0.0000001));

    //matrix.restoreMatrix(false);
    bool ok2 = matrix.solveSystem(&b[0],T(0.0000001),NUM_THREADS);

    // compute residual
    std::vector<T> residual(N,T0);
    matrix.restoreMatrix(false);
    matrix.multiply(&b[0],&residual[0],NUM_THREADS);
    vectorsSubtract<T>(&residual[0],&rhs[0],N);
    T maxres = T0;
    for (auto v : residual)
    {
      maxres = std::max(maxres,v);
    }

    cout << "Max residual " << maxres << endl;

    // compute value of potential and its gradient at centres of elements (they can be computed everywhere)
    std::vector<Tvector> coordinates;
    std::vector<T> potential;
    std::vector<Tvector> gradients;

     cout << "\tX\tY\tZ\tpotential\tgradient" << endl << endl;

    for (auto e : elements)
    {
      Tvector parms(T05,T05,T05);

      Tvector coord = e.coordinate(parms,totalCoordinates);
      coordinates.push_back(coord);

      T pot = e.value(parms,b);
      potential.push_back(pot);

      Tvector grad = e.gradient(parms,totalCoordinates,b);
      gradients.push_back(grad);

      cout << "\t" << coord.X << "\t" << coord.Y << "\t" << coord.Z << "\t" << pot << "\t\t" << grad.X << "\t" << 
        grad.Y << "\t" << grad.Z << endl; 
    }
    
  }

  cout << endl << "Press ENTER to exit..." << endl;
  getchar();
}
