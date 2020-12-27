# FiniteElements

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

  - T - class for real numbers
  - Tvector - class for 3D vector with up to 4 components

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
These formulae are 100% reliable. The class does not contain any transforms.

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

  This is one big conforming finite element with five hanging nodes (5,7,8,9,11) and eight 
non-conforming small finite elements each with 8 nodes.
  The "non-conforming" plane is X = 1.0 to achieve continuous approximation
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
  - Ready 100% reliable code for conforming finite element. Remember that derivatives of shape functions are
discontinuous within the element.
  - Matrix class for small matrices with 100% reliable inversion formulae. If you run into a problem, 
you might have forgotten to transpose a matrix somewhere. Matrix operations are slow, no special attention 
to code speed.
  - Vector class for 4-component vectors Vector4 is accelerated by SIMD, but it cannot help much without
reorganisation of the whole code.
  - FEMVirtMatrix is SIMD accelerated with normal speedup 4 times or more, but multithreading does not help 
very much, maybe vecause of my computer is with just a single hardware set of XMM registers.
  - The last test - solution to a Laplace equation is very very simple and produces correct results. We need
a mesh generation code for proper testing (to come).

  References
  ----------
  - J.J.Connor, C.A.Brebbia. Finite Element Techniques for Fluid Flow (about finite elements).
  - K.-W. Cheng, T.-P. Fries. XFEM with Hanging Nodes for Two-phase Incompressible Flow (shape 
functions for conforming finite element). 
