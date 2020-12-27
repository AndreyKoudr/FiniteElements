/**
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

#pragma once

#include <limits>
#include <algorithm>
#include <functional>
#include <numeric> 
#include <array>
#include <string>
#include <vector>
#include <assert.h>
#include "Types.h"
#include "Matrix.h"
#include "GaussInt.h"

/**
  Conforming linear finite element
  --------------------------------

  K.-W. Cheng, T.-P. Fries. XFEM with Hanging Nodes for Two-phase 
Incompressible Flow.

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
    Due to e.g. T c0a = T1 - std::abs(c0); in shapeFunctionsConforming() - abs value 
  function expressions, derivatives at and near hanging nodes with a parameter 0.5
  are not well defined as abs() is discontinuous, therefore the two functions 
  shapeDerivativesConformingSimple() based on finite differences and  shapeDerivativesExact() - 
  exact expression MAY GIVE VERY DIFFERENT results as the derivative does not exist at x = 0 
  for function y = |x|.
    In addition, first derivatives for a linear approximation must be constant, so the use
  of shapeDerivativesConformingFromSubOctant() which calculates derivatives at centres 
  of sub-octants looks most consistent.
*/

/**
  8 sub-elements :

            7-------------------6 
           /|  (6)         (7) /|    
          / |                 / |     
         /  |                /  |          
        4-------------------5   |    
        |   (4)         (5) |   |     
        |   |               |   |    
        |   |               |   |
        |   3---------------|---2 
        |  /  (2)         (3)  / 
        | /                 | /
        |/ (0)         (1)  |/
        0-------------------1
                                          
*/

static const std::array<std::array<int,8>,8> subElements_ =
{{
  {0,8,20,11,12,21,26,24},
  {8,1,9,20,21,13,22,26},
  {11,20,10,3,24,26,23,15},
  {20,9,2,10,26,22,14,23},
  
  {12,21,26,24,4,16,25,19},
  {21,13,22,26,16,5,17,25},
  {24,26,23,15,19,25,18,7},
  {26,22,14,23,25,17,6,18}
}};

/**
  Element faces :
 
        7-------------------6 
       /|                  /|    
      / |      (5)        / |     
     /  |                /  |          
    4-----------(3)-----5   |    
    |   |               |   |     
    |(0)|               |(1)|    
    |   |               |   |
    |   3----(2)--------|---2 
    |  /                |  / 
    | /        (4)      | /
    |/                  |/
    0-------------------1

  Hex faces nodal numbers, mind the order : 
 -e direction (-x)
 +e direction (+x)
 -n direction (-y)
 +n direction (+y)
 -s direction (-z)
 +s direction (+z)
 First 4 nodes are regular, the other 5 are hanging, last one being face central
                                          
*/

static const std::array<std::array<int,9>,6>  elementFaces_ =
{{
  {3,0,4,7,   11,12,19,15,  24},
  {1,2,6,5,   9,14,17,13,   22},
  {0,1,5,4,   8,13,16,12,   21},
  {2,3,7,6,   10,15,18,14,  23},
  {0,3,2,1,   11,10,9,8,    20},
  {7,4,5,6,   19,16,17,18,  25}
}};

/**
  Every face is divided into 4 subfaces [6][4][4]
*/
static const std::array<std::array<std::array<int,4>,4>,6> subFaceNodes_ =
{{
  {{{20,8,0,11},  {20,11,3,10}, {20,10,2,9},  {20,9,1,8}}},
  {{{21,12,0,8},  {21,8,1,13},  {21,13,5,16}, {21,16,4,12}}}, 
  {{{22,13,1,9},  {22,9,2,14},  {22,14,6,17}, {22,17,5,13}}},
  {{{23,14,2,10}, {23,10,3,15}, {23,15,7,18}, {23,18,6,14}}},
  {{{24,15,3,11}, {24,11,0,12}, {24,12,4,19}, {24,19,7,15}}},
  {{{25,16,5,17}, {25,17,6,18}, {25,18,7,19}, {25,19,4,16}}}
}};

/**
  Element edges, the third number is hanging node [12][3] :

        7--------(6)--------6 
       /|                  /|    
     (7)|                (5)|     
     /  |                /  |          
    4---------(4)-------5  (10)  
    |  (11)             |   |     
    |   |               |   |    
    |   |               |   |
   (8)  3-------(2)----(9)--2 
    |  /                |  / 
    |(3)                |(1) 
    |/                  |/
    0--------(0)--------1
                                          
*/

static const std::array<std::array<int,3>,12> elementEdges_ =
{{
  {0,1,8},
  {1,2,9},
  {2,3,10},
  {3,0,11},

  {4,5,16},
  {5,6,17},
  {6,7,18},
  {7,4,19},

  {0,4,12},
  {1,5,13},
  {2,6,14},
  {3,7,15}
}};

/**
  Two faces for every edge
*/
static const std::array<std::array<int,2>,12> edgeFaces_ =
{{
  {2,4},
  {1,4},
  {3,4},
  {0,4},

  {2,5},
  {1,5},
  {3,5},
  {0,5},

  {0,2},
  {2,1},
  {1,3},
  {3,0}
}};

template <typename T, typename Tvector> class Conforming3D {
public:
  /** Total number of nodes */
  static const size_t NUM_NODES = 26;

  /** Number of basic nodes which always exist */
  static const size_t NUM_BASICNODES = 8;

  /** Number of hanging nodes */
  static const size_t NUM_HANGINGNODES = NUM_NODES - NUM_BASICNODES;

  /** Constructor */
  Conforming3D() = delete;

  /** Constructors */
  //Conforming3D(const std::array<LINT,NUM_BASICNODES> &nodes);
  Conforming3D(const std::array<LINT,NUM_NODES> &nodes);
  Conforming3D(const std::array<LINT,NUM_BASICNODES> &nodes, const std::array<LINT,NUM_HANGINGNODES> &hanging);

  /** Copy constructor */
  Conforming3D(const Conforming3D& copy);

  /** Assignment operator */
  Conforming3D& operator=(const Conforming3D& copy);

  /** Nodes */
  std::array<LINT,NUM_NODES>& nodes()
  {
    return nodes_;
  }

  /** Is this element with hanging nodes? */
  bool hasHangingNodes();

  /** Get min and max node numbers */
  void  nodeMinMax(LINT& min, LINT& max);

  /** Compute linear shape functions, local coordinates (parms) being [0..1] */
  static void shapeFunctions(const Tvector &parms, std::array<T,NUM_NODES> &func);

  /** Compute linear conforming shape functions, local coordinates (parms) [0..1] */
  void shapeFunctionsConforming(const Tvector &parms, std::array<T,NUM_NODES> &func);

  /** Compute derivatives of linear conforming shape functions, local coordinates (parms) [0..1],
  non-conforming 8 nodes only used. 
  */
  void shapeDerivativesNonConforming(const Tvector &parms, std::array<Tvector,NUM_NODES> &func);

  /** Compute derivatives of linear conforming shape functions, local coordinates (parms) [0..1],
    calculation is made by calling shapeFunctionsConforming() to compute finite differences;
    result is approximate. 

      Very important!
      Due to e.g. T c0a = T1 - std::abs(c0); in shapeFunctionsConforming() - abs value 
    function expressions, derivatives at and near hanging nodes with a parameter 0.5
    are not well defined as abs() is discontinuous, therefore the two functions 
    shapeDerivativesConformingSimple() based on finite differences and  shapeDerivativesExact() - 
    exact expression MAY GIVE VERY DIFFERENT results as the derivative does not exist at x = 0 
    for function y = |x|.
      In addition, first derivatives for a linear approximation must be constant, so the use
    of shapeDerivativesConformingFromSubOctant() which calculates derivatives at centres 
    of sub-octants looks most consistent.
  */
  void shapeDerivativesConformingSimple(const Tvector &parms, std::array<Tvector,NUM_NODES> &func);

  /** Compute derivatives of linear conforming shape functions, local coordinates (parms) [0..1],
    calculation is exact.  

      Very important!
      Due to e.g. T c0a = T1 - std::abs(c0); in shapeFunctionsConforming() - abs value 
    function expressions, derivatives at and near hanging nodes with a parameter 0.5
    are not well defined as abs() is discontinuous, therefore the two functions 
    shapeDerivativesConformingSimple() based on finite differences and  shapeDerivativesExact() - 
    exact expression MAY GIVE VERY DIFFERENT results as the derivative does not exist at x = 0 
    for function y = |x|.
      In addition, first derivatives for a linear approximation must be constant, so the use
    of shapeDerivativesConformingFromSubOctant() which calculates derivatives at centres 
    of sub-octants looks most consistent.
  */
  void shapeDerivativesConformingExact(const Tvector &parms, std::array<Tvector,NUM_NODES> &func);

  /** Compute derivatives of linear conforming shape functions, local coordinates 
    (parms) [0..1], parms position is moved to a centre of closest sub-octant and 
    calculation is axact at this point */
  void shapeDerivativesConformingFromSubOctant(const Tvector &parms, std::array<Tvector,NUM_NODES> &func);

  /** Get parameters of all nodes */
  static const std::array<Tvector,NUM_NODES + 1>& nodeParms() 
  {
    return nodeParms_;
  }

  /** Nodes for 8 subelements (icluding central number 26) */
  static const std::array<std::array<int,8>,8>& subElements();

  /** Element 6 faces, see the picture for elemenFaces_ */
  static const std::array<std::array<int,9>,6>& faces();

  /** Every face is divided into 4 subfaces [6][4][4] */
  static const std::array<std::array<std::array<int,4>,4>,6>& subFaceNodes();

  /** Element edges, the third number is hanging node [12][3] */
  static const std::array<std::array<int,3>,12>& edges();

  /** Two faces for every edge [12][2] */
  static const std::array<std::array<int,2>,12>& edgeFaces();

  /** Node exists? */
  bool nodeExists(const size_t i);

  /** Return last error string */
  static std::string errorString()
  {
    return errorString_;
  }


  // Finite element matrices

  /** Get 26 x 3 node coordinate matrix, they are zero for non-existing nodes */
  Matrix<T,NUM_NODES,3> coordMatrix(const std::vector<Tvector> &totalCoordinates)
  {
    Matrix<T,NUM_NODES,3> coords;

    for (int i = 0; i < NUM_NODES; i++)
    {
      if (nodes_[i] >= 0)
      {
        coords[i][0] = totalCoordinates[nodes_[i]].XYZ[0];
        coords[i][1] = totalCoordinates[nodes_[i]].XYZ[1];
        coords[i][2] = totalCoordinates[nodes_[i]].XYZ[2];
      } else
      {
        coords[i][0] = coords[i][1] = coords[i][2] = T(0.0);
      }
    }

    return coords;
  }

  /** Get 26 x 1 shape function matrix, zeroes for non-existing nodes */
  static Matrix<T,NUM_NODES,1> shapeFunctionMatrix(const std::array<T,NUM_NODES> &func)
  {
    Matrix<T,NUM_NODES,1> funcmatrix;
    funcmatrix.fill(T0);

    for (int i = 0; i < NUM_NODES; i++)
    {
      funcmatrix[i][0] = func[i];
    }

    return funcmatrix;
  }

  /** Get 26 x 3 shape function derivative matrix, they are zero for non-existing nodes */
  static Matrix<T,NUM_NODES,3> shapeFunctionDerivativesMatrix(const std::array<Tvector,NUM_NODES> &func)
  {
    Matrix<T,NUM_NODES,3> dermatrix;
    dermatrix.fill(T0);

    for (int i = 0; i < NUM_NODES; i++)
    {
      dermatrix[i][0] = func[i].XYZ[0];
      dermatrix[i][1] = func[i].XYZ[1];
      dermatrix[i][2] = func[i].XYZ[2];
    }

    return dermatrix;
  }


  /**
    Calculate matrix B to get deriveatives on X,Y,Z (real coordinates) by
  {Ux,Uy,Uz} = B[3,NUM_NODES] * { U nodal values [NUM_NODES] }
    See J.J.Connor, C.A.Brebbia Finite Element Techniques for Fluid Flow, page 125, equation (i)
  */
  bool Bmatrix(const std::vector<Tvector> &totalCoordinates, 
    const std::array<Tvector,NUM_NODES> &derivatives, Matrix<T,3,NUM_NODES> &B, T &jacobian);

  /** Calculate element stiffness matrix [K] */
  bool  stiffnessMatrix(const std::vector<Tvector>& totalCoordinates,
    Matrix<T,Conforming3D<T,Tvector>::NUM_NODES,Conforming3D<T,Tvector>::NUM_NODES>& matrix, const int power);

  /** Calculate element mass matrix [M] */
  bool  massMatrix(const std::vector<Tvector>& totalCoordinates,
    Matrix<T,Conforming3D<T,Tvector>::NUM_NODES,Conforming3D<T,Tvector>::NUM_NODES>& matrix, const int power);


  // ===== Solution vector =====

  /** Get nodal values for this element */
  Matrix<T,1,Conforming3D<T,Tvector>::NUM_NODES> values(const std::vector<T>& totalValues);

  /** Calculate coordinate at a point defined by parameter */
  Tvector coordinate(const Tvector& parms, const std::vector<Tvector>& totalCoordinates);

  /** Calculate a value at a point defined by parameter */
  T value(const Tvector& parms, const std::vector<T>& totalValues);

  /** Calculate value gradient at a point defined by parameter */
  Tvector gradient(const Tvector& parms, const std::vector<Tvector> &totalCoordinates,
    const std::vector<T>& totalValues);

private:
  /** Node numbers; for hanging nodes 8..25 a node does not exist if -1. */
  std::array<LINT,NUM_NODES> nodes_;

  /** Parametric values for all nodes plus central node */
  static const std::array<Tvector,27> nodeParms_;

  /** Last error string */
  static std::string errorString_;

  /** Compute derivatives of linear shape functions, local coordinates (parms) being [0..1] */
  void shapeDerivatives(const Tvector &parms, std::array<Tvector,NUM_NODES> &func);

};





template <typename  T, typename Tvector>
  std::string Conforming3D<T,Tvector>::errorString_;

template <typename  T, typename Tvector> Conforming3D<T,Tvector>::
  Conforming3D(const std::array<LINT,NUM_NODES> &nodes)
{
  nodes_ = nodes;
}

template <typename  T, typename Tvector> Conforming3D<T,Tvector>::
  Conforming3D(const std::array<LINT,NUM_BASICNODES> &nodes, const std::array<LINT,NUM_HANGINGNODES> &hanging)
{
  std::copy(nodes.begin(),nodes.end(),nodes_.begin());
  std::copy(hanging.begin(),hanging.end(),nodes_.begin() + NUM_BASICNODES);
}

template <typename  T, typename Tvector> Conforming3D<T,Tvector>::
  Conforming3D(const Conforming3D& copy)
{
  nodes_ = copy.nodes_;
}

template <typename  T, typename Tvector> Conforming3D<T,Tvector>& Conforming3D<T,Tvector>::
  operator = (const Conforming3D& copy)
{
  nodes_ = copy.nodes_;

  return *this;
}

template <typename  T, typename Tvector> bool Conforming3D<T,Tvector>::hasHangingNodes()
{
  return (std::find_if(nodes_.begin() + NUM_BASICNODES,nodes_.end(),[](auto v) { return (v >= 0); })
    != nodes_.end());
}

template <typename  T, typename Tvector> bool Conforming3D<T,Tvector>::nodeExists(const size_t i)
{
  return (nodes_[i] >= 0);
}

template <typename  T, typename Tvector> static const std::array<std::array<int,8>,8>& 
  Conforming3D<T,Tvector>::subElements()
{
  return subElements_;
}

template <typename  T, typename Tvector> const std::array<std::array<int,9>,6>& 
  Conforming3D<T,Tvector>::faces()
{
  return elementFaces_;
}

template <typename  T, typename Tvector> const std::array<std::array<std::array<int,4>,4>,6>& 
  Conforming3D<T,Tvector>::subFaceNodes()
{
  return subFaceNodes_;
}

template <typename  T, typename Tvector> const std::array<std::array<int,3>,12>& 
  Conforming3D<T,Tvector>::edges()
{
  return elementEdges_;
}

template <typename  T, typename Tvector> const std::array<std::array<int,2>,12>& 
  Conforming3D<T,Tvector>::edgeFaces()
{
  return edgeFaces_;
}

template <typename  T, typename Tvector> void Conforming3D<T,Tvector>::shapeFunctions(
  const Tvector &parms, std::array<T,NUM_NODES> &func)
{
#ifndef NDEBUG
  T tol = TOLERANCE(T);
  assert(parms.XYZ[0] >= T0 - tol && parms.XYZ[0] <= T1 + tol);
  assert(parms.XYZ[1] >= T0 - tol && parms.XYZ[1] <= T1 + tol);
  assert(parms.XYZ[2] >= T0 - tol && parms.XYZ[2] <= T1 + tol);
#endif

  std::fill(func.begin(),func.end(),T0);

  T c0 = parms.XYZ[0] * T2 - T1;
  T c1 = parms.XYZ[1] * T2 - T1;
  T c2 = parms.XYZ[2] * T2 - T1;

  T c0m = T1 - c0;
  T c0p = T1 + c0;
  T c1m = T1 - c1;
  T c1p = T1 + c1;
  T c2m = T1 - c2;
  T c2p = T1 + c2;

  func[0] = T0125 * c0m * c1m * c2m;
  func[1] = T0125 * c0p * c1m * c2m;
  func[2] = T0125 * c0p * c1p * c2m;
  func[3] = T0125 * c0m * c1p * c2m;
  func[4] = T0125 * c0m * c1m * c2p;
  func[5] = T0125 * c0p * c1m * c2p;
  func[6] = T0125 * c0p * c1p * c2p;
  func[7] = T0125 * c0m * c1p * c2p;
}

template <class T, class Tvector> void Conforming3D<T,Tvector>::shapeDerivatives(
  const Tvector &parms, std::array<Tvector,NUM_NODES> &func)
{
#ifndef NDEBUG
  T tol = TOLERANCE(T);
  assert(parms.XYZ[0] >= T0 - tol && parms.XYZ[0] <= T1 + tol);
  assert(parms.XYZ[1] >= T0 - tol && parms.XYZ[1] <= T1 + tol);
  assert(parms.XYZ[2] >= T0 - tol && parms.XYZ[2] <= T1 + tol);
#endif

  std::fill(func.begin(),func.end(),Tvector(0,0,0));

  T c0 = parms.XYZ[0] * T2 - T1;
  T c1 = parms.XYZ[1] * T2 - T1;
  T c2 = parms.XYZ[2] * T2 - T1;

  T c0m = T1 - c0;
  T c0p = T1 + c0;
  T c1m = T1 - c1;
  T c1p = T1 + c1;
  T c2m = T1 - c2;
  T c2p = T1 + c2;

  func[0] = Tvector(
    T0125 * (-T1) * c1m * c2m,
    T0125 * c0m * (-T1) * c2m,
    T0125 * c0m * c1m * (-T1)
  );

  func[1] = Tvector( 
    T0125 * (+T1) * c1m * c2m,
    T0125 * c0p * (-T1) * c2m,
    T0125 * c0p * c1m * (-T1)
  );

  func[2] = Tvector( 
    T0125 * (+T1) * c1p * c2m,
    T0125 * c0p * (+T1) * c2m,
    T0125 * c0p * c1p * (-T1)
  );

  func[3] = Tvector( 
    T0125 * (-T1) * c1p * c2m,
    T0125 * c0m * (+T1) * c2m,
    T0125 * c0m * c1p * (-T1)
  );

  func[4] = Tvector( 
    T0125 * (-T1) * c1m * c2p,
    T0125 * c0m * (-T1) * c2p,
    T0125 * c0m * c1m * (+T1)
  );

  func[5] = Tvector( 
    T0125 * (+T1) * c1m * c2p,
    T0125 * c0p * (-T1) * c2p,
    T0125 * c0p * c1m * (+T1)
  );

  func[6] = Tvector( 
    T0125 * (+T1) * c1p * c2p,
    T0125 * c0p * (+T1) * c2p,
    T0125 * c0p * c1p * (+T1)
  );

  func[7] = Tvector( 
    T0125 * (-T1) * c1p * c2p,
    T0125 * c0m * (+T1) * c2p,
    T0125 * c0m * c1p * (+T1)
  );
}

template <class T, class Tvector> void Conforming3D<T,Tvector>::shapeFunctionsConforming(
  const Tvector &parms, std::array<T,NUM_NODES> &func)
{
  // main 8-node shape functions
  shapeFunctions(parms,func);

  if (hasHangingNodes())
  {
    T c0 = parms.XYZ[0] * T2 - T1;
    T c1 = parms.XYZ[1] * T2 - T1;
    T c2 = parms.XYZ[2] * T2 - T1;

    T c0a = T1 - std::abs(c0);
    T c1a = T1 - std::abs(c1);
    T c2a = T1 - std::abs(c2);

    T c0m = T1 - c0;
    T c0p = T1 + c0;
    T c1m = T1 - c1;
    T c1p = T1 + c1;
    T c2m = T1 - c2;
    T c2p = T1 + c2;

    // these node numbers are offset by 8 : func[0] is actually func[8]

    // faces
    T func20 = (nodes_[NUM_BASICNODES + 12] < 0) ? T0 : (T05 * c0a * c1a * c2m);
    T func25 = (nodes_[NUM_BASICNODES + 17] < 0) ? T0 : (T05 * c0a * c1a * c2p);

    T func21 = (nodes_[NUM_BASICNODES + 13] < 0) ? T0 : (T05 * c0a * c1m * c2a);
    T func23 = (nodes_[NUM_BASICNODES + 15] < 0) ? T0 : (T05 * c0a * c1p * c2a);

    T func24 = (nodes_[NUM_BASICNODES + 16] < 0) ? T0 : (T05 * c0m * c1a * c2a);
    T func22 = (nodes_[NUM_BASICNODES + 14] < 0) ? T0 : (T05 * c0p * c1a * c2a);

    // edges (preliminary shape functions)
    T func8p  = (nodes_[NUM_BASICNODES + 0] < 0) ? T0 :  (T025 * c0a * c1m * c2m);
    T func10p = (nodes_[NUM_BASICNODES + 2] < 0) ? T0 :  (T025 * c0a * c1p * c2m);
    T func16p = (nodes_[NUM_BASICNODES + 8] < 0) ? T0 :  (T025 * c0a * c1m * c2p);
    T func18p = (nodes_[NUM_BASICNODES + 10] < 0) ? T0 : (T025 * c0a * c1p * c2p);

    T func9p  = (nodes_[NUM_BASICNODES + 1] < 0) ? T0 :  (T025 * c0p * c1a * c2m);
    T func11p = (nodes_[NUM_BASICNODES + 3] < 0) ? T0 :  (T025 * c0m * c1a * c2m);
    T func17p = (nodes_[NUM_BASICNODES + 9] < 0) ? T0 :  (T025 * c0p * c1a * c2p);
    T func19p = (nodes_[NUM_BASICNODES + 11] < 0) ? T0 : (T025 * c0m * c1a * c2p);

    T func14p = (nodes_[NUM_BASICNODES + 6] < 0) ? T0 :  (T025 * c0p * c1p * c2a);
    T func15p = (nodes_[NUM_BASICNODES + 7] < 0) ? T0 :  (T025 * c0m * c1p * c2a);
    T func13p = (nodes_[NUM_BASICNODES + 5] < 0) ? T0 :  (T025 * c0p * c1m * c2a);
    T func12p = (nodes_[NUM_BASICNODES + 4] < 0) ? T0 :  (T025 * c0m * c1m * c2a);

    T func8 = func8p - T05 * (func20 + func21);
    T func9 = func9p - T05 * (func20 + func22);
    T func10 = func10p - T05 * (func20 + func23);
    T func11 = func11p - T05 * (func20 + func24);
    T func12 = func12p - T05 * (func21 + func24);
    T func13 = func13p - T05 * (func21 + func22);
    T func14 = func14p - T05 * (func22 + func23);
    T func15 = func15p - T05 * (func23 + func24);
    T func16 = func16p - T05 * (func25 + func21);
    T func17 = func17p - T05 * (func25 + func22);
    T func18 = func18p - T05 * (func25 + func23);
    T func19 = func19p - T05 * (func25 + func24);

    func[0] += (-T05 * (func8p + func11p + func12p) + T025 * (func20 + func21 + func24));
    func[1] += (-T05 * (func9p + func8p + func13p) + T025 * (func20 + func22 + func21));
    func[2] += (-T05 * (func10p + func9p + func14p) + T025 * (func20 + func23 + func22));
    func[3] += (-T05 * (func11p + func10p + func15p) + T025 * (func20 + func24 + func23));
    func[4] += (-T05 * (func12p + func16p + func19p) + T025 * (func25 + func24 + func21));
    func[5] += (-T05 * (func13p + func17p + func16p) + T025 * (func25 + func21 + func22));
    func[6] += (-T05 * (func14p + func18p + func17p) + T025 * (func25 + func22 + func23));
    func[7] += (-T05 * (func15p + func19p + func18p) + T025 * (func25 + func23 + func24));

    // additional 18 conforming shape functions for hanging nodes
    // (indices are offset by 8 -> [0] means [8])
    func[NUM_BASICNODES + 0] = func8;
    func[NUM_BASICNODES + 1] = func9;
    func[NUM_BASICNODES + 2] = func10;
    func[NUM_BASICNODES + 3] = func11;
    func[NUM_BASICNODES + 4] = func12;
    func[NUM_BASICNODES + 5] = func13;
    func[NUM_BASICNODES + 6] = func14;
    func[NUM_BASICNODES + 7] = func15;
    func[NUM_BASICNODES + 8] = func16;
    func[NUM_BASICNODES + 9] = func17;
    func[NUM_BASICNODES + 10] = func18;
    func[NUM_BASICNODES + 11] = func19;
    func[NUM_BASICNODES + 12] = func20;
    func[NUM_BASICNODES + 13] = func21;
    func[NUM_BASICNODES + 14] = func22;
    func[NUM_BASICNODES + 15] = func23;
    func[NUM_BASICNODES + 16] = func24;
    func[NUM_BASICNODES + 17] = func25;
  }
}

template <class T, class Tvector> void Conforming3D<T,Tvector>::shapeDerivativesConformingSimple(
  const Tvector &parms, std::array<Tvector,NUM_NODES> &func)
{
  // Use standard shape functions
  if (!hasHangingNodes())
  {
    shapeDerivatives(parms,func);

    // Derivatives are calculated for -1..+1, we need 0..1
    std::transform(func.begin(),func.begin() + NUM_BASICNODES,func.begin(),[](auto v) { return v * T2; } );
  } else
  {
#ifndef NDEBUG
    for (int j = 0; j < 3; j++)
    {
      T d = std::abs(parms.XYZ[j] - T(0.5));
      if (d < 0.005) // five times the finite-difference step
      {
        assert(false && "shapeFunctionsConforming() called near 0.5 cross over discontinuity caused by y = |x|");
      }
    }
#endif

    T step = T(0.001);
    Tvector coords0 = parms - Tvector(step,step,step);
    Tvector coords1 = parms + Tvector(step,step,step);

    Tvector dcoords;
    for (int j = 0; j < 3; j++)
    {
      if (coords0.XYZ[j] < 0)
        coords0.XYZ[j] = 0;
      if (coords1.XYZ[j] > 1)
        coords1.XYZ[j] = 1;
      dcoords.XYZ[j] = coords1.XYZ[j] - coords0.XYZ[j];
      assert(dcoords[j] > 0);
    }

    for (int j = 0; j < 3; j++)
    {
      std::array<T,NUM_NODES> func0;
      Tvector c0 = parms;
      c0.XYZ[j] = coords0[j];

      shapeFunctionsConforming(c0,func0);

      std::array<T,NUM_NODES> func1;
      Tvector c1 = parms;
      c1.XYZ[j] = coords1[j];
        
      shapeFunctionsConforming(c1,func1);

      for (int k = 0; k < NUM_NODES; k++)
      {
        func[k].XYZ[j] = (func1[k] - func0[k]) / dcoords[j];
      }
    }
  } 
}

template <class T, class Tvector> void Conforming3D<T,Tvector>::shapeDerivativesConformingExact(
  const Tvector &parms, std::array<Tvector,NUM_NODES> &func)
{
  // Use standard shape functions
  if (!hasHangingNodes())
  {
    shapeDerivatives(parms,func);

    // derivatives are calculated for -1..+1, we need 0..1
    std::transform(func.begin(),func.begin() + NUM_BASICNODES,func.begin(),[](auto v) { return v * T2; } );
  } else
  {
  #ifndef NDEBUG
    T tol = TOLERANCE(T);
    assert(parms.XYZ[0] >= T0 - tol && parms.XYZ[0] <= T1 + tol);
    assert(parms.XYZ[1] >= T0 - tol && parms.XYZ[1] <= T1 + tol);
    assert(parms.XYZ[2] >= T0 - tol && parms.XYZ[2] <= T1 + tol);
  #endif

    T c0 = parms.XYZ[0] * T2 - T1;
    T c1 = parms.XYZ[1] * T2 - T1;
    T c2 = parms.XYZ[2] * T2 - T1;

    T c0a = T1 - std::abs(c0);
    T c1a = T1 - std::abs(c1);
    T c2a = T1 - std::abs(c2);

    T c0m = T1 - c0;
    T c0p = T1 + c0;
    T c1m = T1 - c1;
    T c1p = T1 + c1;
    T c2m = T1 - c2;
    T c2p = T1 + c2;

    // these node numbers are offset by 8 : func[0] is actually func[8]

    // faces
    Tvector func20 = (nodes_[NUM_BASICNODES + 12] < 0) ? Tvector(0,0,0) : 
      Tvector(
        // T05 * c0a * c1a * c2m
        (c0 >= 0) ? (T05 * (-T1) * c1a * c2m) : (T05 * (+T1) * c1a * c2m),
        (c1 >= 0) ? (T05 * c0a * (-T1) * c2m) : (T05 * c0a * (+T1) * c2m),
        T05 * c0a * c1a * (-T1)
      );
    Tvector func25 = (nodes_[NUM_BASICNODES + 17] < 0) ? Tvector(0,0,0) :
      Tvector(
        // T05 * c0a * c1a * c2p
        (c0 >= 0) ? (T05 * (-T1) * c1a * c2p) : (T05 * (+T1) * c1a * c2p),
        (c1 >= 0) ? (T05 * c0a * (-T1) * c2p) : (T05 * c0a * (+T1) * c2p),
        T05 * c0a * c1a * (+T1)
      );
    Tvector func21 = (nodes_[NUM_BASICNODES + 13] < 0) ? Tvector(0,0,0) : 
      Tvector(
        // T05 * c0a * c1m * c2a
        (c0 >= 0) ? (T05 * (-T1) * c1m * c2a) : (T05 * (+T1) * c1m * c2a),
        T05 * c0a * (-T1) * c2a,
        (c2 >= 0) ? (T05 * c0a * c1m * (-T1)) : (T05 * c0a * c1m * (+T1))
      );
    Tvector func23 = (nodes_[NUM_BASICNODES + 15] < 0) ? Tvector(0,0,0) : 
      Tvector(
        // T05 * c0a * c1p * c2a
        (c0 >= 0) ? (T05 * (-T1) * c1p * c2a) : (T05 * (+T1) * c1p * c2a),
        T05 * c0a * (+T1) * c2a,
        (c2 >= 0) ? (T05 * c0a * c1p * (-T1)) : (T05 * c0a * c1p * (+T1))
      );

    Tvector func24 = (nodes_[NUM_BASICNODES + 16] < 0) ? Tvector(0,0,0) : 
      Tvector(
        // T05 * c0m * c1a * c2a
        T05 * (-T1) * c1a * c2a,
        (c1 >= 0) ? (T05 * c0m * (-T1) * c2a) : (T05 * c0m * (+T1) * c2a),
        (c2 >= 0) ? (T05 * c0m * c1a * (-T1)) : (T05 * c0m * c1a * (+T1))
      );

    Tvector func22 = (nodes_[NUM_BASICNODES + 14] < 0) ? Tvector(0,0,0) : 
      Tvector(
        // T05 * c0p * c1a * c2a
        T05 * (+T1) * c1a * c2a,
        (c1 >= 0) ? (T05 * c0p * (-T1) * c2a) : (T05 * c0p * (+T1) * c2a),
        (c2 >= 0) ? (T05 * c0p * c1a * (-T1)) : (T05 * c0p * c1a * (+T1))
      );

    // edges (preliminary shape functions)
    Tvector func8p  = (nodes_[NUM_BASICNODES + 0] < 0) ? Tvector(0,0,0) :
      Tvector(
        // T025 * c0a * c1m * c2m
        (c0 >= 0) ? (T025 * (-T1) * c1m * c2m) : (T025 * (+T1) * c1m * c2m),
        T025 * c0a * (-T1) * c2m,
        T025 * c0a * c1m * (-T1)
      );

    Tvector func10p = (nodes_[NUM_BASICNODES + 2] < 0) ? Tvector(0,0,0) :
      Tvector(
        // T025 * c0a * c1p * c2m
        (c0 >= 0) ? (T025 * (-T1) * c1p * c2m) : (T025 * (+T1) * c1p * c2m),
        T025 * c0a * (+T1) * c2m,
        T025 * c0a * c1p * (-T1)
      );

    Tvector func16p = (nodes_[NUM_BASICNODES + 8] < 0) ? Tvector(0,0,0) :
      Tvector(
        // T025 * c0a * c1m * c2p
        (c0 >= 0) ? (T025 * (-T1) * c1m * c2p) : (T025 * (+T1) * c1m * c2p),
        T025 * c0a * (-T1) * c2p,
        T025 * c0a * c1m * (+T1)
      );

    Tvector func18p = (nodes_[NUM_BASICNODES + 10] < 0) ? Tvector(0,0,0) :
      Tvector(
        // T025 * c0a * c1p * c2p
        (c0 >= 0) ? (T025 * (-T1) * c1p * c2p) : (T025 * (+T1) * c1p * c2p),
        T025 * c0a * (+T1) * c2p,
        T025 * c0a * c1p * (+T1)
      );

    Tvector func9p = (nodes_[NUM_BASICNODES + 1] < 0) ? Tvector(0,0,0) :
      Tvector(
        // T025 * c0p * c1a * c2m
        T025 * (+T1) * c1a * c2m,
        (c1 >= 0) ? (T025 * c0p * (-T1) * c2m) : (T025 * c0p * (+T1) * c2m),
        T025 * c0p * c1a * (-T1)
      );

    Tvector func11p = (nodes_[NUM_BASICNODES + 3] < 0) ? Tvector(0,0,0) :
      Tvector(
        // T025 * c0m * c1a * c2m
        T025 * (-T1) * c1a * c2m,
        (c1 >= 0) ? (T025 * c0m * (-T1) * c2m) : (T025 * c0m * (+T1) * c2m),
        T025 * c0m * c1a * (-T1)
      );

    Tvector func17p = (nodes_[NUM_BASICNODES + 9] < 0) ? Tvector(0,0,0) :
      Tvector(
        // T025 * c0p * c1a * c2p
        T025 * (+T1) * c1a * c2p,
        (c1 >= 0) ? (T025 * c0p * (-T1) * c2p) : (T025 * c0p * (+T1) * c2p),
        T025 * c0p * c1a * (+T1)
      );

    Tvector func19p = (nodes_[NUM_BASICNODES + 11] < 0) ? Tvector(0,0,0) :
      Tvector(
        // T025 * c0m * c1a * c2p
        T025 * (-T1) * c1a * c2p,
        (c1 >= 0) ? (T025 * c0m * (-T1) * c2p) : (T025 * c0m * (+T1) * c2p),
        T025 * c0m * c1a * (+T1)
      );

    Tvector func14p = (nodes_[NUM_BASICNODES + 6] < 0) ? Tvector(0,0,0) :
      Tvector(
        // T025 * c0p * c1p * c2a
        T025 * (+T1) * c1p * c2a,
        T025 * c0p * (+T1) * c2a,
        (c2 >= 0) ? (T025 * c0p * c1p * (-T1)) : (T025 * c0p * c1p * (+T1))
      );

    Tvector func15p = (nodes_[NUM_BASICNODES + 7] < 0) ? Tvector(0,0,0) :
      Tvector(
        // T025 * c0m * c1p * c2a
        T025 * (-T1) * c1p * c2a,
        T025 * c0m * (+T1) * c2a,
        (c2 >= 0) ? (T025 * c0m * c1p * (-T1)) : (T025 * c0m * c1p * (+T1))
      );

    Tvector func13p = (nodes_[NUM_BASICNODES + 5] < 0) ? Tvector(0,0,0) :
      Tvector(
        // T025 * c0p * c1m * c2a
        T025 * (+T1) * c1m * c2a,
        T025 * c0p * (-T1) * c2a,
        (c2 >= 0) ? (T025 * c0p * c1m * (-T1)) : (T025 * c0p * c1m * (+T1))
      );

    Tvector func12p = (nodes_[NUM_BASICNODES + 4] < 0) ? Tvector(0,0,0) :
      Tvector(
        // T025 * c0m * c1m * c2a
        T025 * (-T1) * c1m * c2a,
        T025 * c0m * (-T1) * c2a,
        (c2 >= 0) ? (T025 * c0m * c1m * (-T1)) : (T025 * c0m * c1m * (+T1))
      );

    Tvector func8 = func8p - (func20 + func21) * T05;
    Tvector func9 = func9p - (func20 + func22) * T05;
    Tvector func10 = func10p - (func20 + func23) * T05;
    Tvector func11 = func11p - (func20 + func24) * T05;
    Tvector func12 = func12p - (func21 + func24) * T05;
    Tvector func13 = func13p - (func21 + func22) * T05;
    Tvector func14 = func14p - (func22 + func23) * T05;
    Tvector func15 = func15p - (func23 + func24) * T05;
    Tvector func16 = func16p - (func25 + func21) * T05;
    Tvector func17 = func17p - (func25 + func22) * T05;
    Tvector func18 = func18p - (func25 + func23) * T05;
    Tvector func19 = func19p - (func25 + func24) * T05;

    // main 8-node shape functions
    shapeDerivatives(parms,func);

    Tvector df0 = (func8p + func11p + func12p) * (-T05) + (func20 + func21 + func24) * T025;
    Tvector df1 = (func9p + func8p + func13p) * (-T05) + (func20 + func22 + func21) * T025;
    Tvector df2 = (func10p + func9p + func14p) * (-T05) + (func20 + func23 + func22) * T025;
    Tvector df3 = (func11p + func10p + func15p) * (-T05) + (func20 + func24 + func23) * T025;
    Tvector df4 = (func12p + func16p + func19p) * (-T05) + (func25 + func24 + func21) * T025;
    Tvector df5 = (func13p + func17p + func16p) * (-T05) + (func25 + func21 + func22) * T025;
    Tvector df6 = (func14p + func18p + func17p) * (-T05) + (func25 + func22 + func23) * T025;
    Tvector df7 = (func15p + func19p + func18p) * (-T05) + (func25 + func23 + func24) * T025;

    func[0] += df0;
    func[1] += df1;
    func[2] += df2;
    func[3] += df3;
    func[4] += df4;
    func[5] += df5;
    func[6] += df6;
    func[7] += df7;

    // additional 18 conforming shape functions for hanging nodes
    // (indices offset by 8 -> [0] means [8])
    func[NUM_BASICNODES + 0] = func8;
    func[NUM_BASICNODES + 1] = func9;
    func[NUM_BASICNODES + 2] = func10;
    func[NUM_BASICNODES + 3] = func11;
    func[NUM_BASICNODES + 4] = func12;
    func[NUM_BASICNODES + 5] = func13;
    func[NUM_BASICNODES + 6] = func14;
    func[NUM_BASICNODES + 7] = func15;
    func[NUM_BASICNODES + 8] = func16;
    func[NUM_BASICNODES + 9] = func17;
    func[NUM_BASICNODES + 10] = func18;
    func[NUM_BASICNODES + 11] = func19;
    func[NUM_BASICNODES + 12] = func20;
    func[NUM_BASICNODES + 13] = func21;
    func[NUM_BASICNODES + 14] = func22;
    func[NUM_BASICNODES + 15] = func23;
    func[NUM_BASICNODES + 16] = func24;
    func[NUM_BASICNODES + 17] = func25;

    // derivatives are calculated for -1..+1, we need 0..1
    std::transform(func.begin(),func.end(),func.begin(),[](auto v) { return v * T2; } );
  }
}

template <class T, class Tvector> void Conforming3D<T,Tvector>::shapeDerivativesNonConforming(
  const Tvector &parms, std::array<Tvector,NUM_NODES> &func)
{
  shapeDerivatives(parms,func);

  // derivatives are calculated for -1..+1, we need 0..1
  std::transform(func.begin(),func.end(),func.begin(),[](auto v) { return v * T2; } );
}

template <class T, class Tvector> void Conforming3D<T,Tvector>::shapeDerivativesConformingFromSubOctant(
  const Tvector &parms, std::array<Tvector,NUM_NODES> &func)
{
  // Use standard shape functions
  if (!hasHangingNodes())
  {
    shapeDerivatives(parms,func);

    // derivatives are calculated for -1..+1, we need 0..1
    std::transform(func.begin(),func.begin() + NUM_BASICNODES,func.begin(),[](auto v) { return v * T2; } );
  } else
  {
    // Calculate derivatives at centres of 8 sub-elements
    Tvector c;
    for (int j = 0; j < 3; j++)
    {
      int sub = static_cast<int>(parms.XYZ[j] * T2);
      if (sub < 0) sub = 0;
      if (sub > 1) sub = 1;

      c.XYZ[j] = static_cast<T>(sub) * T05 + T025;
    }

    shapeDerivativesConformingExact(c,func);
  }
}

template <typename T, typename Tvector> const std::array<Tvector,27> Conforming3D<T,Tvector>::nodeParms_ =
{
  Tvector(0,0,0),
  Tvector(1,0,0),
  Tvector(1,1,0),
  Tvector(0,1,0),

  Tvector(0,0,1),
  Tvector(1,0,1),
  Tvector(1,1,1),
  Tvector(0,1,1),

  Tvector(0.5,0,0),
  Tvector(1,0.5,0),
  Tvector(0.5,1,0),
  Tvector(0,0.5,0),

  Tvector(0,0,0.5),
  Tvector(1,0,0.5),
  Tvector(1,1,0.5),
  Tvector(0,1,0.5),

  Tvector(0.5,0,1),
  Tvector(1,0.5,1),
  Tvector(0.5,1,1),
  Tvector(0,0.5,1),

  Tvector(0.5,0.5,0),
  Tvector(0.5,0,0.5),
  Tvector(1,0.5,0.5),
  Tvector(0.5,1,0.5),
  Tvector(0,0.5,0.5),
  Tvector(0.5,0.5,1),

  Tvector(0.5,0.5,0.5)
};

// Multiply shape functions by values at nodes and sum up
template <typename T, typename Tvalue, size_t N> Tvalue shapeFuncValue(const std::array<T,N> &func,
  const std::array<Tvalue,N> &values, Tvalue zero)
{
  std::array<Tvalue,N> temp;

  std::transform(values.begin(),values.end(),func.begin(),temp.begin(),[](auto a, auto b) { return a * b; } );

  Tvalue sum = std::accumulate(temp.begin(),temp.end(),zero,[](auto a, auto b) { return a + b; } );

  return sum;
}

template <class T, class Tvector> bool Conforming3D<T,Tvector>::
  Bmatrix(const std::vector<Tvector> &totalCoordinates, const std::array<Tvector,NUM_NODES> &derivatives,  
  Matrix<T,3,NUM_NODES> &B, T &jacobian)
{
  // coordinates of this element nodes
  Matrix<T,NUM_NODES,3> coords = coordMatrix(totalCoordinates);

  // matrix for derivatives of shape functions
  Matrix<T,NUM_NODES,3> func = shapeFunctionDerivativesMatrix(derivatives);

  // Jacobi matrix
  Matrix<T,3,NUM_NODES> coordsT = coords.transpose();
  Matrix<T,3,3> Jacobi = coordsT * func;
  // invert matrix
  Matrix<T,3,3> J1 = +Jacobi;

  // check Jacobian value
  jacobian = Jacobi.lastDeterminant();
  if (std::abs(jacobian) < Jacobi.tolerance())
  {
    errorString_ = "Jacobian is close to zero when calculating derivatives on real coordinates";
    return false;
  }

  // fill matrix B[26 x 3] = J-1 * B^ (Connor & Brebbia, Rus edition, page 115, equation 9)
  // or J.J.Connor, C.A.Brebbia Finite Element Techniques for Fluid Flow, page 125, equation (i)
  Matrix<T,3,NUM_NODES> funcT = func.transpose();
  B = J1 * funcT;

  return true;
}

template <class T, class Tvector> bool Conforming3D<T,Tvector>::stiffnessMatrix(
  const std::vector<Tvector>& totalCoordinates,
  Matrix<T,Conforming3D<T,Tvector>::NUM_NODES,Conforming3D<T,Tvector>::NUM_NODES>& matrix, const int power)
{
  // zero matrix
  matrix.fill(T(0.0));

  // shape function derivatives
  for (int i = 0; i < GaussInt[power].numpoints; i++)
  {
    T e = static_cast<T>(GaussInt[power].knots[i]);
    T ew = static_cast<T>(GaussInt[power].weights[i]);
    for (int j = 0; j < GaussInt[power].numpoints; j++)
    {
      T n = static_cast<T>(GaussInt[power].knots[j]);
      T nw = static_cast<T>(GaussInt[power].weights[j]);
      for (int k = 0; k < GaussInt[power].numpoints; k++)
      {
        T c = static_cast<T>(GaussInt[power].knots[k]);
        T cw = static_cast<T>(GaussInt[power].weights[k]);

        Tvector parms((e + T1) * T05,(n + T1) * T05,(c + T1) * T05);
        std::array<Tvector,NUM_NODES> funcderivatives;
#if 1
        shapeDerivativesConformingFromSubOctant(parms,funcderivatives);
#else
        shapeDerivativesConformingExact(parms,funcderivatives);
        // use GAUSSINT_4 for integration over discontinuities
#endif
        Matrix<T,3,NUM_NODES> B;
        T jacobian;
        bool ok = Bmatrix(totalCoordinates,funcderivatives,B,jacobian);
        if (!ok)
          return false;

        for (int l = 0; l < NUM_NODES; l++)
        {
          for (int m = 0; m < NUM_NODES; m++)
          {
            // add to integral
            T value = (B[0][l] * B[0][m] + B[1][l] * B[1][m] + B[2][l] * B[2][m]) * 
              std::abs(jacobian) * ew * nw * cw; 
            matrix[l][m] += value;
          }
        }
      }
    }
  }

  return true;
}

template <class T, class Tvector> bool Conforming3D<T,Tvector>::massMatrix(
  const std::vector<Tvector>& totalCoordinates,
  Matrix<T,Conforming3D<T,Tvector>::NUM_NODES,Conforming3D<T,Tvector>::NUM_NODES>& matrix, const int power)
{
  // zero matrix
  matrix.fill(T(0.0));

  // shape function derivatives
  for (int i = 0; i < GaussInt[power].numpoints; i++)
  {
    T e = static_cast<T>(GaussInt[power].knots[i]);
    T ew = static_cast<T>(GaussInt[power].weights[i]);
    for (int j = 0; j < GaussInt[power].numpoints; j++)
    {
      T n = static_cast<T>(GaussInt[power].knots[j]);
      T nw = static_cast<T>(GaussInt[power].weights[j]);
      for (int k = 0; k < GaussInt[power].numpoints; k++)
      {
        T c = static_cast<T>(GaussInt[power].knots[k]);
        T cw = static_cast<T>(GaussInt[power].weights[k]);

        Tvector parms((e + T1) * T05,(n + T1) * T05,(c + T1) * T05);
        std::array<T,NUM_NODES> func;
        shapeFunctiionsConforming(parms,func);
        std::array<Tvector,NUM_NODES> funcderivatives;
#if 1
        shapeDerivativesConformingFromSubOctant(parms,funcderivatives);
#else
        shapeDerivativesConformingExact(parms,funcderivatives);
        // use GAUSSINT_4 for integration over discontinuities
#endif

        Matrix<T,3,NUM_NODES> B;
        T jacobian;
        bool ok = Bmatrix(totalCoordinates,funcderivatives,B,jacobian);
        if (!ok)
          return false;

        for (int l = 0; l < NUM_NODES; l++)
        {
          for (int m = 0; m < NUM_NODES; m++)
          {
            // add to integral
            T value = func[l] * func[m] * std::abs(jacobian) * ew * nw * cw;
            matrix[l][m] += value;
          }
        }
      }
    }
  }

  return true;
}

template <class T, class Tvector> void Conforming3D<T,Tvector>::nodeMinMax(LINT& min, LINT& max)
{
  min = +std::numeric_limits<LINT>::max();
  max = -std::numeric_limits<LINT>::max();

  for (size_t j = 0; j < NUM_NODES; j++)
  {
    if (nodes_[j] >= 0)
    {
      min = std::min(min,nodes_[j]);
      max = std::max(max,nodes_[j]);
    }
  }
}

template <class T, class Tvector> Matrix<T,1,Conforming3D<T,Tvector>::NUM_NODES> 
  Conforming3D<T,Tvector>::values(const std::vector<T>& totalValues)
{
  Matrix<T,1,Conforming3D<T,Tvector>::NUM_NODES> values;
  values.fill(T0);

  for (size_t i = 0; i < Conforming3D<T,Tvector>::NUM_NODES; i++)
  {
    if (nodes_[i] >= 0)
    {
      values[0][i] = totalValues[nodes_[i]];
    }
  }

  return values;
}

template <class T, class Tvector> Tvector Conforming3D<T,Tvector>::coordinate(const Tvector& parms,
  const std::vector<Tvector>& totalCoordinates)
{
  // nodal coordinate values
  Matrix<T,Conforming3D<T,Tvector>::NUM_NODES,3> coords = coordMatrix(totalCoordinates);

  // shape functions
  std::array<T,NUM_NODES> func;
  shapeFunctionsConforming(parms,func);
  Matrix<T,Conforming3D<T,Tvector>::NUM_NODES,1> funcmatrix = shapeFunctionMatrix(func);

  // multiply
  Matrix<T,3,1> res = coords.transpose() * funcmatrix;

  return Tvector(res[0][0],res[1][0],res[2][0]);
}

template <class T, class Tvector> T Conforming3D<T,Tvector>::value(const Tvector& parms,
  const std::vector<T>& totalValues)
{
  // nodal values
  Matrix<T,1,Conforming3D<T,Tvector>::NUM_NODES> nvalues = values(totalValues);

  // shape functions
  std::array<T,NUM_NODES> func;
  shapeFunctionsConforming(parms,func);
  Matrix<T,Conforming3D<T,Tvector>::NUM_NODES,1> funcmatrix = shapeFunctionMatrix(func);

  // multiply
  Matrix<T,1,1> res = nvalues * funcmatrix;

  return res[0][0];
}

template <class T, class Tvector> Tvector Conforming3D<T,Tvector>::gradient(const Tvector& parms,
  const std::vector<Tvector> &totalCoordinates, const std::vector<T>& totalValues)
{
  // nodal values
  Matrix<T,1,Conforming3D<T,Tvector>::NUM_NODES> nvalues = values(totalValues);

  std::array<Tvector,NUM_NODES> funcderivatives;
#if 1
  shapeDerivativesConformingFromSubOctant(parms,funcderivatives);
#else
  shapeDerivativesConformingExact(parms,funcderivatives);
#endif

  Matrix<T,3,NUM_NODES> B;
  T jacobian;
  bool ok = Bmatrix(totalCoordinates,funcderivatives,B,jacobian);

  Matrix<T,3,1> d = B * nvalues.transpose();

  return Tvector(d[0][0],d[1][0],d[2][0]);
}
