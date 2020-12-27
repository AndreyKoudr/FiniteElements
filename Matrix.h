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

#pragma once

#include <algorithm>
#include <string>

/** Matrix with N rows and M columns */
template <typename T, size_t N, size_t M> class Matrix : public std::array<std::array<T,M>,N>
{
public:
  // Constructors
  Matrix() = default;
  Matrix<T,N,M>(const std::array<T,M> &data);

  // Equal with tolerance
  bool operator == (const Matrix<T,N,M> &other) const;

  // Fill matrix with value
  void fill(const T value);

  // Transpose
  Matrix<T,M,N> transpose() const;

  // Transpose in place (if N == M)
  void transposeSelf();

  // +
  Matrix<T,N,M> operator + (const Matrix<T,N,M> &other) const;
  void operator += (const Matrix<T,N,M> &other);

  // -
  Matrix<T,N,M> operator - (const Matrix<T,N,M> &other) const;
  void operator -= (const Matrix<T,N,M> &other);

  // * scalar
  Matrix<T,N,M> operator * (const T mult) const;
  void operator *= (const T mult);

  // determinant, available for 2 x 2, 3 x 3 and 4 x 4
  T determinant() const;

  // Inverse matrix, 2 x 2, 3 x 3 and 4 x 4 only
  Matrix<T,N,M> operator + ();

  // Is it identity matrix?
  bool isIdentity() const;

  // Set identity matrix (only for square matrix)
  void setIdentity();

  // Set equal diagonal elements
  void setDiag(const T value);

  // Tolerance
  static T tolerance()
  {
    return tolerance_;
  }

  // Set tolerance
  static void setTolerance(const T tolerance)
  {
    tolerance_ = tolerance;
  }

  // Set default tolerance
  static void setDefaultTolerance()
  {
    tolerance_ = std::numeric_limits<T>::epsilon() * static_cast<T>(100.0);
  }

  // Last determinant calculated during matrix inversion
  T lastDeterminant() const;

  // Return last error string
  static std::string errorString()
  {
    return errorString_;
  }

private:
  // determinant after last inversion
  T determinant_ = T(1.0);

  // tolerance
  static T tolerance_;

  // last error string
  static std::string errorString_;
};

template <typename T, size_t N, size_t M>
  T Matrix<T,N,M>::tolerance_ = std::numeric_limits<T>::epsilon() * static_cast<T>(100.0);

template <typename T, size_t N, size_t M>
  std::string Matrix<T,N,M>::errorString_;

template <typename T, size_t N, size_t M>  Matrix<T,N,M>::Matrix(const std::array<T,M> &data)
{
  (*this)[0] = data;
}

template <typename T, size_t N, size_t M> bool Matrix<T,N,M>::operator == (const Matrix<T,N,M> &other) const
{
  for (int i = 0; i < N; i++)
  {
    for (int j = 0; j < M; j++)
    {
      if (std::abs((*this)[i][j] - other[i][j]) > tolerance_)
        return false;
    }
  }

  return true;
}

template <typename T, size_t N, size_t M> void Matrix<T,N,M>::fill(const T value)
{
  std::for_each(this->begin(),this->end(),
    [&](auto &v) { std::fill(v.begin(),v.end(),value); } );
}

template <typename T, size_t N, size_t M> Matrix<T,M,N> Matrix<T,N,M>::transpose() const
{
  Matrix<T,M,N> transposed;

  for (size_t i = 0; i < N; i++)
  {
    for (size_t j = 0; j < M; j++)
    {
      transposed[j][i] = (*this)[i][j];
    }
  }

  return transposed;
}

template <typename T, size_t N, size_t M> void Matrix<T,N,M>::transposeSelf()
{
  static_assert(N == M,"Self-transpose requires square matrix");

  Matrix<T,M,N> transposed = transpose();

  *this = transposed;
}

template <typename T, size_t I, size_t K, size_t J> Matrix<T,I,J> operator * 
  (const Matrix<T,I,K> &m0, const Matrix<T,K,J> &m1)
{
  Matrix<T,I,J> result;

  for (int i = 0; i < I; i++)
  {
    for (int j = 0; j < J; j++)
    {
      result[i][j] = static_cast<T>(0.0);
      for (int k = 0; k < K; k++)
      { 
        result[i][j] += m0[i][k] * m1[k][j];
      }
    }
  }

  return result;
}

template <typename T, size_t N, size_t M> Matrix<T,N,M> Matrix<T,N,M>::operator + 
  (const Matrix<T,N,M> &other) const
{
  Matrix<T,N,M> result;

  for (int i = 0; i < N; i++)
  {
    for (int j = 0; j < M; j++)
    {
      result[i][j] = (*this)[i][j] + other[i][j];
    }
  }

  return result;
}

template <typename T, size_t N, size_t M> void Matrix<T,N,M>::operator += (const Matrix<T,N,M> &other)
{
  *this = *this + other;
}

template <typename T, size_t N, size_t M> Matrix<T,N,M> Matrix<T,N,M>::operator - 
  (const Matrix<T,N,M> &other) const
{
  Matrix<T,N,M> result;

  for (int i = 0; i < N; i++)
  {
    for (int j = 0; j < M; j++)
    {
      result[i][j] = (*this)[i][j] - other[i][j];
    }
  }

  return result;
}

template <typename T, size_t N, size_t M> void Matrix<T,N,M>::operator -= (const Matrix<T,N,M> &other)
{
  *this = *this - other;
}

template <typename T, size_t N, size_t M> Matrix<T,N,M> Matrix<T,N,M>::operator * 
  (const T mult) const
{
  Matrix<T,N,M> result;

  for (int i = 0; i < N; i++)
  {
    for (int j = 0; j < M; j++)
    {
      result[i][j] = (*this)[i][j] * mult;
    }
  }

  return result;
}

template <typename T, size_t N, size_t M> void Matrix<T,N,M>::operator *= (const T mult)
{
  *this = *this * mult;
}

template <typename T, size_t N, size_t M> T Matrix<T,N,M>::determinant() const
{
  static_assert(N == M,"Determinant is defined only for square matrices");

  T value = static_cast<T>(0.0);

  if (N == 2)
  {
    value = (*this)[0][0] * (*this)[1][1] - (*this)[0][1] * (*this)[1][0];
  } else if (N == 3)
  {
    value = 
      ((*this)[0][0]*((*this)[1][1]*(*this)[2][2]-(*this)[2][1]*(*this)[1][2])) - 
      ((*this)[0][1]*((*this)[1][0]*(*this)[2][2]-(*this)[1][2]*(*this)[2][0])) + 
      ((*this)[0][2]*((*this)[1][0]*(*this)[2][1]-(*this)[1][1]*(*this)[2][0]));
  } else if (N == 4)
  {
    value =
    (*this)[0][3]*(*this)[1][2]*(*this)[2][1]*(*this)[3][0] - (*this)[0][2]*(*this)[1][3]*(*this)[2][1]*(*this)[3][0] - (*this)[0][3]*(*this)[1][1]*(*this)[2][2]*(*this)[3][0] + (*this)[0][1]*(*this)[1][3]*(*this)[2][2]*(*this)[3][0]+
    (*this)[0][2]*(*this)[1][1]*(*this)[2][3]*(*this)[3][0] - (*this)[0][1]*(*this)[1][2]*(*this)[2][3]*(*this)[3][0] - (*this)[0][3]*(*this)[1][2]*(*this)[2][0]*(*this)[3][1] + (*this)[0][2]*(*this)[1][3]*(*this)[2][0]*(*this)[3][1]+
    (*this)[0][3]*(*this)[1][0]*(*this)[2][2]*(*this)[3][1] - (*this)[0][0]*(*this)[1][3]*(*this)[2][2]*(*this)[3][1] - (*this)[0][2]*(*this)[1][0]*(*this)[2][3]*(*this)[3][1] + (*this)[0][0]*(*this)[1][2]*(*this)[2][3]*(*this)[3][1]+
    (*this)[0][3]*(*this)[1][1]*(*this)[2][0]*(*this)[3][2] - (*this)[0][1]*(*this)[1][3]*(*this)[2][0]*(*this)[3][2] - (*this)[0][3]*(*this)[1][0]*(*this)[2][1]*(*this)[3][2] + (*this)[0][0]*(*this)[1][3]*(*this)[2][1]*(*this)[3][2]+
    (*this)[0][1]*(*this)[1][0]*(*this)[2][3]*(*this)[3][2] - (*this)[0][0]*(*this)[1][1]*(*this)[2][3]*(*this)[3][2] - (*this)[0][2]*(*this)[1][1]*(*this)[2][0]*(*this)[3][3] + (*this)[0][1]*(*this)[1][2]*(*this)[2][0]*(*this)[3][3]+
    (*this)[0][2]*(*this)[1][0]*(*this)[2][1]*(*this)[3][3] - (*this)[0][0]*(*this)[1][2]*(*this)[2][1]*(*this)[3][3] - (*this)[0][1]*(*this)[1][0]*(*this)[2][2]*(*this)[3][3] + (*this)[0][0]*(*this)[1][1]*(*this)[2][2]*(*this)[3][3];
  } else
  {
    assert(false && "Determinant is calculated for square matrices 2 x 2, 3 x 3 and 4 x 4 only");
  }

  return value;
}

template <typename T, size_t N, size_t M> Matrix<T,N,M> Matrix<T,N,M>::operator + ()
{
  static_assert(N == M,"Only square matrix can be inverted");

  determinant_ = determinant();
  assert(std::abs(determinant_) >= std::numeric_limits<T>::epsilon() * static_cast<T>(10.0));

  Matrix<T,N,M> res;

  // failure
  if (std::abs(determinant_) < tolerance_)
  {
    errorString_ = "Matrix inversion failed due to close-to-zero determinant";
    res.fill(T(0.0));
    return res;
  }

  T invdet = static_cast<T>(1.0) / determinant_;

  if (N == 2)
  {
    res[0][0] =  +(*this)[1][1];
    res[0][1] =  -(*this)[0][1];
    res[1][0] =  -(*this)[1][0];
    res[1][1] =  +(*this)[0][0];
  } else if (N == 3)
  {
    res[0][0] =  ((*this)[1][1]*(*this)[2][2]-(*this)[1][2]*(*this)[2][1]);
    res[0][1] = -((*this)[0][1]*(*this)[2][2]-(*this)[2][1]*(*this)[0][2]);
    res[0][2] =  ((*this)[0][1]*(*this)[1][2]-(*this)[0][2]*(*this)[1][1]);
    res[1][0] = -((*this)[1][0]*(*this)[2][2]-(*this)[2][0]*(*this)[1][2]);
    res[1][1] =  ((*this)[0][0]*(*this)[2][2]-(*this)[2][0]*(*this)[0][2]);
    res[1][2] = -((*this)[0][0]*(*this)[1][2]-(*this)[0][2]*(*this)[1][0]);
    res[2][0] =  ((*this)[1][0]*(*this)[2][1]-(*this)[2][0]*(*this)[1][1]);
    res[2][1] = -((*this)[0][0]*(*this)[2][1]-(*this)[0][1]*(*this)[2][0]);
    res[2][2] =  ((*this)[0][0]*(*this)[1][1]-(*this)[0][1]*(*this)[1][0]);
  } else if (N == 4)
  {
    res[0][0] = (*this)[1][2]*(*this)[2][3]*(*this)[3][1] - (*this)[1][3]*(*this)[2][2]*(*this)[3][1] + (*this)[1][3]*(*this)[2][1]*(*this)[3][2] - (*this)[1][1]*(*this)[2][3]*(*this)[3][2] - (*this)[1][2]*(*this)[2][1]*(*this)[3][3] + (*this)[1][1]*(*this)[2][2]*(*this)[3][3];
    res[0][1] = (*this)[0][3]*(*this)[2][2]*(*this)[3][1] - (*this)[0][2]*(*this)[2][3]*(*this)[3][1] - (*this)[0][3]*(*this)[2][1]*(*this)[3][2] + (*this)[0][1]*(*this)[2][3]*(*this)[3][2] + (*this)[0][2]*(*this)[2][1]*(*this)[3][3] - (*this)[0][1]*(*this)[2][2]*(*this)[3][3];
    res[0][2] = (*this)[0][2]*(*this)[1][3]*(*this)[3][1] - (*this)[0][3]*(*this)[1][2]*(*this)[3][1] + (*this)[0][3]*(*this)[1][1]*(*this)[3][2] - (*this)[0][1]*(*this)[1][3]*(*this)[3][2] - (*this)[0][2]*(*this)[1][1]*(*this)[3][3] + (*this)[0][1]*(*this)[1][2]*(*this)[3][3];
    res[0][3] = (*this)[0][3]*(*this)[1][2]*(*this)[2][1] - (*this)[0][2]*(*this)[1][3]*(*this)[2][1] - (*this)[0][3]*(*this)[1][1]*(*this)[2][2] + (*this)[0][1]*(*this)[1][3]*(*this)[2][2] + (*this)[0][2]*(*this)[1][1]*(*this)[2][3] - (*this)[0][1]*(*this)[1][2]*(*this)[2][3];
    res[1][0] = (*this)[1][3]*(*this)[2][2]*(*this)[3][0] - (*this)[1][2]*(*this)[2][3]*(*this)[3][0] - (*this)[1][3]*(*this)[2][0]*(*this)[3][2] + (*this)[1][0]*(*this)[2][3]*(*this)[3][2] + (*this)[1][2]*(*this)[2][0]*(*this)[3][3] - (*this)[1][0]*(*this)[2][2]*(*this)[3][3];
    res[1][1] = (*this)[0][2]*(*this)[2][3]*(*this)[3][0] - (*this)[0][3]*(*this)[2][2]*(*this)[3][0] + (*this)[0][3]*(*this)[2][0]*(*this)[3][2] - (*this)[0][0]*(*this)[2][3]*(*this)[3][2] - (*this)[0][2]*(*this)[2][0]*(*this)[3][3] + (*this)[0][0]*(*this)[2][2]*(*this)[3][3];
    res[1][2] = (*this)[0][3]*(*this)[1][2]*(*this)[3][0] - (*this)[0][2]*(*this)[1][3]*(*this)[3][0] - (*this)[0][3]*(*this)[1][0]*(*this)[3][2] + (*this)[0][0]*(*this)[1][3]*(*this)[3][2] + (*this)[0][2]*(*this)[1][0]*(*this)[3][3] - (*this)[0][0]*(*this)[1][2]*(*this)[3][3];
    res[1][3] = (*this)[0][2]*(*this)[1][3]*(*this)[2][0] - (*this)[0][3]*(*this)[1][2]*(*this)[2][0] + (*this)[0][3]*(*this)[1][0]*(*this)[2][2] - (*this)[0][0]*(*this)[1][3]*(*this)[2][2] - (*this)[0][2]*(*this)[1][0]*(*this)[2][3] + (*this)[0][0]*(*this)[1][2]*(*this)[2][3];
    res[2][0] = (*this)[1][1]*(*this)[2][3]*(*this)[3][0] - (*this)[1][3]*(*this)[2][1]*(*this)[3][0] + (*this)[1][3]*(*this)[2][0]*(*this)[3][1] - (*this)[1][0]*(*this)[2][3]*(*this)[3][1] - (*this)[1][1]*(*this)[2][0]*(*this)[3][3] + (*this)[1][0]*(*this)[2][1]*(*this)[3][3];
    res[2][1] = (*this)[0][3]*(*this)[2][1]*(*this)[3][0] - (*this)[0][1]*(*this)[2][3]*(*this)[3][0] - (*this)[0][3]*(*this)[2][0]*(*this)[3][1] + (*this)[0][0]*(*this)[2][3]*(*this)[3][1] + (*this)[0][1]*(*this)[2][0]*(*this)[3][3] - (*this)[0][0]*(*this)[2][1]*(*this)[3][3];
    res[2][2] = (*this)[0][1]*(*this)[1][3]*(*this)[3][0] - (*this)[0][3]*(*this)[1][1]*(*this)[3][0] + (*this)[0][3]*(*this)[1][0]*(*this)[3][1] - (*this)[0][0]*(*this)[1][3]*(*this)[3][1] - (*this)[0][1]*(*this)[1][0]*(*this)[3][3] + (*this)[0][0]*(*this)[1][1]*(*this)[3][3];
    res[2][3] = (*this)[0][3]*(*this)[1][1]*(*this)[2][0] - (*this)[0][1]*(*this)[1][3]*(*this)[2][0] - (*this)[0][3]*(*this)[1][0]*(*this)[2][1] + (*this)[0][0]*(*this)[1][3]*(*this)[2][1] + (*this)[0][1]*(*this)[1][0]*(*this)[2][3] - (*this)[0][0]*(*this)[1][1]*(*this)[2][3];
    res[3][0] = (*this)[1][2]*(*this)[2][1]*(*this)[3][0] - (*this)[1][1]*(*this)[2][2]*(*this)[3][0] - (*this)[1][2]*(*this)[2][0]*(*this)[3][1] + (*this)[1][0]*(*this)[2][2]*(*this)[3][1] + (*this)[1][1]*(*this)[2][0]*(*this)[3][2] - (*this)[1][0]*(*this)[2][1]*(*this)[3][2];
    res[3][1] = (*this)[0][1]*(*this)[2][2]*(*this)[3][0] - (*this)[0][2]*(*this)[2][1]*(*this)[3][0] + (*this)[0][2]*(*this)[2][0]*(*this)[3][1] - (*this)[0][0]*(*this)[2][2]*(*this)[3][1] - (*this)[0][1]*(*this)[2][0]*(*this)[3][2] + (*this)[0][0]*(*this)[2][1]*(*this)[3][2];
    res[3][2] = (*this)[0][2]*(*this)[1][1]*(*this)[3][0] - (*this)[0][1]*(*this)[1][2]*(*this)[3][0] - (*this)[0][2]*(*this)[1][0]*(*this)[3][1] + (*this)[0][0]*(*this)[1][2]*(*this)[3][1] + (*this)[0][1]*(*this)[1][0]*(*this)[3][2] - (*this)[0][0]*(*this)[1][1]*(*this)[3][2];
    res[3][3] = (*this)[0][1]*(*this)[1][2]*(*this)[2][0] - (*this)[0][2]*(*this)[1][1]*(*this)[2][0] + (*this)[0][2]*(*this)[1][0]*(*this)[2][1] - (*this)[0][0]*(*this)[1][2]*(*this)[2][1] - (*this)[0][1]*(*this)[1][0]*(*this)[2][2] + (*this)[0][0]*(*this)[1][1]*(*this)[2][2];
  } else
  {
    assert(false && "Only square matrices 2 x 2, 3 x 3 and 4 x 4 can be inverted");
    res.fill(static_cast<T>(0.0));
  }

  res *= invdet;

  return res;
}

template <typename T, size_t N, size_t M> bool Matrix<T,N,M>::isIdentity() const
{ 
  if (N != M)
    return false;

  for (int i = 0; i < N; i++)
  {
    for (int j = 0; j < M; j++)
    {
      if (i == j)
      {
        if (std::abs((*this)[i][j] - static_cast<T>(1.0)) > tolerance_)
          return false;
      } else
      {
        if (std::abs((*this)[i][j]) > tolerance_)
          return false;
      }
    }
  }

  return true;
}

template <typename T, size_t N, size_t M> void Matrix<T,N,M>::setDiag(const T value)
{
  static_assert(N == M,"Setting diagonal requires square matrix");

  for (int i = 0; i < N; i++)
  {
    for (int j = 0; j < M; j++)
    {
      if (i == j)
        (*this)[i][j] = value;
    }
  }
}

template <typename T, size_t N, size_t M> void Matrix<T,N,M>::setIdentity()
{
  static_assert(N == M,"Setting identity requires square matrix");

  fill(static_cast<T>(0.0));

  setDiag(static_cast<T>(1.0));
}

template <typename T, size_t N, size_t M> T Matrix<T,N,M>::lastDeterminant() const
{
  return determinant_;
}
