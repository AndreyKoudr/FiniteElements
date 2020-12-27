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
                              // intrinsics
#include <xmmintrin.h>
#include <emmintrin.h>
#include <limits>
#include <assert.h>
#include "Types.h"

/**
  These are two basic classes for 3D vectors with operators which make
possible to use simple operations on vectors like v3 = v1 + v2 etc.

  The first class is a traditional C without SIMD acceleration; the second is
based on basic SSE2 Intel intrinsics 
(see https://software.intel.com/sites/landingpage/IntrinsicsGuide/). The intrinsics
together with the optimiser may produce very good acceleration of the code
(e.g. https://github.com/AndreyKoudr/ThreadedGaussElimination),
but it is a little bit tricky.

  To make a fast code, you must follow the two principles
  (1) more calculations (in registers)
  (2) less memory access because it is much slower

It means that such constructions without temporary variables like these

    v3.data = _mm_sub_ps(
        _mm_mul_ps(_mm_shuffle_ps(data, data, _MM_SHUFFLE(3, 0, 2, 1)),
          _mm_shuffle_ps(v2.data, v2.data, _MM_SHUFFLE(3, 1, 0, 2))),
        _mm_mul_ps(_mm_shuffle_ps(data, data, _MM_SHUFFLE(3, 1, 0, 2)),
          _mm_shuffle_ps(v2.data, v2.data, _MM_SHUFFLE(3, 0, 2, 1))));
or 

v3[i] = 6 * (v3[i] + v4[i]) / (!(v3[i] + v4[i])) + (+(v4[i] - v3[i])) * static_cast<float>(0.333); 

are preferable.
  
  The first template class TVector<T> works with both float-s and double-s; the second
SIMD-base class Vector4 is written for 4-byte floats only.

  The sets of operations are identical for both classes.

  The fourth homogeneous coordinate is set to zero everywhere but normally it should be 
initialised to 1.0.
*/

/**
  4-component vector template, no SIMD. The fourth component is a homogeneous 
coordinate.
  It does not make much sense to make a template for a N-component (e.g. 3) 
vector as currently there are such things as SIMD enabling to make operations
on all 4 components at once.
 
  Operations
  ----------
  "float" below is 4 or 8 byte float

  float = V[AxisX]  - get vector coordinate
  float = !V        - get vector length                         
  v3 = v1 ^ v2      - cross product                     
  v2 = +v1          - normalisation
  float = v1 * v2   - dot product, W is ignored
  (v1 > v2)         - vectors co-directed?    
  bool(v1==v2)      - vectors equal?                            
  v2 = -v1          - change sign of components                 
  v3 = v1 + v2      - addition
  v1 += v2;         - addition
  v3 = v1 - v2      - subtraction
  v1 -= v2;         - subtraction
  v3 = v1 * float   - multiply by scalar
  v3 = float * v    - multiply by scalar
  v1 *= float       - multiply by scalar
  v3 = v1 / float   - divide by scalar
  v1 /= float       - divide by scalar
*/

                              // axes
enum Axes {AxisX,AxisY,AxisZ,AxisW};                  

template <class T> class TVector {
public:
                              // data
  union {
    struct { T X,Y,Z,W; };
    T XYZW[4];
    T XYZ[3];
  };
                              // constructors
  TVector() 
  { 
                              // better to avoid any initialisation here as it takes time
    init();
  }
                              // initialise to zero vector
  void init()
  {
    X = Y = Z = 0.0; 
    W = 0.0;
  }

  TVector(T x, T y = 0.0, T z = 0.0, T w = 0.0) 
  { 
    X = x; Y = y; Z = z; W = w; 
  }

  TVector(const TVector<T>& other) 
  { 
    X = other.X; 
    Y = other.Y; 
    Z = other.Z; 
    W = other.W; 
  }

  TVector operator=(const TVector<T>& other) 
  { 
    X = other.X; 
    Y = other.Y; 
    Z = other.Z; 
    W = other.W; 

    return *this;
  }
                              // get coordinate
  T& operator[](Axes A) {
    switch (A) {
    case AxisX : return X ;
    case AxisY : return Y ;
    case AxisZ : return Z ;
    case AxisW : return W ;
    default : assert(false && "Wrong axis value");
    };
  }
                              // get coordinate
  T& operator[](size_t index) {
    assert(index >= 0 && index <= 3);
    switch (index) {
    case 0 : return X ;
    case 1 : return Y ;
    case 2 : return Z ;
    case 3 : return W ;
    default : return X;
    };
  }
                              // operator ! - get vector length
  T operator!() const
  {
    return static_cast<T>(sqrt(X * X + Y * Y + Z * Z));
  }
                              // length
  T length() const
  {
    return static_cast<T>(sqrt(X * X + Y * Y + Z * Z));
  }
                              // length squared
  T length2() const
  {
    return X * X + Y * Y + Z * Z;
  }
                              // operator ^ for two vectors (cross product)
  TVector operator^(TVector v2) const
  {
    TVector v3;
		v3.X = Y*v2.Z-Z*v2.Y;
		v3.Y = Z*v2.X-X*v2.Z;
		v3.Z = X*v2.Y-Y*v2.X;
    return v3;
  }
                              // operator + - normalization
  TVector operator+() {
    TVector v3;
    T Len = static_cast<T>(sqrt(X*X+Y*Y+Z*Z));

		if (Len > tolerance_)
		{
			v3.X = X/Len;
			v3.Y = Y/Len;
			v3.Z = Z/Len;
		} else
		{
      v3.X = 0;
      v3.Y = 0;
      v3.Z = 0;
		}
    return v3;
  }
                              // dot product, W is ignored
  T operator*(TVector v2) const {
    return X * v2.X + Y * v2.Y + Z * v2.Z;
  }
                              // operator > - two vectors co-directed?
  bool operator>(TVector v2) {
    return ((*this * v2) >= 0);
  }
                              // operator == - two vectors equal?
  bool operator==(TVector v2) const {
    T dlen = !(*this - v2);
    return (dlen < tolerance_);
  }
                              // operator - (change direction)
  TVector operator-() {
    return (*this) * (-1.0);
  }
                              // addition
  TVector operator+=(TVector v2) {
    for (int i = 0; i < 4; i++) XYZW[i] += v2.XYZW[i];
    return *this;
  }

  TVector inline operator+(TVector v2) const {
    TVector v3;
    for (int i = 0; i < 4; i++) v3.XYZW[i] = XYZW[i] + v2.XYZW[i];
    return v3;
  }
			                        // subtraction
  TVector operator-(TVector v2) const {
    TVector v3;
    for (int i = 0; i < 4; i++) v3.XYZW[i] = XYZW[i] - v2.XYZW[i];
    return v3;
  }
                              // subtract
  TVector operator-=(TVector v2) {
    for (int i = 0; i < 4; i++) XYZW[i] -= v2.XYZW[i];
    return *this;
  }
                              // operator *scalar
  TVector operator*(T multiplier) const {
    TVector v3;
    v3.X = X*multiplier;
    v3.Y = Y*multiplier;
    v3.Z = Z*multiplier;
    v3.W = W*multiplier;
    return v3;
  }
                              // operator *= scalar
  TVector operator*=(T multiplier) {
    X *= multiplier;
    Y *= multiplier;
    Z *= multiplier;
    W *= multiplier;
    return *this;
  }
                              // operator / scalar
  TVector operator/(T divider) const {
    TVector v3;
    assert(std::abs(divider) > tolerance_ && "Division by zero in TVector");
    v3.X = X/divider;
    v3.Y = Y/divider;
    v3.Z = Z/divider;
    v3.W = W/divider;
    return v3;
  }
                              // operator /= scalar
  TVector operator/=(T divider) {
    assert(std::abs(divider) > tolerance_ && "Division by zero in TVector");
    X /= divider;
    Y /= divider;
    Z /= divider;
    W /= divider;
    return *this;
  }
                              // operator /= vector
  TVector operator/=(TVector divider) {
    assert((std::abs(divider.X) > tolerance_ && std::abs(divider.Y) > tolerance_ &&
      std::abs(divider.Z) > tolerance_ && std::abs(divider.W) > tolerance_) && "Division by zero in TVector");
    X /= divider.X;
    Y /= divider.Y;
    Z /= divider.Z;
    W /= divider.W;
    return *this;
  }

  static TVector zero()
  {
    return TVector(0,0,0,0);
  }

//!!!
															// get angle in radians between this vector and v, 
                              // measured counter-clockwise; plane normal is
                              // plane with normal to make rotation comparisons
	T GetAngleSigned(TVector planennormal, TVector v)
	{
		TVector p = (+(*this)) ^ (+v);
    T plen = !p;
    //assert(plen <= T1);
    LIMIT_MAX(plen,T1);
    T a = fabs(asin(plen));

    if (!(p > planennormal))
      a = -a;

		return a;
	}

  // tolerance
  static T tolerance()
  {
    return tolerance_;
  }

  static void setTolerance(const T tolerance)
  {
    tolerance_ = tolerance;
  }

  static void setDefaultTolerance(const T tolerance)
  {
    tolerance_ = std::numeric_limits<T>::epsilon() * static_cast<T>(100.0);
  }

private:
  // tolerance
  static T tolerance_;
};

template <typename T>
T TVector<T>::tolerance_ = std::numeric_limits<T>::epsilon() * static_cast<T>(100.0);

                              // scalar * vector
template <class T> TVector<T> operator*(const double scalar, const TVector<T> v)
{
  return v * static_cast<T>(scalar);
}



/**
  4-component vector on 4-byte floats driven by SIMD. The fourth component is a 
homogeneous coordinate.
 
  Operations
  ----------
  "float" below is a 4-byte float

  float = V[AxisX]  - get vector coordinate
  float = !V        - get vector length                         
  v3 = v1 ^ v2      - cross product                     
  v2 = +v1          - normalisation
  float = v1 * v2   - dot product, W is ignored
  (v1 > v2)         - vectors co-directed?    
  bool(v1==v2)      - vectors equal?                            
  v2 = -v1          - change sign of components                 
  v3 = v1 + v2      - addition
  v1 += v2;         - addition
  v3 = v1 - v2      - subtraction
  v1 -= v2;         - subtraction
  v3 = v1 * float   - multiply by scalar
  v3 = float * v    - multiply by scalar
  v1 *= float       - multiply by scalar
  v3 = v1 / float   - divide by scalar
  v1 /= float       - divide by scalar
*/

#pragma pack(push,4)

class Vector4 {
public:
                              // data
  union {
    __m128 data;
    struct { float X,Y,Z,W; };
    float XYZW[4];
    float XYZ[3];
  };
                              // default constructor
  Vector4()
  { 
                              // better to avoid any initialisation here as it takes time
    init();
  }                             
                              // initialise to zero vector
  void init()
  {
    data = _mm_setzero_ps();
  }
                              // constructor
  Vector4(float x, float y = 0, float z = 0, float w = 0) 
  { 
    data = _mm_setr_ps(x,y,z,w);
  }
                              // destructor, does nothing
  ~Vector4() = default;
                              // copy constructor
  Vector4(const Vector4& other) 
  { 
    data = other.data;
  }
                              // assignment
  Vector4 operator=(const Vector4& other) 
  { 
    data = other.data;
    return *this;
  }
                              // get coordinate
  float& operator[](Axes A) {
    switch (A) {
    case AxisX : return X;
    case AxisY : return Y ;
    case AxisZ : return Z ;
    case AxisW : return W ;
    default : assert(false && "Wrong axis value");
    };
  }
                              // get coordinate
  float& operator[](size_t index) {
    switch (index) {
    case 0 : return X ;
    case 1 : return Y ;
    case 2 : return Z ;
    case 3 : return W ;
    default : { assert(false && "Wrong axis value"); return X; };
    };
  }
                              // operator ! - get vector length
  float operator!() const {
    __m128 c1,temp;

    // x * x, y * y, z * z, 0 * 0 (X,Y,Z,W)
    c1 = _mm_mul_ps(data,data);

    // Shift two times and add 0-value to c1
    temp = _mm_shuffle_ps(c1,c1,0x21);
    c1 = _mm_add_ss(c1,temp);
    temp = _mm_shuffle_ps(temp,temp,0x02);
    c1 = _mm_add_ss(c1,temp);

    // Spread 0-value over all register
    c1 = _mm_sqrt_ps(_mm_shuffle_ps(c1,c1,0x00));

    return c1.m128_f32[0];
  }
                              // length
  float length() const
  {
    return !(*this);
  }
                              // length squared
  float length2() const
  {
    return (*this) * (*this);
  }
                              // operator ^ for two vectors (cross product)
  Vector4 operator^(Vector4 v2) const
  {
    Vector4 v3;

    v3.data = _mm_sub_ps(
        _mm_mul_ps(_mm_shuffle_ps(data, data, _MM_SHUFFLE(3, 0, 2, 1)),
          _mm_shuffle_ps(v2.data, v2.data, _MM_SHUFFLE(3, 1, 0, 2))),
        _mm_mul_ps(_mm_shuffle_ps(data, data, _MM_SHUFFLE(3, 1, 0, 2)),
          _mm_shuffle_ps(v2.data, v2.data, _MM_SHUFFLE(3, 0, 2, 1))));

    return v3;
  }
                              // operator + - normalization
  Vector4 operator+() const
  {
    Vector4 v3;

    __m128 c1,temp;
                              // x * x, y * y, z * z, 0 * 0 (X,Y,Z,W)
    c1 = _mm_mul_ps(data,data);
                              // shift two times and add 0-value to c1
    temp = _mm_shuffle_ps(c1,c1,0x21);
    c1 = _mm_add_ss(c1,temp);
    temp = _mm_shuffle_ps(temp,temp,0x02);
    c1 = _mm_add_ss(c1,temp);
                              // spread 0-value over all register and get sqrt
                              // divide by length
    c1 = _mm_div_ps(data,_mm_sqrt_ps(_mm_shuffle_ps(c1,c1,0x00)));
                              // set result to zero
    v3.data = _mm_setzero_ps();
                              // if division is not NaN or INF, add to result
    v3.data = _mm_and_ps(_mm_cmpeq_ps(c1,c1),c1);

    return v3;
  }
                              // dot product, W is ignored
  float operator*(Vector4 v2) const 
  {
     __m128 c1,temp;
                              // x * x, y * y, z * z, 0 * 0 (X,Y,Z,W)
    c1 = _mm_mul_ps(data,v2.data);
                              // shift two times and add 0-value to c1
    temp = _mm_shuffle_ps(c1,c1,0x21);
    c1 = _mm_add_ss(c1,temp);
    temp = _mm_shuffle_ps(temp,temp,0x02);
    c1 = _mm_add_ss(c1,temp);
                              // spread 0-value over all register
    c1 = _mm_shuffle_ps(c1,c1,0x00);

    return c1.m128_f32[0];
  }
                              // operator > - two vectors co-directed?
  bool operator>(Vector4 v2) 
  {
    return ((*this * v2) >= 0);
  }
                              // operator == - two vectors equal?
  bool operator==(Vector4 v2) const 
  {
    float dlen = !(*this - v2);
    return (dlen < TOLERANCE4);
  }
                              // operator - (change direction)
  Vector4 operator-() {
    return (*this) * (-1);
  }
                              // addition
  Vector4 inline operator+(Vector4 v2) const 
  {
    Vector4 v3;
    v3.data = _mm_add_ps(data,v2.data);
    return v3;
  }
                              // addition
  Vector4 operator+=(Vector4 v2) 
  {
    data = _mm_add_ps(data,v2.data);
    return *this;
  }
                              // subtraction
  Vector4 operator-(Vector4 v2) const 
  {
    Vector4 v3;
    v3.data = _mm_sub_ps(data,v2.data);
    return v3;
  }
                              // subtraction
  Vector4 operator-=(Vector4 v2) 
  {
    data = _mm_sub_ps(data,v2.data);
    return *this;
  }
                              // operator *scalar
  Vector4 operator*(float scalar) const 
  {
    Vector4 v3;
    v3.data = _mm_mul_ps(data,_mm_set1_ps(scalar));
    return v3;
  }
                              // operator *= scalar
  Vector4 operator*=(float scalar) 
  {
    data = _mm_mul_ps(data,_mm_set1_ps(scalar));
    return *this;
  }
                              // operator /scalar
  Vector4 operator/(float scalar) const 
  {
    assert(std::abs(scalar) > TOLERANCE4 && "Division by zero in Vector4");

    Vector4 v3;
    v3.data = _mm_div_ps(data,_mm_set1_ps(scalar));
    return v3;
  }
                              // operator /= scalar
  Vector4 operator/=(float scalar) 
  {
    assert(std::abs(scalar) > TOLERANCE4 && "Division by zero in Vector4");

    data = _mm_div_ps(data,_mm_set1_ps(scalar));
    return *this;
  }
                              // operator /= vector
  Vector4 operator/=(Vector4 divider) {
    assert((std::abs(divider.X) > TOLERANCE4 && std::abs(divider.Y) > TOLERANCE4 &&
      std::abs(divider.Z) > TOLERANCE4 && std::abs(divider.W) > TOLERANCE4) && "Division by zero in Vector4");

    data = _mm_div_ps(data,divider.data);

    return *this;
  }

  static Vector4 zero()
  {
    return Vector4(float(0.0),float(0.0),float(0.0),float(0.0));
  }

  static float tolerance()
  {
    return tolerance_;
  }

  static void setTolerance(const float tolerance)
  {
    tolerance_ = tolerance;
  }

  static void setDefaultTolerance(const float tolerance)
  {
    tolerance_ = std::numeric_limits<float>::epsilon() * static_cast<float>(100.0);
  }

private:
  // tolerance
  static float tolerance_;
};

float Vector4::tolerance_ = std::numeric_limits<float>::epsilon() * static_cast<float>(100.0);

                              // scalar * vector
Vector4 operator*(const double scalar, const Vector4 v);


#pragma pack(pop)