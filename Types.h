#pragma once

// LINT is always 64-bit
#define LINT long long int

// undefined
#define NOT_DEFINED -1
#define NOT_FOUND NOT_DEFINED

// useful macros
#define LIMIT(x,xmin,xmax) if (x < xmin) x = xmin; if (x > xmax) x = xmax
#define LIMIT_MIN(x,xmin) if (x < xmin) x = xmin
#define LIMIT_MAX(x,xmax) if (x > xmax) x = xmax

// tolerance
#define TOLERANCE(T) std::numeric_limits<T>::epsilon() * static_cast<T>(10.0)
#define TOLERANCE4 TOLERANCE(float)
#define TOLERANCE8 TOLERANCE(double)

// round
#define ROUND(x)  (int) (floor(x + .5))

//===== PI - associated ========================================================
#ifndef M_PI
  #define M_PI    3.14159265358979323846
#endif

#ifndef PI05
  #define PI05 (M_PI * 0.5)
#endif

#ifndef PI20
  #define PI20 (M_PI * 2.0)
#endif

#ifndef PI10
  #define PI10 M_PI
#endif

#ifndef PI40
  #define PI40 (M_PI * 4.0)
#endif


#ifndef PCI
  #define PCI (180.0 / M_PI)
#endif

#ifndef CPI
  #define CPI (M_PI / 180.0)
#endif

#ifndef MI_RAD
  #define MI_RAD  (M_PI / 10800.)
#endif

#ifndef RAD_MI
  #define RAD_MI  (10800.0 / M_PI)
#endif

#ifndef MI_45
  #define	MI_45	  2700.0
#endif

// swap
#define SWAP(T,x1,x2) { T temp = x1; x1 = x2; x2 = temp; }

//===== constants ==============================================================

#define T0 static_cast<T>(0.0)
#define T1 static_cast<T>(1.0)
#define T2 static_cast<T>(2.0)
#define T3 static_cast<T>(3.0)
#define T6 static_cast<T>(6.0)
#define T9 static_cast<T>(9.0)
#define T10 static_cast<T>(10.0)
#define T12 static_cast<T>(12.0)

#define T05 static_cast<T>(0.5)
#define T025 static_cast<T>(0.25)
#define T0125 static_cast<T>(0.125)
#define T13 static_cast<T>(0.33333333333333)
#define T23 static_cast<T>(0.66666666666667)


