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

#include <assert.h>


//===== Gauss numerical integration ============================================
// exists in doubles only

                              // for rectangles
#define MAX_GAUSSPOINTS 20

template <class T> struct TGaussInt {
  int numpoints;
  T knots[MAX_GAUSSPOINTS];
  T weights[MAX_GAUSSPOINTS];
};

typedef TGaussInt<double> SGaussInt;

enum {
  GAUSSINT_1,
  GAUSSINT_2,
  GAUSSINT_4,
  GAUSSINT_8,
  GAUSSINT_20,
  GAUSSINT_4_CONFORMING,
  GAUSSINT_TOTAL,
  OTHER_INTEGRATION = GAUSSINT_TOTAL
};

static SGaussInt GaussInt[GAUSSINT_TOTAL] =
{
  {1,
    {0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0},  
    {2.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0}},
  {2,
    {-0.333333333333333,+0.333333333333333,0.0,0.0,0.0,0.0,0.0,0.0},
    {1.0,1.0,0.0,0.0,0.0,0.0,0.0,0.0}},
  {4,
    {-0.861136311594053,-0.339981043584856,+0.339981043584856,+0.861136311594053,0.0,0.0,0.0,0.0},
    {+0.347854845137454,+0.652145154862546,+0.652145154862546,+0.347854845137454,0.0,0.0,0.0,0.0}},
  {8,
    {-0.9602898564975363,-0.7966664774136267,-0.5255324099163290,-0.1834346424956498,
     +0.1834346424956498,+0.5255324099163290,+0.7966664774136267,+0.9602898564975363},
    {+0.1012285362903763,+0.2223810344533745,+0.3137066458778873,+0.3626837833783620,
     +0.3626837833783620,+0.3137066458778873,+0.2223810344533745,+0.1012285362903763}},
  {20,
    {
    -0.0765265211334973,
    0.0765265211334973,
    -0.2277858511416451,
    0.2277858511416451,
    -0.3737060887154195,
    0.3737060887154195,
    -0.5108670019508271,
    0.5108670019508271,
    -0.6360536807265150,
    0.6360536807265150,
    -0.7463319064601508,
    0.7463319064601508,
    -0.8391169718222188,
    0.8391169718222188,
    -0.9122344282513259,
    0.9122344282513259,
    -0.9639719272779138,
    0.9639719272779138,
    -0.9931285991850949,
    0.9931285991850949
    },

    {
    0.1527533871307258,
    0.1527533871307258,
    0.1491729864726037,
    0.1491729864726037,
    0.1420961093183820,
    0.1420961093183820,
    0.1316886384491766,
    0.1316886384491766,
    0.1181945319615184,
    0.1181945319615184,
    0.1019301198172404,
    0.1019301198172404,
    0.0832767415767048,
    0.0832767415767048,
    0.0626720483341091,
    0.0626720483341091,
    0.0406014298003869,
    0.0406014298003869,
    0.0176140071391521,
    0.0176140071391521
    }},

  {4,
    {-0.416666666665,-0.083333333335,+0.08333333335,0.4166666666665,0.0,0.0,0.0,0.0},
    {0.5,0.5,0.5,0.5,0.0,0.0,0.0,0.0}}
};

                              // for triangles
#define MAX_HAMMERPOINTS 16

template <class T> struct THammerInt {
  int numpoints;
                              // e3 = 1 - e1 - e2
  T e1e2w[MAX_HAMMERPOINTS][3];
};

typedef THammerInt<double> SHammerInt;

enum {
  HAMMERINT_1,
  HAMMERINT_3,
  HAMMERINT_4,
  HAMMERINT_7,
  HAMMERINT_16,
  HAMMERINT_TOTAL
};

static SHammerInt HammerInt[HAMMERINT_TOTAL] =
{
  {1,{
    0.33333333333333, 0.33333333333333, 1.00000000000000
  }},

  {3,{
    0.16666666666667, 0.16666666666667, 0.33333333333333,
    0.16666666666667, 0.66666666666667, 0.33333333333333,
    0.66666666666667, 0.16666666666667, 0.33333333333333,
  }},

  {4,{
    0.33333333333333, 0.33333333333333, -0.5625000000000,
    0.20000000000000, 0.20000000000000, 0.52083333333333,
    0.20000000000000, 0.60000000000000, 0.52083333333333,
    0.60000000000000, 0.20000000000000, 0.52083333333333
  }},

  {7,{
    0.33333333333333, 0.33333333333333, 0.22500000000000,
    0.47014206410511, 0.47014206410511, 0.13239415278851,
    0.47014206410511, 0.05971587178977, 0.13239415278851,
    0.05971587178977, 0.47014206410511, 0.13239415278851,
    0.10128650732346, 0.10128650732346, 0.12593918054483,
    0.10128650732346, 0.79742698535309, 0.12593918054483,
    0.79742698535309, 0.10128650732346, 0.12593918054483
  }},

  {16,{
    0.33333333333333, 0.33333333333333, 0.14431560767779,
    0.45929258829272, 0.45929258829272, 0.09509163426728,
    0.45929258829272, 0.08141482341455, 0.09509163426728,
    0.08141482341455, 0.45929258829272, 0.09509163426728,
    0.17056930775176, 0.17056930775176, 0.10321737053472,
    0.17056930775176, 0.65886138449648, 0.10321737053472,
    0.65886138449648, 0.17056930775176, 0.10321737053472,
    0.05054722831703, 0.05054722831703, 0.03245849762320,
    0.05054722831703, 0.89890554336594, 0.03245849762320,
    0.89890554336594, 0.05054722831703, 0.03245849762320,
    0.26311282963464, 0.72849239295540, 0.02723031417443,
    0.72849239295540, 0.00839477740996, 0.02723031417443,
    0.00839477740996, 0.26311282963464, 0.02723031417443,
    0.72849239295540, 0.26311282963464, 0.02723031417443,
    0.26311282963464, 0.00839477740996, 0.02723031417443,
    0.00839477740996, 0.72849239295540, 0.02723031417443
  }}
};

                              // for triangles with singularity at U = 1
enum {
  SINGULAR_1,
  SINGULAR_3,
  SINGULAR_4,
  SINGULAR_TOTAL
};

static SHammerInt SingularInt[SINGULAR_TOTAL] =
{
  {1,{
    0.5, 0.25, 1.24645048
  }},

  {3,{
    0.66666667, 0.16666667, 0.93483790,
    0.00000000, 0.81742619, 0.15580629,
    0.00000000, 0.18257381, 0.15580629
  }},

  {4,{
    0.78857548, 0.16385495, 0.31161231,
    0.21132509, 0.61114353, 0.31161293,
    0.78857548, 0.04756957, 0.31161231,
    0.21132509, 0.17753138, 0.31161293
  }}

};

