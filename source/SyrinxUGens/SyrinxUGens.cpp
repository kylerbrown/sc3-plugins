/*
  SuperCollider real time audio synthesis system
  Copyright (c) 2002 James McCartney. All rights reserved.
  http://www.audiosynth.com
  This program is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; either version 2 of the License, or
  (at your option) any later version.
  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
  GNU General Public License for more details.
  You should have received a copy of the GNU General Public License
  along with this program; if not, write to the Free Software
  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA

  Syrinx ugen created by Kyler Brown, 2014, based off LorenzL from
  the Chaos Ugens.

  Chaos Ugens created by Lance Putnam on Mon Jul 19 2004.
*/
#include "SC_PlugIn.h"
#define ONESIXTH 0.1666666666666667
static InterfaceTable *ft;
struct NonLinear : public Unit {
  double x0, y0, xn, yn, xnm1, ynm1;
  float counter;
  //bool stable;
};

struct SyrinxN : public NonLinear {
};

struct SyrinxL : public SyrinxN {
  double frac;
};

extern "C" {
  void SyrinxL_next(SyrinxL *unit, int inNumSamples);
  void SyrinxL_Ctor(SyrinxL *unit);
}

////////////////////////////////////////////////////////////////////////////////
void SyrinxL_next(SyrinxL *unit, int inNumSamples)
{
  float *out = ZOUT(0);
  float freq = ZIN0(0);
  double a = ZIN0(1);
  double b = ZIN0(2);
  double g = ZIN0(3);
  double h = ZIN0(4);
  double x0 = ZIN0(5);
  double y0 = ZIN0(6);
  double xn = unit->xn;
  double yn = unit->yn;
  float counter = unit->counter;
  double xnm1 = unit->xnm1;
  double ynm1 = unit->ynm1;
  double frac = unit->frac;
  float samplesPerCycle;
  double slope;
  if(freq < unit->mRate->mSampleRate){
    samplesPerCycle = unit->mRate->mSampleRate / sc_max(freq, 0.001f);
    slope = 1.f / samplesPerCycle;
  }
  else {
    samplesPerCycle = 1.f;
    slope = 1.f;
  }
  if((unit->x0 != x0) || (unit->y0 != y0)){
    xnm1 = xn;
    ynm1 = yn;
    unit->x0 = xn = x0;
    unit->y0 = yn = y0;
  }
  double dx = xn - xnm1;
  for (int i=0; i<inNumSamples; ++i) {
    if(counter >= samplesPerCycle){
      counter -= samplesPerCycle;
      frac = 0.f;
      xnm1 = xn;
      ynm1 = yn;
      double k1x, k2x, k3x, k4x,
	k1y, k2y, k3y, k4y,
	kxHalf, kyHalf;
      double gg = g * g;
      double agg = a * g * g;
      double bgg = b * g * g;
      // 4th order Runge-Kutta
      k1x = h * ynm1;
      k1y = h * (-agg - bgg*xnm1 - gg*xnm1*xnm1*xnm1 - g*xnm1*xnm1*ynm1
		 + gg*xnm1*xnm1 - g*xnm1*ynm1);
      kxHalf = k1x * 0.5;
      kyHalf = k1y * 0.5;
      //k2x = hTimesS * (ynm1 + kyHalf - xnm1 - kxHalf);
      k2x = h * (ynm1 + kyHalf);
      //k2y = h * ((xnm1 + kxHalf) * (r - znm1 - kzHalf) - (ynm1 + kyHalf));
      k2y = h * (-agg - bgg*(xnm1+kxHalf) - gg*(xnm1+kxHalf)*(xnm1+kxHalf)*(xnm1+kxHalf)
		 - g*(xnm1+kxHalf)*(xnm1+kxHalf)*(ynm1+kyHalf)
		 + gg*(xnm1+kxHalf)*(xnm1+kxHalf) - g*(xnm1+kxHalf)*(ynm1+kyHalf));
      kxHalf = k2x * 0.5;
      kyHalf = k2y * 0.5;
      //k3x = hTimesS * (ynm1 + kyHalf - xnm1 - kxHalf);
      k3x = h * (ynm1 + kyHalf);
      //k3y = h * ((xnm1 + kxHalf) * (r - znm1 - kzHalf) - (ynm1 + kyHalf));
      k3y = h * (-agg - bgg*(xnm1+kxHalf) - gg*(xnm1+kxHalf)*(xnm1+kxHalf)*(xnm1+kxHalf)
		 - g*(xnm1+kxHalf)*(xnm1+kxHalf)*(ynm1+kyHalf)
		 + gg*(xnm1+kxHalf)*(xnm1+kxHalf) - g*(xnm1+kxHalf)*(ynm1+kyHalf));
      //k4x = hTimesS * (ynm1 + k3y - xnm1 - k3x);
      k4x = h * (ynm1 + k3y);
      //k4y = h * ((xnm1 + k3x) * (r - znm1 - k3z) - (ynm1 + k3y));
      k4y = h * (-agg - bgg*(xnm1+k3x) - gg*(xnm1+k3x)*(xnm1+k3x)*(xnm1+k3x)
		 - g*(xnm1+k3x)*(xnm1+k3x)*(ynm1+k3y)
		 + gg*(xnm1+k3x)*(xnm1+k3x) - g*(xnm1+k3x)*(ynm1+k3y));
      xn = xn + (k1x + 2.0*(k2x + k3x) + k4x) * ONESIXTH;
      yn = yn + (k1y + 2.0*(k2y + k3y) + k4y) * ONESIXTH;
      // Euler's method
      // xn = xnm1 + h * (s * (ynm1 - xnm1));
      // yn = ynm1 + h * (xnm1 * (r - znm1) - ynm1);
      // zn = znm1 + h * (xnm1 * ynm1 - b * znm1);
      dx = xn - xnm1;
    }
    counter++;
    ZXP(out) = (xnm1 + dx * frac) * 0.04f;
    frac += slope;
  }
  unit->xn = xn;
  unit->yn = yn;
  unit->counter = counter;
  unit->xnm1 = xnm1;
  unit->ynm1 = ynm1;
  unit->frac = frac;
}
void SyrinxL_Ctor(SyrinxL* unit){
  SETCALC(SyrinxL_next);
  unit->x0 = unit->xn = unit->xnm1 = ZIN0(5);
  unit->y0 = unit->yn = unit->ynm1 = ZIN0(6);
  unit->counter = 0.f;
  unit->frac = 0.f;
  SyrinxL_next(unit, 1);
}
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
PluginLoad(SyrinxUGens)
{
  ft = inTable;
  DefineSimpleUnit(SyrinxL);
}
