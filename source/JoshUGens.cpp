/*
 *  JoshUGens.cpp
 *  xSC3plugins
 *
 *  Created by Josh Parmenter on 2/4/05.
 *  Copyright 2005 __MyCompanyName__. All rights reserved.
 *
 */
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
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program; if not, write to the Free Software
    Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
*/


#include "SC_PlugIn.h"

#define ATS_BAND_WIDTH { 100.0, 100.0, 100.0, 100.0, 110.0, 120.0, 140.0, 150.0, 160.0, 190.0, 210.0, 240.0, 280.0, 320.0, 380.0, 450.0, 550.0, 700.0, 900.0, 1100.0, 1300.0, 1800.0, 2500.0, 3500.0, 4500.0}
#define ATS_CENTER_FREQ {50.0, 150.0, 250.0, 350.0, 455.0, 570.0, 700.0, 845.0, 1000.0, 1175.0, 1375.0, 1600.0, 1860.0, 2160.0, 2510.0, 2925.0, 3425.0, 4050.0, 4850.0, 5850.0, 7050.0, 8600.0, 10750.0, 13750.0, 17750.0} 

// band edges from atsa
//{ 0.0, 100.0, 200.0, 300.0, 400.0, 510.0, 630.0, 770.0, 920.0,  1080.0, 1270.0, 1480.0, 1720.0, 2000.0, 2320.0, 2700.0, 3150.0, 3700.0, 4400.0, 5300.0, 6400.0, 7700.0, 9500.0, 12000.0, 15500.0, 20000.0}

// macros to put rgen state in registers
#define RGET \
	RGen& rgen = *unit->mParent->mRGen; \
	uint32 s1 = rgen.s1; \
	uint32 s2 = rgen.s2; \
	uint32 s3 = rgen.s3; 

#define RPUT \
	rgen.s1 = s1; \
	rgen.s2 = s2; \
	rgen.s3 = s3;
	
static InterfaceTable *ft;

struct A2B : public Unit
{
};

struct B2A : public Unit
{
};

struct B2UHJ : public Unit
{
	float m_wy1[12];
	float m_xy1[12];
	float m_yy1[12];
	float m_coefs[12];
};

struct UHJ2B : public Unit
{
	float m_lsy1[12];
	float m_rsy1[12];
	float m_coefs[12];

};

struct LPCVals : public Unit
{
	SndBuf *m_buf;
	float m_fbufnum, m_framestart;
};

struct LPCSynth : public Unit
{
	SndBuf *m_buf;
	float m_fbufnum;
	int m_valindex; 
	float m_framestart;
	float m_storevals[100]; // this will handle up to 50 poles


};
struct AtsSynth : public Unit
{
	int32 *m_phase; 
	int32 m_phaseoffset, m_lomask;
	float m_numPartials, m_partialStart, m_partialSkip;
	float m_sinePct, m_noisePct, m_freqShift;
	float m_fbufnum, m_framestart;
	double m_cpstoinc, m_radtoinc;
	SndBuf *m_buf;
};

struct PVSynth : public Unit
{
	int32 *m_phase; 
	int32 m_phaseoffset, m_lomask;
	float m_numPartials, m_partialStart, m_partialSkip;
	float m_sinePct, m_noisePct, m_freqShift;
	float m_fbufnum;
	double m_cpstoinc, m_radtoinc;
	SndBuf *m_buf;
};

struct PVInfo : public Unit
{
	float m_partialNum;
	float m_numPartials, m_partialStart, m_partialSkip;
	float m_fbufnum, m_framestart;
	SndBuf *m_buf;
};

struct AtsNoiSynth : public Unit
{
	int32 *m_phase; 
	int32  m_lomask;
	float m_numPartials, m_partialStart, m_partialSkip;
	float m_numBands, m_bandStart, m_bandSkip;
	float m_sinePct, m_noisePct, m_freqShift;
	float m_fbufnum,  m_framestart;
	double m_cpstoinc, m_radtoinc;
	// for noise resynth
	int32 m_sinphase[25];
	int32 m_sinphaseoffset, m_sinlomask;
	double m_sincpstoinc, m_sinradtoinc;
	float m_atsBandFreq[25]; 
	float m_atsCenterFreq[25];
	float m_level[25];
	float m_slope[25];
	int32 m_counter[25];
	uint32 m_atss1[25];
	uint32 m_atss2[25];
	uint32 m_atss3[25];
	SndBuf *m_buf;	
};

struct AtsFreq : public Unit
{
	float m_partialNum;
	float m_numPartials, m_partialStart, m_partialSkip;
	float m_fbufnum, m_framestart;
	SndBuf *m_buf;
};

struct AtsAmp : public Unit
{
	float m_partialNum;
	float m_numPartials, m_partialStart, m_partialSkip, m_framestart;
	float m_fbufnum;
	SndBuf *m_buf;
};

struct AtsNoise : public Unit
{
	float m_bandNum;
	float m_numPartials, m_partialStart, m_partialSkip, m_framestart;
	float m_fbufnum;
	SndBuf *m_buf;
};

struct AtsParInfo : public Unit
{
	float m_partialNum;
	float m_numPartials, m_partialStart, m_partialSkip;
	float m_fbufnum, m_framestart;
	SndBuf *m_buf;
};

struct SinTone : public Unit
{
	double m_phase;
};

struct BFEncode1 : public Unit
{
	float m_azimuth, m_elevation, m_rho, m_level, m_W_amp, m_X_amp, m_Y_amp, m_Z_amp;
};

struct BFEncode2 : public Unit
{
	float m_point_x, m_point_y, m_elevation, m_level, m_W_amp, m_X_amp, m_Y_amp, m_Z_amp;
}; 
/*
struct BFEncodeSter : public Unit
{
	float m_azimuthl, m_azimuthr, m_elevation, m_rho, m_level, m_W_amp, m_X_amp, m_Y_amp, m_Z_amp;
}; 
*/
struct BFDecode1 : public Unit
{
	float m_cosa, m_sina, m_sinb;
};

struct BFManipulate : public Unit
{
	float m_rotate, m_tilt, m_tumble;
};

struct B2Ster : public Unit
{
};

struct Maxamp : public Unit
{
	float m_maxamp, m_numSamps, m_remainingSamps, m_time, m_totalMaxamp, m_totalMaxampTime;
};

struct Grain
{
	double phase, rate;
	double b1, y1, y2; // envelope
	int counter;
	int bufnum;
	int interp;
};

const int kMaxGrains = 64;

struct Warp1 : public Unit
{
	float mPrevTrig;
	int mNumActive;
	uint32 ms1, ms2, ms3;
	int mNextGrain;
	Grain mGrains[kMaxGrains];
};

struct MonoGrain : public Unit
{
	int mNumActive;
	uint32 ms1, ms2, ms3;
	int m_NextGrain;
	Grain mGrains[kMaxGrains];
};

struct Grainbf
{
	double phase, rate;
	double b1, y1, y2; // envelope
	int counter;
	int bufnum;
	int interp;
	double loc, elevation, rho;
	float m_azimuth, m_elevation, m_rho, m_level, m_W_amp, m_X_amp, m_Y_amp, m_Z_amp;

};

struct MonoGrainBF : public Unit
{
	int mNumActive;
	uint32 ms1, ms2, ms3;
	int m_NextGrain;
	Grainbf mGrains[kMaxGrains];
};

struct AudioMSG : public Unit
{
	float m_index;
};

struct Balance : public Unit
{
	float m_scaler, m_pastTestSig;
	float *m_pastTestSums;
	int m_numTurns, m_sampbuf, m_numPeriods;
};

struct MoogVCF : public Unit
{
	float m_fco;
	float m_res;
	float m_xnm1, m_y1nm1, m_y2nm1, m_y3nm1;
	float m_y1n, m_y2n, m_y3n, m_y4n;
};

struct PosRatio : public Unit
{

	int maxsamples, posvals, numvals, testposvals, testnumvals;
	float period, lastperiod, ratio, testratio, outval, lastsample;
};

struct DelayUnit : public Unit
{
	float *m_dlybuf;

	float m_dsamp, m_fdelaylen;
	float m_delaytime, m_maxdelaytime;
	long m_iwrphase, m_idelaylen, m_mask;
	long m_numoutput;
};

struct FeedbackDelay : public DelayUnit
{
	float m_feedbk, m_decaytime;
};

struct Pluck : public FeedbackDelay
{	float m_lastsamp, m_prevtrig, m_coef;
	long m_inputsamps;
};

struct CombLP : public FeedbackDelay
{	float m_lastsamp, m_prevtrig, m_coef;
	long m_inputsamps;
};

extern "C"
{
	void load(InterfaceTable *inTable);

	void A2B_next(A2B *unit, int inNumSamples);
	void A2B_Ctor(A2B* unit);
	
	// B2A and A2B are the same transformation... B2A will call A2B_next
	void B2A_Ctor(A2B* unit);

	void B2UHJ_next(B2UHJ *unit, int inNumSamples);
	void B2UHJ_Ctor(B2UHJ* unit);
	
	void UHJ2B_next(UHJ2B *unit, int inNumSamples);
	void UHJ2B_Ctor(UHJ2B* unit);
	
	void Maxamp_next(Maxamp *unit, int inNumSamples);
	void Maxamp_Ctor(Maxamp* unit);
	
	// being a little lazy for the time being... treat PVSynth as
	// a type 1 ATS file
	
	void PVSynth_Ctor(AtsSynth* unit);
	void PVSynth_Dtor(AtsSynth* unit);
	
	void PVInfo_next(AtsParInfo *unit, int inNumSamples);
	void PVInfo_Ctor(AtsParInfo* unit);
	
	void AtsSynth_next(AtsSynth *unit, int inNumSamples);
	void AtsSynth_Ctor(AtsSynth* unit);
	void AtsSynth_Dtor(AtsSynth* unit);
	
	void AtsNoiSynth_next(AtsNoiSynth *unit, int inNumSamples);
	void AtsNoiSynth_Ctor(AtsNoiSynth* unit);
	void AtsNoiSynth_Dtor(AtsNoiSynth* unit);
	
	void AtsFreq_next(AtsFreq *unit, int inNumSamples);
	void AtsFreq_Ctor(AtsFreq* unit);

	void AtsAmp_next(AtsAmp *unit, int inNumSamples);
	void AtsAmp_Ctor(AtsAmp* unit);

	void AtsNoise_next(AtsNoise *unit, int inNumSamples);
	void AtsNoise_Ctor(AtsNoise* unit);
	
	void AtsParInfo_next(AtsParInfo *unit, int inNumSamples);
	void AtsParInfo_Ctor(AtsParInfo* unit);
			    
	void SinTone_next(SinTone *unit, int inNumSamples);
	void SinTone_Ctor(SinTone* unit);
	
	void BFEncode1_next(BFEncode1 *unit, int inNumSamples);
	void BFEncode1_Ctor(BFEncode1* unit);

	void BFEncode2_next(BFEncode2 *unit, int inNumSamples);
	void BFEncode2_Ctor(BFEncode2* unit); 
/*
	void BFEncodeSter_next(BFEncodeSter *unit, int inNumSamples);
	void BFEncodeSter_Ctor(BFEncodeSter* unit);
*/
	void BFDecode1_next(BFDecode1 *unit, int inNumSamples);
	void BFDecode1_Ctor(BFDecode1* unit);

	void BFManipulate_next(BFManipulate *unit, int inNumSamples);
	void BFManipulate_Ctor(BFManipulate* unit);

	void B2Ster_next(B2Ster *unit, int inNumSamples);
	void B2Ster_Ctor(B2Ster* unit);
	
	void Warp1_next(Warp1 *unit, int inNumSamples);
	void Warp1_Ctor(Warp1* unit);

	void MonoGrain_next(MonoGrain *unit, int inNumSamples);
	void MonoGrain_Ctor(MonoGrain* unit);

	void MonoGrainBF_next(MonoGrainBF *unit, int inNumSamples);
	void MonoGrainBF_Ctor(MonoGrainBF* unit);
			
	void LPCVals_next_k(LPCVals *unit,  int inNumSamples);
	void LPCVals_next_a(LPCVals *unit,  int inNumSamples);
	void LPCVals_Ctor(LPCVals* unit);

	void LPCSynth_next_k(LPCSynth *unit, int inNumSamples);
	void LPCSynth_next_a(LPCSynth *unit, int inNumSamples);
	void LPCSynth_Ctor(LPCSynth* unit);
	
	void AudioMSG_Ctor(AudioMSG* unit);
	void AudioMSG_next_a(AudioMSG *unit, int inNumSamples);
	void AudioMSG_next_k(AudioMSG *unit, int inNumSamples);
	
	void Balance_Ctor(Balance* unit);
	void Balance_Dtor(Balance* unit);
	void Balance_next_a(Balance *unit, int inNumSamples);
	void Balance_next_za(Balance *unit, int inNumSamples);
	void Balance_next_k(Balance *unit, int inNumSamples);
	void Balance_next_zk(Balance *unit, int inNumSamples);
		
	void MoogVCF_Ctor(MoogVCF* unit);
	void MoogVCF_next_kk(MoogVCF *unit, int inNumSamples);
	void MoogVCF_next_ka(MoogVCF *unit, int inNumSamples);
	void MoogVCF_next_ak(MoogVCF *unit, int inNumSamples);
	void MoogVCF_next_aa(MoogVCF *unit, int inNumSamples);
	
	void PosRatio_Ctor(PosRatio* unit);
	void PosRatio_next(PosRatio *unit, int inNumSamples);

	void Pluck_Ctor(Pluck* unit);
	void Pluck_next_aa(Pluck *unit, int inNumSamples);
	void Pluck_next_aa_z(Pluck *unit, int inNumSamples);
	void Pluck_next_kk(Pluck *unit, int inNumSamples);
	void Pluck_next_kk_z(Pluck *unit, int inNumSamples);
	void Pluck_next_ka(Pluck *unit, int inNumSamples);
	void Pluck_next_ka_z(Pluck *unit, int inNumSamples);
	void Pluck_next_ak(Pluck *unit, int inNumSamples);
	void Pluck_next_ak_z(Pluck *unit, int inNumSamples);

	void CombLP_Ctor(CombLP* unit);
	void CombLP_next_aa(CombLP *unit, int inNumSamples);
	void CombLP_next_aa_z(CombLP *unit, int inNumSamples);
	void CombLP_next_kk(CombLP *unit, int inNumSamples);
	void CombLP_next_kk_z(CombLP *unit, int inNumSamples);
	void CombLP_next_ka(CombLP *unit, int inNumSamples);
	void CombLP_next_ka_z(CombLP *unit, int inNumSamples);
	void CombLP_next_ak(CombLP *unit, int inNumSamples);
	void CombLP_next_ak_z(CombLP *unit, int inNumSamples);
}


//////////////////////////////////////////////////////////////////////////////////////////////////

void AudioMSG_Ctor(AudioMSG* unit)
{
    if (INRATE(1) == calc_FullRate) {
	SETCALC(AudioMSG_next_a);
	AudioMSG_next_a(unit, 1);
	   } else {
	SETCALC(AudioMSG_next_k);
	unit->m_index = IN0(1);
	AudioMSG_next_k(unit, 1);
	 };
}

void AudioMSG_next_a(AudioMSG *unit, int inNumSamples)
{
    float* in = IN(0);
    float* index = IN(1);
    float* out = OUT(0);

    for (int i = 0; i < inNumSamples; ++i)
	{ 
	    out[i] = sin((in[i] * pi2) + (index[i] * sin(in[i] * twopi)));
	}
}

void AudioMSG_next_k(AudioMSG *unit, int inNumSamples)
{
    float* in = IN(0);
    float endindex = IN0(1);
    float* out = OUT(0);

    float index = unit->m_index;

    float indexinc = CALCSLOPE(endindex, index);

    for (int i = 0; i < inNumSamples; ++i) {
	    out[i] = sin((in[i] * pi2) + (index * sin(in[i] * twopi)));
	    index += indexinc;
	}

    unit->m_index = endindex;
	
}

void Maxamp_Ctor(Maxamp* unit)
{
	SETCALC(Maxamp_next);
	unit->m_numSamps = unit->m_remainingSamps = ZIN0(1);
	unit->m_maxamp = 0.0;
	unit->m_time = 0.0;
	unit->m_totalMaxampTime = 0.0;
	unit->m_totalMaxamp = 0.0;
	Maxamp_next(unit, 1);
}

void Maxamp_next(Maxamp *unit, int inNumSamples)
{
	float* in = IN(0);
	float* out = OUT(0);
	int numSamps = (int)unit->m_numSamps;
	int remainingSamps = (int)unit->m_remainingSamps;
	float currenttime = unit->m_time;
	float totalMaxamp = unit->m_totalMaxamp;
	float maxamp = unit->m_maxamp;
	float totalMaxampTime = unit->m_totalMaxampTime;
	float curamp, lasttime;
	
	for (int j = 0; j < inNumSamples; ++j){
	    curamp = fabsf(in[j]);
	    if (curamp > maxamp) maxamp = unit->m_maxamp = curamp;
	    if (curamp > totalMaxamp){ 
		totalMaxamp = unit->m_totalMaxamp = curamp;
		totalMaxampTime = unit->m_totalMaxampTime = currenttime + ((numSamps - remainingSamps) * SAMPLEDUR);
		}
	    if (remainingSamps == 0){
		lasttime = currenttime;
		currenttime = unit->m_time = currenttime + (SAMPLEDUR * numSamps);     
		Print("Maxamp from %g to %g: %g, Overall Maxamp at %g: %g\n", lasttime, currenttime, maxamp, totalMaxampTime, totalMaxamp);
		unit->m_remainingSamps = remainingSamps = numSamps;
		unit->m_maxamp = maxamp = 0.0;
		}
		else{
		    unit->m_remainingSamps = --remainingSamps;
		    }
	    out[j] =  in[j];
	}   
}
	
void AtsSynth_Ctor(AtsSynth* unit)
{
	SETCALC(AtsSynth_next);
	int tableSizeSin = ft->mSineSize;
	unit->m_numPartials = ZIN0(1);
	unit->m_partialStart = ZIN0(2);
	unit->m_partialSkip = ZIN0(3);
	unit->m_framestart = ZIN0(4);
	unit->m_lomask = (tableSizeSin - 1) << 3; 
	unit->m_radtoinc = tableSizeSin * (rtwopi * 65536.); 
	unit->m_cpstoinc = tableSizeSin * SAMPLEDUR * 65536.;     
	unit->m_phase = (int32*)RTAlloc(unit->mWorld, (int)unit->m_numPartials * sizeof(int32));
	unit->m_fbufnum = -1e9f;
	ClearUnitOutputs(unit, 1);
}

void AtsSynth_Dtor(AtsSynth* unit)
{	
	RTFree(unit->mWorld, unit->m_phase);
}


void AtsSynth_next(AtsSynth *unit, int inNumSamples)
{
	float fbufnum = IN0(0);
	if (fbufnum != unit->m_fbufnum) {
		uint32 bufnum = (int)fbufnum;
		World *world = unit->mWorld;
		if (bufnum >= world->mNumSndBufs) bufnum = 0;
		unit->m_fbufnum = fbufnum;
		unit->m_buf = world->mSndBufs + bufnum;
	}
	SndBuf *buf = unit->m_buf;
	float *bufData __attribute__((__unused__)) = buf->data; 
	float *out = OUT(0);
	float freqMul = IN0(5);
	float freqAdd = IN0(6);
	float* filePartials = bufData + 1;
	float* fileFrames = bufData + 2;
	int ampDataStart = 4;
	int freqDataStart = ((int)filePartials[0] * (int)fileFrames[0]) + 4;
	//leaving phase out for now
//	int phaseDataStart = (((int)filePartials[0] * (int)fileFrames[0]) * 2) + 4;
	int numFrames = (int)fileFrames[0];
	float framestart = unit->m_framestart * (numFrames - 1);
	float frame = IN0(4);
	float frameend = frame * (numFrames - 1);
	// get the frame vals for the freq data
	int frameround1 = (int)frameend;
	int frameround2 = 1;
	// get frames vals for amp data (interpolated over the k-rate)
	int startframeround1 = (int)framestart;
	int startframeround2 = 1;
	
	float framepct = frameend - (float)frameround1;
	float startframepct = framestart - (float)startframeround1;
	
	float numPartials = unit->m_numPartials;
	float partialStart = unit->m_partialStart;
	float partialSkip = unit->m_partialSkip;
	int totalPartials = (int)numPartials;

	int *partials = new int[(int)numPartials]; //this will be the partials to synthesize
	
	startframeround2 = 1;

	if ((frameround1 + 1) > numFrames) frameround2 = 0;	

	if ((startframeround1 + 1) > numFrames) startframeround2 = 0;

	if (!bufData) { 
                unit->mDone = true; 
		ClearUnitOutputs(unit, inNumSamples); 
		return;
		};
	
	//fill partials with the correct partials, check if partials requested are actually present
	
	for (int j = 0; j < (int)numPartials; ++j){
	    int aPartial = (int)partialStart + ((int)partialSkip * j);
	    if (aPartial < (int)filePartials[0]) 
		partials[j] = aPartial;
		else
		--totalPartials;
	}    

	ClearUnitOutputs(unit, inNumSamples);

	float *table0 = ft->mSineWavetable;
	float *table1 = table0 + 1;
	int32 lomask = unit->m_lomask;
	
	int ampData, freqData;
	float startampin, endampin, ampslope, freqin;
	int32 freq, phase;
	float* startatsAmp1;
	float* endatsAmp1;
	float* atsFreq1;
	
	// this loop calulates the partial sample values for the each of the sins, stores them in an array
	for (int i = 0; i< totalPartials; ++i){
	    ampData = ampDataStart + (numFrames * partials[i]);
	    freqData = freqDataStart + (numFrames * partials[i]);
	    
	    startatsAmp1 = bufData + (int)(ampData + startframeround1); //this is the position in the file of the amp

	    endatsAmp1 = bufData + (int)(ampData + frameround1); 
	    	    
	    startampin = lininterp(startframepct, startatsAmp1[0], startatsAmp1[startframeround2]); 
	    endampin = lininterp(framepct, endatsAmp1[0], endatsAmp1[frameround2]); 
	    ampslope = CALCSLOPE(endampin, startampin);
	    
	    atsFreq1 = bufData + (int)(freqData + frameround1); //this is the position in the file of the freq

	    freqin = (lininterp(framepct, atsFreq1[0], atsFreq1[frameround2]) * freqMul) + freqAdd; 
	    freq = (int32)(unit->m_cpstoinc * freqin);
	    
	    phase = unit->m_phase[i];
	
	    for (int j = 0; j < inNumSamples; ++j){
		out[j] += (lookupi1(table0, table1, phase, lomask) * startampin);
		phase += freq;
		startampin += ampslope;
	    };
	    unit->m_phase[i] = phase;
	    }
	
	unit->m_framestart = frame;

	delete [] partials;
}


/////////////////////////////////////////////////////////////////

void AtsNoiSynth_Ctor(AtsNoiSynth* unit)
{
	SETCALC(AtsNoiSynth_next);
	int i;
	int tableSizeSin = ft->mSineSize;
	unit->m_numPartials = ZIN0(1);
	unit->m_partialStart = ZIN0(2);
	unit->m_partialSkip = ZIN0(3);
	unit->m_framestart = ZIN0(4);
	unit->m_numBands = ZIN0(9);
	unit->m_bandStart = ZIN0(10);
	unit->m_bandSkip = ZIN0(11);
	int numPartials = (int)unit->m_numPartials;
	//int numBands = (int)unit->m_numBands;
	unit->m_lomask = (tableSizeSin - 1) << 3; 
	unit->m_radtoinc = tableSizeSin * (rtwopi * 65536.); 
	unit->m_cpstoinc = tableSizeSin * SAMPLEDUR * 65536.; 
	
	// this allocates space to save phase data for each partial... 
	
	unit->m_phase = (int32*)RTAlloc(unit->mWorld, numPartials * sizeof(int32));
	
	unit->m_fbufnum = -1e9f;

	float atsBandFreq[25] = ATS_BAND_WIDTH;
	float atsCenterFreq[25] = ATS_CENTER_FREQ;

	for (i = 0; i < 25; i++){
	    unit->m_atsBandFreq[i] = atsBandFreq[i];
	    unit->m_atsCenterFreq[i] = atsCenterFreq[i];
	    unit->m_counter[i] = 0;
	    unit->m_level[i] = unit->mParent->mRGen->frand2();
	    unit->m_slope[i] = 0.f;
	    }
		
	unit->m_sinlomask = (tableSizeSin - 1) << 3; 
	unit->m_sinradtoinc = tableSizeSin * (rtwopi * 65536.); 
	unit->m_sincpstoinc = tableSizeSin * SAMPLEDUR * 65536.; 
	
	ClearUnitOutputs(unit, 1);

}

void AtsNoiSynth_Dtor(AtsNoiSynth* unit)
{	
	RTFree(unit->mWorld, unit->m_phase);
}

void AtsNoiSynth_next(AtsNoiSynth *unit, int inNumSamples)
{
	float fbufnum  = ZIN0(0);
	if (fbufnum != unit->m_fbufnum) {
		uint32 bufnum = (int)fbufnum;
		World *world = unit->mWorld;
		if (bufnum >= world->mNumSndBufs) bufnum = 0;
		unit->m_fbufnum = fbufnum;
		unit->m_buf = world->mSndBufs + bufnum;
	}
	SndBuf *buf = unit->m_buf;
	float *bufData __attribute__((__unused__)) = buf->data; 
	float *out = OUT(0);
//	float* filetype = bufData;
	float* filePartials = bufData + 1;
	float* fileFrames = bufData + 2;
	float* winSize = bufData + 3;
	int ampDataStart = 4;
	int freqDataStart = ((int)filePartials[0] * (int)fileFrames[0]) + 4;
	// leave phase out for now
	//int phaseDataStart = (((int)filePartials[0] * (int)fileFrames[0]) * 2) + 4;
	int noiseDataStart = (((int)filePartials[0] * (int)fileFrames[0]) * 3) + 4;
	int numFrames = (int)fileFrames[0];
	float framestart = unit->m_framestart * (numFrames - 1);
	float frame = IN0(4);
	float frameend = frame * (numFrames - 1);
	// get the frame vals for the freq data
	int frameround1 = (int)frameend;
	int frameround2 = 1;
	// get frames vals for amp data (interpolated over the k-rate)
	int startframeround1 = (int)framestart;
	int startframeround2 = 1;
	
	float framepct = frameend - (float)frameround1;
	float startframepct = framestart - (float)startframeround1;
	
	float numPartials = unit->m_numPartials;
	float partialStart = unit->m_partialStart;
	float partialSkip = unit->m_partialSkip;
	float numBands = unit->m_numBands;
	float bandStart = unit->m_bandStart;
	float bandSkip = unit->m_bandSkip;
	float sinePct = IN0(5);
	float noisePct = IN0(6);
	float freqMul = IN0(7);
	float freqAdd = IN0(8);
	int totalPartials = (int)numPartials;
	int totalBands = (int)numBands;
	
	int *partials = new int[totalPartials]; //this will be the partials to synthesize
	int *bands = new int[totalBands]; //this will be the bands to synthesize
	
	//frameround2 = frameround1 + 1;
	//startframeround2 = startframeround1 + 1;

	if ((frameround1 + 1) > numFrames) frameround2 = 0;	
	if ((startframeround1 + 1) > numFrames) startframeround2 = 0;
	
	if (!bufData) { 
                unit->mDone = true; 
		ClearUnitOutputs(unit, inNumSamples); 
		return;
		};
		
	//fill partials with the correct partials, check if partials requested are actually present
	
	for (int i = 0; i < (int)numPartials; ++i){
	    int aPartial = (int)partialStart + ((int)partialSkip * i);
	    if (aPartial < (int)filePartials[0]) 
		partials[i] = aPartial;
		else
		--totalPartials;
	}    
	
    	ClearUnitOutputs(unit, inNumSamples);

	float *table0 = ft->mSineWavetable;
	float *table1 = table0 + 1;
	int32 lomask = unit->m_lomask;
	    
	// this loop calulates the partial sample values for the sins, stores them in an array
	for (int i = 0; i< totalPartials; ++i){
	    
	    int ampData = ampDataStart + (numFrames * partials[i]);
	    int freqData = freqDataStart + (numFrames * partials[i]);
	    
	    float* startatsAmp1 = bufData + (int)(ampData + startframeround1); //this is the position in the file of the amp

	    float* endatsAmp1 = bufData + (int)(ampData + frameround1); 
	    	    
	    float startampin = lininterp(startframepct, startatsAmp1[0], startatsAmp1[startframeround2]) * sinePct;
	    float endampin = lininterp(framepct, endatsAmp1[0], endatsAmp1[frameround2]) * sinePct;
	    float ampslope = CALCSLOPE(endampin, startampin);
	    
	    float* atsFreq1 = bufData + (int)(freqData + frameround1); //this is the position in the file of the freq

	    float freqin = ((lininterp(framepct, atsFreq1[0], atsFreq1[frameround2]) * freqMul) + freqAdd);
	    int32 freq = (int32)(unit->m_cpstoinc * freqin);
	    
	    int32 phase = unit->m_phase[i];
	
	    for (int j = 0; j < inNumSamples; ++j){
		out[j] += (lookupi1(table0, table1, phase, lomask) * startampin);
		phase += freq;
		startampin += ampslope;
	    };
	    unit->m_phase[i] = phase;
	    }
	
	for (int i = 0; i < (int)numBands; ++i){
	    int aBand = (int)bandStart + ((int)bandSkip * i);
	    if (aBand < 25) 
		bands[i] = aBand;
		else
		--totalBands;
	}    
	
	float scale = (winSize[0] * 0.33166610955984); 
	
	//calculate the noise
	for (int i = 0; i < totalBands; i++){
	    int bandnum = bands[i];
		
	    int noiseData = noiseDataStart + ((int)numFrames * bandnum);
	    float* startatsnoise1 = bufData + (int)(noiseData + startframeround1); //this is the position in the file of the noise
	    float* endatsnoise1 = bufData + (int)(noiseData + frameround1); 
	    float startnoisein = sqrt((lininterp(startframepct, startatsnoise1[0], startatsnoise1[startframeround2]) * noisePct) / scale);
	    float endnoisein = sqrt((lininterp(framepct, endatsnoise1[0], endatsnoise1[frameround2]) * noisePct) / scale);
	    float noiseslope = CALCSLOPE(endnoisein, startnoisein); 
	    
	    int32 sinphase = unit->m_sinphase[i];
		
	    float noisefreq = unit->m_atsBandFreq[bandnum];
	    float sinfreq = unit->m_atsCenterFreq[bandnum];
	    int32 freq = (int32)(unit->m_sincpstoinc * sinfreq); 

	    float level = unit->m_level[i];
	    float slope = unit->m_slope[i];
	    int32 counter = unit->m_counter[i];
		
	    RGET
		
	    for (int j = 0; j < inNumSamples; j++){
		if (counter<=0) {
		    counter = (int32)(unit->mRate->mSampleRate / sc_max(noisefreq, .001f));
		    counter = sc_max(1, counter);
		    float nextlevel = frand2(s1, s2, s3);
		    slope = (nextlevel - level) / counter;
		} else --counter;
		out[j] += lookupi1(table0, table1, sinphase, lomask) * startnoisein * level;
		startnoisein += noiseslope;
		level += slope;
    		sinphase += freq;
		}
		
	    RPUT
	    unit->m_level[i] = level;
	    unit->m_slope[i] = slope;
	    unit->m_counter[i] = counter;
	    unit->m_sinphase[i] = sinphase;
	    }
	
    	unit->m_framestart = frame;
		delete [] partials;
		delete [] bands;
}

/////////////////////////////////////////////////////////////////

void AtsFreq_Ctor(AtsFreq* unit)
{
	SETCALC(AtsFreq_next);
	unit->m_partialNum = ZIN0(1);
	unit->m_fbufnum = -1e9f;
	unit->m_framestart = ZIN0(2);
	ClearUnitOutputs(unit, 1);

}

void AtsFreq_next(AtsFreq *unit, int inNumSamples)
{

	float fbufnum  = ZIN0(0);
	if (fbufnum != unit->m_fbufnum) {
		uint32 bufnum = (int)fbufnum;
		World *world = unit->mWorld;
		if (bufnum >= world->mNumSndBufs) bufnum = 0;
		unit->m_fbufnum = fbufnum;
		unit->m_buf = world->mSndBufs + bufnum;
	    }
	SndBuf *buf = unit->m_buf;
	float *bufData __attribute__((__unused__)) = buf->data; 
	float *out = OUT(0);
	float* filePartials = bufData + 1;
	float* fileFrames = bufData + 2;
	int freqDataStart = ((int)filePartials[0] * (int)fileFrames[0]) + 4;
	float numFrames = fileFrames[0];
	float framestart = unit->m_framestart * ((int)numFrames - 1);
	float frame = IN0(2);
	float frameend = frame * ((int)numFrames - 1);
	int startframeround1 = (int)framestart;
	int startframeround2;
	float startframepct = framestart-(float)startframeround1;
	int endframeround1 = (int)frameend;
	int endframeround2;
	float endframepct = frameend-(float)endframeround1;
		
	float partialNum = unit->m_partialNum;
	
	startframeround2 = startframeround1 + 1;
	endframeround2 = endframeround1 + 1;
	
	if (partialNum >= filePartials[0]){ 
                unit->mDone = true; 
		ClearUnitOutputs(unit, inNumSamples); 
		return;
		};
	    
	if (startframeround2 > numFrames) startframeround2 = startframeround1;
	if (endframeround2 > numFrames) endframeround2 = endframeround1;
	
	if (!bufData) { 
                unit->mDone = true; 
		ClearUnitOutputs(unit, inNumSamples); 
		return;
		};
	
	int freqData = freqDataStart + ((int)numFrames * (int)partialNum);
	int startfreqData1 = freqData + startframeround1;
	int startfreqData2 = freqData  + startframeround2;
	float* startatsFreq1 = bufData + startfreqData1; //this is the position in the file of the freq
	float* startatsFreq2 = bufData + startfreqData2;
	float startfreqin = (startatsFreq1[0] + ((startatsFreq2[0] - startatsFreq1[0]) * startframepct));
	
	int endfreqData1 = freqData + endframeround1;
	int endfreqData2 = freqData + endframeround2;
	float* endatsFreq1 = bufData + endfreqData1; //this is the position in the file of the freq
	float* endatsFreq2 = bufData + endfreqData2;
	float endfreqin = (endatsFreq1[0] + ((endatsFreq2[0] - endatsFreq1[0]) * endframepct));
		
	float freqslope = CALCSLOPE(endfreqin, startfreqin);
	
	for (int j = 0; j < inNumSamples; ++ j){
	    out[j] = startfreqin;
	    startfreqin += freqslope;
	    }
	
	unit->m_framestart = frame;
}
/////////////////////////////////////////////////////////////////

void AtsAmp_Ctor(AtsAmp* unit)
{
	SETCALC(AtsAmp_next);
	unit->m_partialNum = ZIN0(1);
	unit->m_fbufnum = -1e9f;
	ClearUnitOutputs(unit, 1);
	unit->m_framestart = ZIN0(2);

}

void AtsAmp_next(AtsAmp *unit, int inNumSamples)
{

	float fbufnum  = ZIN0(0);
	if (fbufnum != unit->m_fbufnum) {
		uint32 bufnum = (int)fbufnum;
		World *world = unit->mWorld;
		if (bufnum >= world->mNumSndBufs) bufnum = 0;
		unit->m_fbufnum = fbufnum;
		unit->m_buf = world->mSndBufs + bufnum;
	}
	SndBuf *buf = unit->m_buf;
	float *bufData __attribute__((__unused__)) = buf->data; 
	float *out = OUT(0);
	float* filePartials = bufData + 1;
	float* fileFrames = bufData + 2;
	int ampDataStart = 4;
	float numFrames = fileFrames[0];
	float framestart = unit->m_framestart * ((int)numFrames - 1);
	float frame = IN0(2);
	float frameend = frame * ((int)numFrames - 1);
	int startframeround1 = (int)framestart;
	int startframeround2;
	float startframepct = framestart-(float)startframeround1;
	int endframeround1 = (int)frameend;
	int endframeround2;
	float endframepct = frameend - (float)endframeround1;
	float partialNum = unit->m_partialNum;

	
	startframeround2 = startframeround1 + 1;
	endframeround2 = endframeround1 + 1;
	
	if (partialNum >= filePartials[0]){ 
		unit->mDone = true; 
		ClearUnitOutputs(unit, inNumSamples); 
		return;
		};
		
	if (startframeround2 > numFrames) startframeround2 = startframeround1;
	if (endframeround2 > numFrames) endframeround2 = endframeround1;
	
	if (!bufData) { 
                unit->mDone = true; 
		ClearUnitOutputs(unit, inNumSamples); 
		return;
		};
		
	int ampData = ampDataStart + ((int)numFrames * (int)partialNum);
	int startampData1 = ampData + startframeround1;
	int startampData2 = ampData + startframeround2;
	float* startatsAmp1 = bufData + startampData1; //this is the position in the file of the amp
	float* startatsAmp2 = bufData + startampData2;
	float startampin = (startatsAmp1[0] + ((startatsAmp2[0]- startatsAmp1[0]) * startframepct));

	int endampData1 = ampData + endframeround1;
	int endampData2 = ampData + endframeround2;
	float* endatsAmp1 = bufData + endampData1; //this is the position in the file of the amp
	float* endatsAmp2 = bufData + endampData2;
	float endampin = (endatsAmp1[0] + ((endatsAmp2[0]- endatsAmp1[0]) * endframepct));
		   
	float ampin = startampin;
	float ampslope = CALCSLOPE(endampin, startampin);
	for (int j = 0; j < inNumSamples; ++ j){
	    out[j] = ampin;
	    ampin += ampslope;
	    }
	
	unit->m_framestart = frame;
}


void AtsNoise_Ctor(AtsNoise* unit)
{
	SETCALC(AtsNoise_next);
	unit->m_bandNum = IN0(1);
	unit->m_fbufnum = -1e9f;
	ClearUnitOutputs(unit, 1);
	unit->m_framestart = IN0(2);

}

void AtsNoise_next(AtsNoise *unit, int inNumSamples)
{

	float fbufnum  = ZIN0(0);
	if (fbufnum != unit->m_fbufnum) {
		uint32 bufnum = (int)fbufnum;
		World *world = unit->mWorld;
		if (bufnum >= world->mNumSndBufs) bufnum = 0;
		unit->m_fbufnum = fbufnum;
		unit->m_buf = world->mSndBufs + bufnum;
	}
	SndBuf *buf = unit->m_buf;
	float *bufData __attribute__((__unused__)) = buf->data; 
	float *out = OUT(0);
	float* filePartials = bufData + 1;
	float* fileFrames = bufData + 2;
	int noiseDataStart = (((int)filePartials[0] * (int)fileFrames[0]) * 3) + 4;
	float numFrames = fileFrames[0];
	float framestart = unit->m_framestart * ((int)numFrames - 1);
	float frame = IN0(2);
	float frameend = frame * ((int)numFrames - 1);
	int startframeround1 = (int)framestart;
	int startframeround2;
	float startframepct = framestart-(float)startframeround1;
	int endframeround1 = (int)frameend;
	int endframeround2;
	float endframepct = frameend - (float)endframeround1;
	float bandNum = unit->m_bandNum;

	
	startframeround2 = startframeround1 + 1;
	endframeround2 = endframeround1 + 1;
	
	if (bandNum >= 25){ 
		unit->mDone = true; 
		ClearUnitOutputs(unit, inNumSamples); 
		return;
		};
		
	if (startframeround2 > numFrames) startframeround2 = startframeround1;
	if (endframeround2 > numFrames) endframeround2 = endframeround1;
	
	if (!bufData) { 
                unit->mDone = true; 
		ClearUnitOutputs(unit, inNumSamples); 
		return;
		};
	
	int noiseData = noiseDataStart + ((int)numFrames * (int)bandNum);		
	int startnoiseData1 =  noiseData + startframeround1;
	int startnoiseData2 = noiseData + startframeround2;
	float* startatsnoise1 = bufData + startnoiseData1; //this is the position in the file of the noise
	float* startatsnoise2 = bufData + startnoiseData2;
	float startnoisein = (startatsnoise1[0] + ((startatsnoise2[0]- startatsnoise1[0]) * startframepct));

	int endnoiseData1 = noiseData + endframeround1;
	int endnoiseData2 = noiseData + endframeround2;
	float* endatsnoise1 = bufData + endnoiseData1; //this is the position in the file of the noise
	float* endatsnoise2 = bufData + endnoiseData2;
	float endnoisein = (endatsnoise1[0] + ((endatsnoise2[0]- endatsnoise1[0]) * endframepct));
		   
	float noisein = startnoisein;
	float noiseslope = CALCSLOPE(endnoisein, startnoisein);
	for (int j = 0; j < inNumSamples; ++ j){
	    out[j] = noisein;
	    noisein += noiseslope;
	    }
	
	unit->m_framestart = frame;
}


void AtsParInfo_Ctor(AtsParInfo* unit)
{
	SETCALC(AtsParInfo_next);
	unit->m_partialNum = ZIN0(1);
	unit->m_framestart = ZIN0(2);
	unit->m_fbufnum = -1e9f;
	ClearUnitOutputs(unit, 1);

}

void AtsParInfo_next(AtsParInfo *unit, int inNumSamples)
{

	float fbufnum  = ZIN0(0);
	if (fbufnum != unit->m_fbufnum) {
		uint32 bufnum = (int)fbufnum;
		World *world = unit->mWorld;
		if (bufnum >= world->mNumSndBufs) bufnum = 0;
		unit->m_fbufnum = fbufnum;
		unit->m_buf = world->mSndBufs + bufnum;
	}
	SndBuf *buf = unit->m_buf;
	float *bufData __attribute__((__unused__)) = buf->data; 
	float *out1 = OUT(0);
	float *out2 = OUT(1);
	float* filePartials = bufData + 1;
	float* fileFrames = bufData + 2;
	int ampDataStart = 4;
	int freqDataStart = ((int)filePartials[0] * (int)fileFrames[0]) + 4;
	float numFrames = fileFrames[0];
	float framestart = unit->m_framestart * ((int)numFrames - 1);
	float frame = IN0(2);
	float frameend = frame * ((int)numFrames - 1);
	int startframeround1 = (int)framestart;
	int startframeround2;
	float startframepct = framestart-(float)startframeround1;

	int endframeround1 = (int)frameend;
	int endframeround2;
	float endframepct = frameend-(float)endframeround1;
		
	float partialNum = unit->m_partialNum;

	
	startframeround2 = startframeround1 + 1;
	endframeround2 = endframeround1 + 1;
	
	if (partialNum >= filePartials[0]){ 
		unit->mDone = true; 
		ClearUnitOutputs(unit, inNumSamples); 
		return;
		};
		
	if (startframeround2 > numFrames) startframeround2 = startframeround1;
	if (endframeround2 > numFrames) endframeround2 = endframeround1;
	
	if (!bufData) { 
                unit->mDone = true; 
		ClearUnitOutputs(unit, inNumSamples); 
		return;
		};
		
	int startampData1 = ampDataStart + (((int)numFrames * (int)partialNum) + startframeround1);
	int startampData2 = ampDataStart + (((int)numFrames * (int)partialNum) + startframeround2);
	float* startatsAmp1 = bufData + startampData1; //this is the position in the file of the amp
	float* startatsAmp2 = bufData + startampData2;
	float startampin = (startatsAmp1[0] + ((startatsAmp2[0] - startatsAmp1[0]) * startframepct));

	int endampData1 = ampDataStart + (((int)numFrames * (int)partialNum) + endframeround1);
	int endampData2 = ampDataStart + (((int)numFrames * (int)partialNum) + endframeround2);
	float* endatsAmp1 = bufData + endampData1; //this is the position in the file of the amp
	float* endatsAmp2 = bufData + endampData2;
	float endampin = (endatsAmp1[0] + ((endatsAmp2[0] - endatsAmp1[0]) * endframepct));
	
	float ampslope = CALCSLOPE(endampin, startampin);
	
	int startfreqData1 = freqDataStart + (((int)numFrames * (int)partialNum) + startframeround1);
	int startfreqData2 = freqDataStart + (((int)numFrames * (int)partialNum) + startframeround2);
	float* startatsFreq1 = bufData + startfreqData1; //this is the position in the file of the amp
	float* startatsFreq2 = bufData + startfreqData2;
	float startfreqin = (startatsFreq1[0] + ((startatsFreq2[0] - startatsFreq1[0]) * startframepct));

	int endfreqData1 = freqDataStart + (((int)numFrames * (int)partialNum) + endframeround1);
	int endfreqData2 = freqDataStart + (((int)numFrames * (int)partialNum) + endframeround2);
	float* endatsFreq1 = bufData + endfreqData1; //this is the position in the file of the amp
	float* endatsFreq2 = bufData + endfreqData2;
	float endfreqin = (endatsFreq1[0] + ((endatsFreq2[0]- endatsFreq1[0]) * endframepct));
			   
	float freqslope = CALCSLOPE(endfreqin, startfreqin);
	
	for (int j = 0; j < inNumSamples; ++ j){
	    out1[j] = startampin;
	    out2[j] = startfreqin;
	    startampin += ampslope;
	    startfreqin += freqslope;
	    }
	
	unit->m_framestart = frame;
}


/////////////////////////////////////////////////////////////////

void SinTone_Ctor(SinTone* unit)
{
    SETCALC(SinTone_next);
    unit->m_phase = IN0(1);
    SinTone_next(unit, 1);

}

void SinTone_next(SinTone *unit, int inNumSamples)
{
	float freq = IN0(0);
	double phase = unit->m_phase;
	float *valout = OUT(0);
	float incr = (freq / SAMPLERATE) * 6.2831853071796;
	
	for (int i = 0; i< inNumSamples; ++i)
	    {
	    valout[i] = sin(phase);
	    phase += incr;
	    };
	
	unit->m_phase = phase;
}


//////////////////////////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////////////////////////////////

// A2B and B2A conversions
/* See Ambisonics Ugens (ATK) */
void A2B_Ctor(A2B *unit)
{
	SETCALC(A2B_next);
	A2B_next(unit, 1);
}

void B2A_Ctor(A2B *unit)
{
	SETCALC(A2B_next);
	A2B_next(unit, 1);
}

// four channels in, four channels out... 
// regardless... call in w, x, y, z and out a, b, c, d
void A2B_next(A2B *unit, int inNumSamples)
{       
	float *Aout = ZOUT(0);
	float *Bout = ZOUT(1);
	float *Cout = ZOUT(2);
	float *Dout = ZOUT(3);
	
	float *Win = ZIN(0);
	float *Xin = ZIN(1);
	float *Yin = ZIN(2);
	float *Zin = ZIN(3);

	float lf, rf, rr, lr;
	float oneHalf = 0.5;
	
	LOOP(inNumSamples,
	    lf = ZXP(Win) * oneHalf;
	    rf = ZXP(Xin) * oneHalf;
	    lr = ZXP(Yin) * oneHalf;
	    rr = ZXP(Zin) * oneHalf;
	    
	    ZXP(Aout) = lf + rf + lr + rr;
	    ZXP(Bout) = lf + rf + (lr * -1) + (rr * -1);
	    ZXP(Cout) = lf + (rf * -1) + lr + (rr * -1);
	    ZXP(Dout) = lf + (rf * -1) + (lr * -1) + rr;
	    );
}

///////////////////////////////////////////////////////////////////////////////////////////////////////
// B to UHJ - transform a 3 channel 2d ambisonic signal into a 2 channel UHJ
//
//		L = 0.5 * (0.9397*W + 0.1856*X - j*0.342*W + j*0.5099*X + 0.655*Y)
//		R = 0.5 * (0.9397*W+ 0.1856*X + j*0.342*W - j*0.5099*X - 0.655*Y)
//
// Where j  is a 90 shift

void B2UHJ_Ctor(B2UHJ *unit)
{
    SETCALC(B2UHJ_next);
    float gamconst = (15.0 * pi) / SAMPLERATE;
    float gamma01 = gamconst * 0.3609;
    float gamma02 = gamconst * 2.7412;
    float gamma03 = gamconst * 11.1573;
    float gamma04 = gamconst * 44.7581;
    float gamma05 = gamconst * 179.6242;
    float gamma06 = gamconst * 798.4578;
    float gamma07 = gamconst * 1.2524;
    float gamma08 = gamconst * 5.5671;
    float gamma09 = gamconst * 22.3423;
    float gamma10 = gamconst * 89.6271;
    float gamma11 = gamconst * 364.7914;
    float gamma12 = gamconst * 2770.1114;
    unit->m_coefs[0] = (gamma01 - 1) / (gamma01 + 1);
    unit->m_coefs[1] = (gamma02 - 1) / (gamma02 + 1);
    unit->m_coefs[2] = (gamma03 - 1) / (gamma03 + 1);
    unit->m_coefs[3] = (gamma04 - 1) / (gamma04 + 1);
    unit->m_coefs[4] = (gamma05 - 1) / (gamma05 + 1);
    unit->m_coefs[5] = (gamma06 - 1) / (gamma06 + 1);
    unit->m_coefs[6] = (gamma07 - 1) / (gamma07 + 1);
    unit->m_coefs[7] = (gamma08 - 1) / (gamma08 + 1);
    unit->m_coefs[8] = (gamma09 - 1) / (gamma09 + 1);
    unit->m_coefs[9] = (gamma10 - 1) / (gamma10 + 1);
    unit->m_coefs[10] = (gamma11 - 1) / (gamma11 + 1);
    unit->m_coefs[11] = (gamma12 - 1) / (gamma12 + 1);
    // clear out the storage variables for the filters
    for(int i = 0; i < 12; ++i){
	unit->m_wy1[i] = 0.0;
	unit->m_xy1[i] = 0.0;
	unit->m_yy1[i] = 0.0;
	}
    ClearUnitOutputs(unit, 1);
//    B2UHJ_next(unit, 1);
}

void B2UHJ_next(B2UHJ *unit, int inNumSamples)
{
    float *Win = IN(0);
    float *Xin = IN(1);
    float *Yin = IN(2);
    float *lsout = OUT(0);
    float *rsout = OUT(1);

    float wy1[12];
    float xy1[12];
    float yy1[12];

    // each filter's last sample
    for(int i = 0; i < 12; ++i){
	wy1[i] = unit->m_wy1[i];
	xy1[i] = unit->m_xy1[i];
	yy1[i] = unit->m_yy1[i];
	}

    for (int i = 0; i < inNumSamples; ++i){
	float wsig, w, wj, xsig, x, xj, ysig, y;
	
	float way1, way2, way3, way4, way5, way6;
	float way7, way8, way9, way10, way11, way12;

	float xay1, xay2, xay3, xay4, xay5, xay6;
	float xay7, xay8, xay9, xay10, xay11, xay12;    

	float yay1, yay2, yay3, yay4, yay5, yay6;
	
	float wy0_1, wy0_2, wy0_3, wy0_4, wy0_5, wy0_6;
	float wy0_7, wy0_8, wy0_9, wy0_10, wy0_11, wy0_12;

	float xy0_1, xy0_2, xy0_3, xy0_4, xy0_5, xy0_6;
	float xy0_7, xy0_8, xy0_9, xy0_10, xy0_11, xy0_12;    

	float yy0_1, yy0_2, yy0_3, yy0_4, yy0_5, yy0_6;
	
	wsig = Win[i];
	xsig = Xin[i];
	ysig = Yin[i];
	
	// pull out some values with a Shift90 - 6 cascading First order sections
	wy0_1 = wsig - (unit->m_coefs[0]) * wy1[0];
	way1 = unit->m_coefs[0] * wy0_1 + 1 * wy1[0]; 
	wy1[0] = wy0_1;
	wy0_2 = way1 - (unit->m_coefs[1]) * wy1[1];
	way2 = unit->m_coefs[1] * wy0_2 + 1 * wy1[1]; 
	wy1[1] = wy0_2;
	wy0_3 = way2 - (unit->m_coefs[2]) * wy1[2];
	way3 = unit->m_coefs[2] * wy0_3 + 1 * wy1[2]; 
	wy1[2] = wy0_3;	
	wy0_4 = way3 - (unit->m_coefs[3]) * wy1[3];
	way4 = unit->m_coefs[3] *  wy0_4 + 1 * wy1[3]; 
	wy1[3] = wy0_4;
	wy0_5 = way4 - (unit->m_coefs[4]) * wy1[4];
	way5 = unit->m_coefs[4] * wy0_5 + 1 * wy1[4]; 
	wy1[4] = wy0_5;	
	wy0_6 = way5 - (unit->m_coefs[5]) * wy1[5];
	w = way6 = unit->m_coefs[5] * wy0_6 + 1 * wy1[5]; 
	wy1[5] = wy0_6;
	
	wy0_7 = wsig - (unit->m_coefs[6]) * wy1[6];
	way7 = unit->m_coefs[6] * wy0_7 + 1 * wy1[6]; 
	wy1[6] = wy0_7;
	wy0_8 = way7 - (unit->m_coefs[7]) * wy1[7];
	way8 = unit->m_coefs[7] * wy0_8 + 1 * wy1[7]; 
	wy1[7] = wy0_8;
	wy0_9 = way8 - (unit->m_coefs[8]) * wy1[8];
	way9 = unit->m_coefs[8] * wy0_9 + 1 * wy1[8]; 
	wy1[8] = wy0_9;
	wy0_10 = way9 - (unit->m_coefs[9]) * wy1[9];
	way10 = unit->m_coefs[9] * wy0_10 + 1 * wy1[9]; 
	wy1[9] = wy0_10;
	wy0_11 = way10 - (unit->m_coefs[10]) * wy1[10];
	way11 = unit->m_coefs[10] * wy0_11  + 1 * wy1[10]; 
	wy1[10] = wy0_11;
	wy0_12 = way11 - (unit->m_coefs[11]) * wy1[11];
	wj = way12 = unit->m_coefs[11] * wy0_12 + 1 * wy1[11]; 
	wy1[11] = wy0_12;
	
	xy0_1 = xsig - (unit->m_coefs[0]) * xy1[0];
	xay1 = unit->m_coefs[0] * xy0_1 + 1 * xy1[0]; 
	xy1[0] = xy0_1;
	xy0_2 = xay1 - (unit->m_coefs[1]) * xy1[1];
	xay2 = unit->m_coefs[1] * xy0_2 + 1 * xy1[1]; 
	xy1[1] = xy0_2;
	xy0_3 = xay2 - (unit->m_coefs[2]) * xy1[2];
	xay3 = unit->m_coefs[2] * xy0_3 + 1 * xy1[2]; 
	xy1[2] = xy0_3;	
	xy0_4 = xay3 - (unit->m_coefs[3]) * xy1[3];
	xay4 = unit->m_coefs[3] * xy0_4 + 1 * xy1[3]; 
	xy1[3] = xy0_4;
	xy0_5 = xay4 - (unit->m_coefs[4]) * xy1[4];
	xay5 = unit->m_coefs[4] * xy0_5 + 1 * xy1[4]; 
	xy1[4] = xy0_5;
	xy0_6 = xay5 - (unit->m_coefs[5]) * xy1[5];
	x = xay6 = unit->m_coefs[5] * xy0_6 + 1 * xy1[5]; 
	xy1[5] = xy0_6;
	
	xy0_7 = xsig - (unit->m_coefs[6]) * xy1[6];
	xay7 = unit->m_coefs[6] * xy0_7 + 1 * xy1[6]; 
	xy1[6] = xy0_7;
	xy0_8 = xay7 - (unit->m_coefs[7]) * xy1[7];
	xay8 = unit->m_coefs[7] * xy0_8 + 1 * xy1[7]; 
	xy1[7] = xy0_8;
	xy0_9 = xay8 - (unit->m_coefs[8]) * xy1[8];
	xay9 = unit->m_coefs[8] * xy0_9 + 1 * xy1[8]; 
	xy1[8] = xy0_9;
	xy0_10 = xay9 - (unit->m_coefs[9]) * xy1[9];
	xay10 = unit->m_coefs[9] * xy0_10  + 1 * xy1[9]; 
	xy1[9] = xy0_10;
	xy0_11 = xay10 - (unit->m_coefs[10]) * xy1[10];
	xay11 = unit->m_coefs[10] * xy0_11 + 1 * xy1[10]; 
	xy1[10] = xy0_11;
	xy0_12 = xay11 - (unit->m_coefs[11]) * xy1[11];
	xj = xay12 = unit->m_coefs[11] * xy0_12 + 1 * xy1[11]; 
	xy1[11] = xy0_12;
	
	yy0_1 = ysig - (unit->m_coefs[0]) * yy1[0];
	yay1 = unit->m_coefs[0] * yy0_1 + 1 * yy1[0]; 
	yy1[0] = yy0_1;
	yy0_2 = yay1 - (unit->m_coefs[1]) * yy1[1];
	yay2 = unit->m_coefs[1] * yy0_2 + 1 * yy1[1]; 
	yy1[1] = yy0_2;
	yy0_3 = yay2 - (unit->m_coefs[2]) * yy1[2];
	yay3 = unit->m_coefs[2] * yy0_3 + 1 * yy1[2]; 
	yy1[2] = yy0_3;
	yy0_4 = yay3 - (unit->m_coefs[3]) * yy1[3];
	yay4 = unit->m_coefs[3] * yy0_4 + 1 * yy1[3]; 
	yy1[3] = yy0_4;
	yy0_5 = yay4 - (unit->m_coefs[4]) * yy1[4];
	yay5 = unit->m_coefs[4] * yy0_5 + 1 * yy1[4]; 
	yy1[4] = yy0_5;
	yy0_6 = yay5 - (unit->m_coefs[5]) * yy1[5];
	y = yay6 = unit->m_coefs[5] * yy0_6 + 1 * yy1[5]; 
	yy1[5] = yy0_6;	
		
	w = 0.9397 * w;
	wj = 0.342 * wj;
	x = 0.1856 * x;
	xj = 0.5099 * xj;
	y = y * 0.655;
		
	lsout[i] = ((((w + x) - wj) + xj) + y) * 0.5;
	rsout[i] = ((((w + x) + wj) - xj) - y) * 0.5;
	}
	
    for(int i = 0; i < 12; ++i){
	unit->m_wy1[i] = zapgremlins(wy1[i]);
	unit->m_xy1[i] = zapgremlins(xy1[i]);
	unit->m_yy1[i] = zapgremlins(yy1[i]);
	}

}


void UHJ2B_Ctor(UHJ2B *unit)
{
    SETCALC(UHJ2B_next);
    float gamconst = (15.0 * pi) / SAMPLERATE;
    float gamma01 = gamconst * 0.3609;
    float gamma02 = gamconst * 2.7412;
    float gamma03 = gamconst * 11.1573;
    float gamma04 = gamconst * 44.7581;
    float gamma05 = gamconst * 179.6242;
    float gamma06 = gamconst * 798.4578;
    float gamma07 = gamconst * 1.2524;
    float gamma08 = gamconst * 5.5671;
    float gamma09 = gamconst * 22.3423;
    float gamma10 = gamconst * 89.6271;
    float gamma11 = gamconst * 364.7914;
    float gamma12 = gamconst * 2770.1114;
    unit->m_coefs[0] = (gamma01 - 1) / (gamma01 + 1);
    unit->m_coefs[1] = (gamma02 - 1) / (gamma02 + 1);
    unit->m_coefs[2] = (gamma03 - 1) / (gamma03 + 1);
    unit->m_coefs[3] = (gamma04 - 1) / (gamma04 + 1);
    unit->m_coefs[4] = (gamma05 - 1) / (gamma05 + 1);
    unit->m_coefs[5] = (gamma06 - 1) / (gamma06 + 1);
    unit->m_coefs[6] = (gamma07 - 1) / (gamma07 + 1);
    unit->m_coefs[7] = (gamma08 - 1) / (gamma08 + 1);
    unit->m_coefs[8] = (gamma09 - 1) / (gamma09 + 1);
    unit->m_coefs[9] = (gamma10 - 1) / (gamma10 + 1);
    unit->m_coefs[10] = (gamma11 - 1) / (gamma11 + 1);
    unit->m_coefs[11] = (gamma12 - 1) / (gamma12 + 1);
    // clear out the storage variables for the filters
    for(int i = 0; i < 12; ++i){
	unit->m_lsy1[i] = 0.0;
	unit->m_rsy1[i] = 0.0;}
    UHJ2B_next(unit, 1);

}

void UHJ2B_next(UHJ2B *unit, int inNumSamples)
{
    float *Lsin = IN(0);
    float *Rsin = IN(1);
    float *Wout = OUT(0);
    float *Xout = OUT(1);
    float *Yout = OUT(2);

    float lsy1[12];
    float rsy1[12];
    float lsig, l, lj, rsig, r, rj;
	
    float lay1, lay2, lay3, lay4, lay5, lay6;
    float lay7, lay8, lay9, lay10, lay11, lay12;

    float ray1, ray2, ray3, ray4, ray5, ray6;
    float ray7, ray8, ray9, ray10, ray11, ray12;    
	
    float ly0_1, ly0_2, ly0_3, ly0_4, ly0_5, ly0_6;
    float ly0_7, ly0_8, ly0_9, ly0_10, ly0_11, ly0_12;

    float ry0_1, ry0_2, ry0_3, ry0_4, ry0_5, ry0_6;
    float ry0_7, ry0_8, ry0_9, ry0_10, ry0_11, ry0_12;    

    // each filter's last sample
    for(int i = 0; i < 12; ++i){
	lsy1[i] = unit->m_lsy1[i];
	rsy1[i] = unit->m_rsy1[i];
	}

    for (int i = 0; i < inNumSamples; ++i){

	lsig = Lsin[i];
	rsig = Rsin[i];
	
	// pull out some values with a Shift90 - 6 cascading First order sections
	ly0_1 = lsig - (unit->m_coefs[0]) * lsy1[0];
	lay1 = unit->m_coefs[0] * ly0_1 + 1 * lsy1[0]; 
	lsy1[0] = ly0_1;
	ly0_2 = lay1 - (unit->m_coefs[1]) * lsy1[1];
	lay2 = unit->m_coefs[1] * ly0_2 + 1 * lsy1[1]; 
	lsy1[1] = ly0_2;
	ly0_3 = lay2 - (unit->m_coefs[2]) * lsy1[2];
	lay3 = unit->m_coefs[2] * ly0_3 + 1 * lsy1[2]; 
	lsy1[2] = ly0_3;	
	ly0_4 = lay3 - (unit->m_coefs[3]) * lsy1[3];
	lay4 = unit->m_coefs[3] *  ly0_4 + 1 * lsy1[3]; 
	lsy1[3] = ly0_4;
	ly0_5 = lay4 - (unit->m_coefs[4]) * lsy1[4];
	lay5 = unit->m_coefs[4] * ly0_5 + 1 * lsy1[4]; 
	lsy1[4] = ly0_5;	
	ly0_6 = lay5 - (unit->m_coefs[5]) * lsy1[5];
	l = lay6 = unit->m_coefs[5] * ly0_6 + 1 * lsy1[5]; 
	lsy1[5] = ly0_6;
	
	ly0_7 = lsig - (unit->m_coefs[6]) * lsy1[6];
	lay7 = unit->m_coefs[6] * ly0_7 + 1 * lsy1[6]; 
	lsy1[6] = ly0_7;
	ly0_8 = lay7 - (unit->m_coefs[7]) * lsy1[7];
	lay8 = unit->m_coefs[7] * ly0_8 + 1 * lsy1[7]; 
	lsy1[7] = ly0_8;
	ly0_9 = lay8 - (unit->m_coefs[8]) * lsy1[8];
	lay9 = unit->m_coefs[8] * ly0_9 + 1 * lsy1[8]; 
	lsy1[8] = ly0_9;
	ly0_10 = lay9 - (unit->m_coefs[9]) * lsy1[9];
	lay10 = unit->m_coefs[9] * ly0_10 + 1 * lsy1[9]; 
	lsy1[9] = ly0_10;
	ly0_11 = lay10 - (unit->m_coefs[10]) * lsy1[10];
	lay11 = unit->m_coefs[10] * ly0_11  + 1 * lsy1[10]; 
	lsy1[10] = ly0_11;
	ly0_12 = lay11 - (unit->m_coefs[11]) * lsy1[11];
	lj = lay12 = unit->m_coefs[11] * ly0_12 + 1 * lsy1[11]; 
	lsy1[11] = ly0_12;

	ry0_1 = rsig - (unit->m_coefs[0]) * rsy1[0];
	ray1 = unit->m_coefs[0] * ry0_1 + 1 * rsy1[0]; 
	rsy1[0] = ry0_1;
	ry0_2 = ray1 - (unit->m_coefs[1]) * rsy1[1];
	ray2 = unit->m_coefs[1] * ry0_2 + 1 * rsy1[1]; 
	rsy1[1] = ry0_2;
	ry0_3 = ray2 - (unit->m_coefs[2]) * rsy1[2];
	ray3 = unit->m_coefs[2] * ry0_3 + 1 * rsy1[2]; 
	rsy1[2] = ry0_3;	
	ry0_4 = ray3 - (unit->m_coefs[3]) * rsy1[3];
	ray4 = unit->m_coefs[3] *  ry0_4 + 1 * rsy1[3]; 
	rsy1[3] = ry0_4;
	ry0_5 = ray4 - (unit->m_coefs[4]) * rsy1[4];
	ray5 = unit->m_coefs[4] * ry0_5 + 1 * rsy1[4]; 
	rsy1[4] = ry0_5;	
	ry0_6 = ray5 - (unit->m_coefs[5]) * rsy1[5];
	r = ray6 = unit->m_coefs[5] * ry0_6 + 1 * rsy1[5]; 
	rsy1[5] = ry0_6;
	
	ry0_7 = rsig - (unit->m_coefs[6]) * rsy1[6];
	ray7 = unit->m_coefs[6] * ry0_7 + 1 * rsy1[6]; 
	rsy1[6] = ry0_7;
	ry0_8 = ray7 - (unit->m_coefs[7]) * rsy1[7];
	ray8 = unit->m_coefs[7] * ry0_8 + 1 * rsy1[7]; 
	rsy1[7] = ry0_8;
	ry0_9 = ray8 - (unit->m_coefs[8]) * rsy1[8];
	ray9 = unit->m_coefs[8] * ry0_9 + 1 * rsy1[8]; 
	rsy1[8] = ry0_9;
	ry0_10 = ray9 - (unit->m_coefs[9]) * rsy1[9];
	ray10 = unit->m_coefs[9] * ry0_10 + 1 * rsy1[9]; 
	rsy1[9] = ry0_10;
	ry0_11 = ray10 - (unit->m_coefs[10]) * rsy1[10];
	ray11 = unit->m_coefs[10] * ry0_11  + 1 * rsy1[10]; 
	rsy1[10] = ry0_11;
	ry0_12 = ray11 - (unit->m_coefs[11]) * rsy1[11];
	rj = ray12 = unit->m_coefs[11] * ry0_12 + 1 * rsy1[11]; 
	rsy1[11] = ry0_12;	
		
	Wout[i] = ((0.982 * l) + (0.982 * r) + (0.164 * lj) - (0.164 * rj)) * 0.5;
	Xout[i] = ((0.419 * l) + (0.419 * r) - (0.828 * lj) + (0.828 * rj)) * 0.5;
	Yout[i] = ((0.763 * l) - (0.763 * r) + (0.385 * lj) + (0.385 * rj)) * 0.5;
	}
	
    for(int i = 0; i < 12; ++i){
	unit->m_lsy1[i] = zapgremlins(lsy1[i]);
	unit->m_rsy1[i] = zapgremlins(rsy1[i]);
	}

}


////////////////////////////////////////////////////////////////////////

void BFManipulate_Ctor(BFManipulate *unit)
{
	SETCALC(BFManipulate_next);
	unit->m_rotate = ZIN0(4);
	unit->m_tilt = ZIN0(5);
	unit->m_tumble = ZIN0(6);
	BFManipulate_next(unit, 1);
}

void BFManipulate_next(BFManipulate *unit, int inNumSamples)
{       
	float *Wout = ZOUT(0);
	float *Xout = ZOUT(1);
	float *Yout = ZOUT(2);
	float *Zout = ZOUT(3);
	
	float *Win = ZIN(0);
	float *Xin = ZIN(1);
	float *Yin = ZIN(2);
	float *Zin = ZIN(3);
	
	float w, x, y, z;
	float sina, sinb, sinc, cosa, cosb, cosc;
	
	float rotate = unit->m_rotate;
	float tilt = unit->m_tilt;
	float tumble = unit->m_tumble;
	
	float endrotate = ZIN0(4);
	float endtilt = ZIN0(5);
	float endtumble = ZIN0(6);
	
	float rotate_slope = CALCSLOPE(endrotate, rotate);
	float tilt_slope = CALCSLOPE(endtilt, tilt);
	float tumble_slope = CALCSLOPE(endtumble, tumble);

	LOOP(inNumSamples,
	    	sina = sin(rotate);
		sinb = sin(tilt);
		sinc = sin(tumble);
	
		cosa = cos(rotate);
		cosb = cos(tilt);
		cosc = cos(tumble);
		
		w = ZXP(Win);
		x = ZXP(Xin);
		y = ZXP(Yin);
		z = ZXP(Zin);
		
		ZXP(Wout) = w;
		ZXP(Xout) = ((x * cosa * cosc) - (y * sina) - (z * sinc));
		ZXP(Yout) = ((x * sina) + (y * cosa * cosb) - (z * sinb));
		ZXP(Zout) = ((x * sinc) + (y * sinb) + (z * cosb * cosc));
		
		rotate += rotate_slope;
		tilt += tilt_slope;
		tumble += tumble_slope;
		);
		
	unit->m_rotate = endrotate;
	unit->m_tilt = endtilt;
	unit->m_tumble = endtumble;
}

///////////////////////////////////////////////////////////////////////////////////////////////////////

void B2Ster_Ctor(B2Ster *unit)
{
    SETCALC(B2Ster_next);
    B2Ster_next(unit, 1);
}

// http://www.cyber.rdg.ac.uk/P.Sharkey/WWW/icdvrat/WWW96/Papers/keating.htm
// L = W + 0.35X + 0.61Y
// R = W + 0.35X - 0.61Y
		
void B2Ster_next(B2Ster *unit, int inNumSamples)
{
    float *Win = IN(0);
    float *Xin = IN(1);
    float *Yin = IN(2);
    float *lout = OUT(0);
    float *rout = OUT(1);

    float w, x, y;
    for(int i = 0; i < inNumSamples; ++i){
	w = Win[i];
	x = Xin[i] * 0.35;
	y = Yin[i] * 0.61;
	
	lout[i] = w + x + y;
	rout[i] = w + x - y;
    }
}

///////////////////////////////////////////////////////////////////////////////////

void BFEncode2_Ctor(BFEncode2 *unit)
{
	SETCALC(BFEncode2_next);
	float sinint, cosint, azimuth, rho;
	float point_x = unit->m_point_x = ZIN0(1);
	float point_y = unit->m_point_y = ZIN0(2);
	float elevation = unit->m_elevation = ZIN0(3);
	float level = unit->m_level = ZIN0(4);

	azimuth = atan2(point_x, point_y);
	rho = hypot(point_x, point_y);
	
	float sina = sin(azimuth);
	float sinb = sin(elevation);
	
	float cosa = cos(azimuth);
	float cosb = cos(elevation);
	
	if(rho >= 1) {
		sinint = (rsqrt2 * (sin(0.78539816339745 * 1.0))) / pow(rho, (float)1.5);
		cosint =  (rsqrt2 * (cos(0.78539816339745 * 1.0))) / pow(rho, (float)1.5);
		}
	    else 
		{
		sinint = rsqrt2 * (sin(0.78539816339745 * rho));
		cosint = rsqrt2 * (cos(0.78539816339745 * rho));
		};

	
	unit->m_W_amp = level * cosint;
	unit->m_X_amp = cosa * cosb * level * sinint;
	unit->m_Y_amp = sina * cosb * level * sinint;
	unit->m_Z_amp = sinb * level * sinint;

	BFEncode2_next(unit, 1);
}

void BFEncode2_next(BFEncode2 *unit, int inNumSamples)
{       
	float sinint, cosint, azimuth, rho;
	float *Wout = ZOUT(0);
	float *Xout = ZOUT(1);
	float *Yout = ZOUT(2);
	float *Zout = ZOUT(3);
	
	float *in = ZIN(0);
	float point_x = ZIN0(1);
	float point_y = ZIN0(2);
	float elevation = ZIN0(3);
	float level = ZIN0(4);

	azimuth = atan2(point_x, point_y);
	rho = hypot(point_x, point_y);

	float W_amp = unit->m_W_amp;
	float X_amp = unit->m_X_amp;
	float Y_amp = unit->m_Y_amp;
	float Z_amp = unit->m_Z_amp;

	if (point_x != unit->m_point_x || point_y != unit->m_point_y ||elevation != unit->m_elevation || level != unit->m_level) {
		unit->m_point_x = point_x;
		unit->m_point_y = point_y;
		unit->m_elevation = elevation;
		unit->m_level = level;

		float sina =  sin(azimuth);
		float sinb =  sin(elevation);

		float cosa = cos(azimuth);
		float cosb = cos(elevation);
		if(rho >= 1) {
			sinint = (rsqrt2 * (sin(0.78539816339745))) / pow(rho, (float)1.5);
			cosint =  (rsqrt2 * (cos(0.78539816339745))) / pow(rho, (float)1.5);
			}
		    else 
			{sinint = rsqrt2 * (sin(0.78539816339745 * rho));
			cosint = rsqrt2 * (cos(0.78539816339745 * rho));
			};
		
		// vary w according to x, y and z
//		float next_W_amp = rsqrt2 * level * cosint;
		float next_W_amp = level * cosint;
		float next_X_amp = cosa * cosb * level * sinint;
		float next_Y_amp = sina * cosb * level * sinint;
		float next_Z_amp = sinb * level * sinint;
		
		float W_slope = CALCSLOPE(next_W_amp, W_amp);
		float X_slope = CALCSLOPE(next_X_amp, X_amp);
		float Y_slope = CALCSLOPE(next_Y_amp, Y_amp);
		float Z_slope = CALCSLOPE(next_Z_amp, Z_amp);
		
		LOOP(inNumSamples, 
			float z = ZXP(in);
			// vary w according to x, y and z
			ZXP(Wout) = z * (W_amp * (1 - (0.293 * ((X_amp * X_amp) + (Y_amp * Y_amp) + (Z_amp * Z_amp)))));
			ZXP(Xout) = z * X_amp;
			ZXP(Yout) = z * Y_amp;
			ZXP(Zout) = z * Z_amp;
			W_amp += W_slope;
			X_amp += X_slope;
			Y_amp += Y_slope;
			Z_amp += Z_slope;
		);
		unit->m_W_amp = W_amp;
		unit->m_X_amp = X_amp;
		unit->m_Y_amp = Y_amp;
		unit->m_Z_amp = Z_amp;
	} else {
		LOOP(inNumSamples, 
			float z = ZXP(in);
			// vary w according to x, y and z
			ZXP(Wout) = z * (W_amp * (1 - (0.293 * ((X_amp * X_amp) + (Y_amp * Y_amp) + (Z_amp * Z_amp)))));
			ZXP(Xout) = z * X_amp;
			ZXP(Yout) = z * Y_amp;
			ZXP(Zout) = z * Z_amp;
		);
	}
}

////////////////////////////////////////////////////////////////////////////////////////////////////////

void BFEncode1_Ctor(BFEncode1 *unit)
{
	SETCALC(BFEncode1_next);
	float sinint, cosint;
	float azimuth = unit->m_azimuth = ZIN0(1);
	float elevation = unit->m_elevation = ZIN0(2);
 	float rho = unit->m_rho = ZIN0(3);
	float level = unit->m_level = ZIN0(4);

	float sina = sin(azimuth);
	float sinb = sin(elevation);
	
	float cosa = cos(azimuth);
	float cosb = cos(elevation);
	
	if(rho >= 1) 
		{sinint = (rsqrt2 * (sin(0.78539816339745))) / pow(rho, (float)1.5);
		cosint =  (rsqrt2 * (cos(0.78539816339745))) / pow(rho, (float)1.5);}
	    else 
		{
		sinint = rsqrt2 * (sin(0.78539816339745 * rho));
		cosint = rsqrt2 * (cos(0.78539816339745 * rho));};
	
	unit->m_W_amp = level * cosint;
	unit->m_X_amp = cosa * cosb * level * sinint;
	unit->m_Y_amp = sina * cosb * level * sinint;
	unit->m_Z_amp = sinb * level * sinint;

	BFEncode1_next(unit, 1);
}

void BFEncode1_next(BFEncode1 *unit, int inNumSamples)
{       
	float sinint, cosint;
	float *Wout = ZOUT(0);
	float *Xout = ZOUT(1);
	float *Yout = ZOUT(2);
	float *Zout = ZOUT(3);
	
	float *in = ZIN(0);
	float azimuth = ZIN0(1);
	float elevation = ZIN0(2);
	float rho = ZIN0(3);
	float level = ZIN0(4);

	float W_amp = unit->m_W_amp;
	float X_amp = unit->m_X_amp;
	float Y_amp = unit->m_Y_amp;
	float Z_amp = unit->m_Z_amp;

	if (azimuth != unit->m_azimuth || rho != unit->m_rho || elevation != unit->m_elevation || level != unit->m_level) {
		unit->m_azimuth = azimuth;
		unit->m_elevation = elevation;
		unit->m_level = level;
		unit->m_rho = rho;
		
		float sina = sin(azimuth);
		float sinb = sin(elevation);
		
		float cosa = cos(azimuth);
		float cosb = cos(elevation);
		
		if(rho >= 1) {
			sinint = (rsqrt2 * (sin(0.78539816339745))) / pow(rho, (float)1.5);
			cosint =  (rsqrt2 * (cos(0.78539816339745))) / pow(rho, (float)1.5);}
		    else { 
			sinint = rsqrt2 * (sin(0.78539816339745 * rho));
			cosint = rsqrt2 * (cos(0.78539816339745 * rho));
			};
		
		float next_W_amp = level * cosint;
		float next_X_amp = cosa * cosb * level * sinint;
		float next_Y_amp = sina * cosb * level * sinint;
		float next_Z_amp = sinb * level * sinint;
		
		float W_slope = CALCSLOPE(next_W_amp, W_amp);
		float X_slope = CALCSLOPE(next_X_amp, X_amp);
		float Y_slope = CALCSLOPE(next_Y_amp, Y_amp);
		float Z_slope = CALCSLOPE(next_Z_amp, Z_amp);
		
		LOOP(inNumSamples, 
			float z = ZXP(in);
			// vary w according to x, y and z
			ZXP(Wout) = z * (W_amp * (1 - (0.293 * ((X_amp * X_amp) + (Y_amp * Y_amp) + (Z_amp * Z_amp)))));			
			ZXP(Xout) = z * X_amp;
			ZXP(Yout) = z * Y_amp;
			ZXP(Zout) = z * Z_amp;
			W_amp += W_slope;
			X_amp += X_slope;
			Y_amp += Y_slope;
			Z_amp += Z_slope;
		);
		unit->m_W_amp = W_amp;
		unit->m_X_amp = X_amp;
		unit->m_Y_amp = Y_amp;
		unit->m_Z_amp = Z_amp;
	} else {
		LOOP(inNumSamples, 
			float z = ZXP(in);
			// vary w according to x, y and z
			ZXP(Wout) = z * (W_amp * (1 - (0.293 * ((X_amp * X_amp) + (Y_amp * Y_amp) + (Z_amp * Z_amp)))));
			ZXP(Xout) = z * X_amp;
			ZXP(Yout) = z * Y_amp;
			ZXP(Zout) = z * Z_amp;
		);
	}
}


////////////////////////////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////////////////////////////
/*
void BFEncodeSter_Ctor(BFEncodeSter *unit)
{
	SETCALC(BFEncodeSter_next);
	float sinint, cosint;
	float spread = unit->m_spread = IN0(2);
	float azimuth = unit->m_azimuth = IN0(3);
	float elevation = unit->m_elevation = IN0(4);
 	float rho = unit->m_rho = IN0(5);
	float level = unit->m_level = IN0(6);

	float sina = sin(azimuth);
	float sinb = sin(elevation);
	
	float cosa = cos(azimuth);
	float cosb = cos(elevation);
	
	if(rho >= 1) 
		{sinint = (rsqrt2 * (sin(0.78539816339745))) / pow(rho, 1.5);
		cosint =  (rsqrt2 * (cos(0.78539816339745))) / pow(rho, 1.5);}
	    else 
		{
		sinint = rsqrt2 * (sin(0.78539816339745 * rho));
		cosint = rsqrt2 * (cos(0.78539816339745 * rho));};
	
	unit->m_W_amp = level * cosint;
	unit->m_X_amp = cosa * cosb * level * sinint;
	unit->m_Y_amp = sina * cosb * level * sinint;
	unit->m_Z_amp = sinb * level * sinint;

	BFEncodeSter_next(unit, 1);
}

void BFEncodeSter_next(BFEncodeSter *unit, int inNumSamples)
{       
	float sinint, cosint, curr, curl;
	float *Wout = OUT(0);
	float *Xout = OUT(1);
	float *Yout = OUT(2);
	float *Zout = OUT(3);
	
	float *l = IN(0); // the left part of the signal
	float *r = IN(1); // the right part
	float spread = IN(2); // the spread (in radians) between the signals
	float azimuth = IN0(3); // the center point of the signals
	float elevation = IN0(4);
	float rho = IN0(5);
	float level = IN0(6);

	float W_ampL = unit->m_W_ampL;
	float X_ampL = unit->m_X_ampL;
	float Y_ampL = unit->m_Y_ampL;
	float Z_ampL = unit->m_Z_ampL;

	if (spread != unit->m_spread, azimuth != unit->m_azimuth || rho != unit->m_rho || elevation != unit->m_elevation || 
			level != unit->m_level) {
		unit->m_spread = spread
		unit->m_azimuth = azimuth;
		unit->m_elevation = elevation;
		unit->m_level = level;
		unit->m_rho = rho;
		
		float sinaL = sin(azimuth - (spread * 0.5));
		float sinbL = sin(elevation);
		
		float cosaL = cos(azimuth + (spread * 0.5));
		float cosbL = cos(elevation);
		
		if(rho >= 1) {
			sinint = (rsqrt2 * (sin(0.78539816339745))) / pow(rho, 1.5);
			cosint =  (rsqrt2 * (cos(0.78539816339745))) / pow(rho, 1.5);}
		    else { 
			sinint = rsqrt2 * (sin(0.78539816339745 * rho));
			cosint = rsqrt2 * (cos(0.78539816339745 * rho));
			};
		
		float next_W_amp = level * cosint;
		float next_X_amp = cosa * cosb * level * sinint;
		float next_Y_amp = sina * cosb * level * sinint;
		float next_Z_amp = sinb * level * sinint;
		
		float W_slope = CALCSLOPE(next_W_amp, W_amp);
		float X_slope = CALCSLOPE(next_X_amp, X_amp);
		float Y_slope = CALCSLOPE(next_Y_amp, Y_amp);
		float Z_slope = CALCSLOPE(next_Z_amp, Z_amp);
		
		for (int i = 0; i < inNumSamples; ++i){
			curl = l[i];
			curr = r[i]; /// 'z' was the mono input sig
			// vary w according to x, y and z
			Wout[i] = curl * (W_amp * (1 - (0.293 * ((X_amp * X_amp) + (Y_amp * Y_amp) + (Z_amp * Z_amp)))));			
			Xout[i] = curl * X_amp;
			Yout[i] = curl * Y_amp;
			Zout[i] = curl * Z_amp;
			W_amp += W_slope;
			X_amp += X_slope;
			Y_amp += Y_slope;
			Z_amp += Z_slope;
		};
		unit->m_W_amp = W_amp;
		unit->m_X_amp = X_amp;
		unit->m_Y_amp = Y_amp;
		unit->m_Z_amp = Z_amp;
	} else {
		for (int i = 0; i < inNumSamples; ++i){
			float z = ZXP(in);
			// vary w according to x, y and z
			ZXP(Wout) = z * (W_amp * (1 - (0.293 * ((X_amp * X_amp) + (Y_amp * Y_amp) + (Z_amp * Z_amp)))));
			ZXP(Xout) = z * X_amp;
			ZXP(Yout) = z * Y_amp;
			ZXP(Zout) = z * Z_amp;
		};
	}
}
*/
////////////////////////////////////////////////////////////////////////////////////////////////////////

void BFDecode1_Ctor(BFDecode1 *unit)
{
	SETCALC(BFDecode1_next);

	BFDecode1_next(unit, 1);
	
	float angle = ZIN0(4);
	float elevation = ZIN0(5);
	
	unit->m_cosa = cos(angle);
	unit->m_sina = sin(angle);
	unit->m_sinb = sin(elevation);
}

void BFDecode1_next(BFDecode1 *unit, int inNumSamples)
{
	float *Win0 = ZIN(0);
	float *Xin0 = ZIN(1);
	float *Yin0 = ZIN(2);
	float *Zin0 = ZIN(3);
	float *out = ZOUT(0);

	float cosa = unit->m_cosa;
	float sina = unit->m_sina;
	float sinb = unit->m_sinb;
	
	float W_amp = 1.0; //0.7071067811865476;
	float X_amp = cosa; //0.5 * cosa;
	float Y_amp = sina; //0.5 * sina;
	float Z_amp = sinb; //0.5 * sinb;

	LOOP(inNumSamples, 
		ZXP(out) = (ZXP(Win0) * W_amp) +    
			(ZXP(Xin0) * X_amp) + 
			(ZXP(Yin0) * Y_amp) + 
			(ZXP(Zin0) * Z_amp);
	)
}


////////////////////////////////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////////////////////////

static float cubicinterp(float x, float y0, float y1, float y2, float y3)
{
	// 4-point, 3rd-order Hermite (x-form)
	float c0 = y1;
	float c1 = 0.5f * (y2 - y0);
	float c2 = y0 - 2.5f * y1 + 2.f * y2 - 0.5f * y3;
	float c3 = 0.5f * (y3 - y0) + 1.5f * (y1 - y2);

	return ((c3 * x + c2) * x + c1) * x + c0;
}

inline float IN_AT(Unit* unit, int index, int offset) 
{
	if (INRATE(index) == calc_FullRate) return IN(index)[offset];
	return ZIN0(index);
}

inline double sc_gloop(double in, double hi) 
{
	// avoid the divide if possible
	if (in >= hi) {
		in -= hi;
		if (in < hi) return in;
	} else if (in < 0.) {
		in += hi;
		if (in >= 0.) return in;
	} else return in;
	
	return in - hi * floor(in/hi); 
}

#define GRAIN_BUF \
	SndBuf *buf = bufs + bufnum; \
	float *bufData __attribute__((__unused__)) = buf->data; \
	uint32 bufChannels __attribute__((__unused__)) = buf->channels; \
	uint32 bufSamples __attribute__((__unused__)) = buf->samples; \
	uint32 bufFrames = buf->frames; \
	int guardFrame __attribute__((__unused__)) = bufFrames - 2; \

#define GRAIN_LOOP_BODY_4 \
		float amp = y1 * y1; \
		phase = sc_gloop(phase, loopMax); \
		int32 iphase = (int32)phase; \
		float* table1 = bufData + iphase; \
		float* table0 = table1 - 1; \
		float* table2 = table1 + 1; \
		float* table3 = table1 + 2; \
		if (iphase == 0) { \
			table0 += bufSamples; \
		} else if (iphase >= guardFrame) { \
			if (iphase == guardFrame) { \
				table3 -= bufSamples; \
			} else { \
				table2 -= bufSamples; \
				table3 -= bufSamples; \
			} \
		} \
		float fracphase = phase - (double)iphase; \
		float a = table0[0]; \
		float b = table1[0]; \
		float c = table2[0]; \
		float d = table3[0]; \
		float outval = amp * cubicinterp(fracphase, a, b, c, d); \
		ZXP(out1) += outval; \
		double y0 = b1 * y1 - y2; \
		y2 = y1; \
		y1 = y0; \


#define GRAIN_LOOP_BODY_2 \
		float amp = y1 * y1; \
		phase = sc_gloop(phase, loopMax); \
		int32 iphase = (int32)phase; \
		float* table1 = bufData + iphase; \
		float* table2 = table1 + 1; \
		if (iphase > guardFrame) { \
			table2 -= bufSamples; \
		} \
		float fracphase = phase - (double)iphase; \
		float b = table1[0]; \
		float c = table2[0]; \
		float outval = amp * (b + fracphase * (c - b)); \
		ZXP(out1) += outval; \
		double y0 = b1 * y1 - y2; \
		y2 = y1; \
		y1 = y0; \
		

#define GRAIN_LOOP_BODY_1 \
		float amp = y1 * y1; \
		phase = sc_gloop(phase, loopMax); \
		int32 iphase = (int32)phase; \
		float outval = amp * bufData[iphase]; \
		ZXP(out1) += outval; \
		double y0 = b1 * y1 - y2; \
		y2 = y1; \
		y1 = y0; \
	
	
void Warp1_next(Warp1 *unit, int inNumSamples)
{
	ClearUnitOutputs(unit, inNumSamples);
	float *out;
	out = ZOUT(0);
	int nextGrain = unit->mNextGrain;
	
	World *world = unit->mWorld;
	SndBuf *bufs = world->mSndBufs;
	uint32 numBufs = world->mNumSndBufs;
	
	for (int i=0; i < unit->mNumActive; ) {
		Grain *grain = unit->mGrains + i;
		uint32 bufnum = grain->bufnum;
		
		GRAIN_BUF
		
		if (bufChannels != 1) {
			 ++i;
			 continue;
		}

		double loopMax = (double)bufFrames;

		double rate = grain->rate;
		double phase = grain->phase;
		double b1 = grain->b1;
		double y1 = grain->y1;
		double y2 = grain->y2;
		
		float *out1 = out;

		int nsmps = sc_min(grain->counter, inNumSamples);
		if (grain->interp >= 4) {
			for (int j=0; j<nsmps; ++j) {
				GRAIN_LOOP_BODY_4;
				phase += rate;
			}
		} else if (grain->interp >= 2) {
			for (int j=0; j<nsmps; ++j) {
				GRAIN_LOOP_BODY_2;
				phase += rate;
			}
		} else {
			for (int j=0; j<nsmps; ++j) {
				GRAIN_LOOP_BODY_1;
				phase += rate;
			}
		}
		
		grain->phase = phase;
		grain->y1 = y1;
		grain->y2 = y2;

		grain->counter -= nsmps;
		if (grain->counter <= 0) {
			// remove grain
			*grain = unit->mGrains[--unit->mNumActive];
		} else ++i;
	}
	
	for (int i=0; i<inNumSamples; ++i) {
		--nextGrain;
		if (nextGrain == 0) {
			// start a grain
			if (unit->mNumActive+1 >= kMaxGrains) break;
			uint32 bufnum = (uint32)IN_AT(unit, 0, i);
			if (bufnum >= numBufs) continue;
			
			GRAIN_BUF
			
			if (bufChannels != 1) continue;

			float bufSampleRate = buf->samplerate;
			float bufRateScale = bufSampleRate * SAMPLEDUR;
			double loopMax = (double)bufFrames;

			Grain *grain = unit->mGrains + unit->mNumActive++;
			grain->bufnum = bufnum;
			
			RGET
			
			s1 = trand(s1, s2, s3);
			s2 = trand(s2, s3, s1);
			s3 = trand(s3, s1, s2);
			

			double overlaps = IN_AT(unit, 4, i);
			double counter = IN_AT(unit, 3, i) * SAMPLERATE;
			double winrandamt = frand2(s1, s2, s3) * IN_AT(unit, 5, i);
			counter = floor(counter + (counter * winrandamt));
			counter = sc_max(4., counter);
			grain->counter = (int)counter;
			
			nextGrain = (int)(counter / overlaps);
			
			unit->mNextGrain = nextGrain;
			
			double rate = grain->rate = IN_AT(unit, 2, i) * bufRateScale;
			double phase = IN_AT(unit, 1, i) * bufFrames;
			grain->interp = (int)IN_AT(unit, 6, i);
			
			double w = pi / counter;
			double b1 = grain->b1 = 2. * cos(w);
			double y1 = sin(w);
			double y2 = 0.;
			
			float *out1 = out + i;
			
			int nsmps = sc_min(grain->counter, inNumSamples - i);
			if (grain->interp >= 4) {
				for (int j=0; j<nsmps; ++j) {
					GRAIN_LOOP_BODY_4;
					phase += rate;
				}
			} else if (grain->interp >= 2) {
				for (int j=0; j<nsmps; ++j) {
					GRAIN_LOOP_BODY_2;
					phase += rate;
				}
			} else {
				for (int j=0; j<nsmps; ++j) {
					GRAIN_LOOP_BODY_1;
					phase += rate;
				}
			}
			
			grain->phase = phase;
			grain->y1 = y1;
			grain->y2 = y2;
			// store random values
			RPUT
			// end change
			grain->counter -= nsmps;
			if (grain->counter <= 0) {
				// remove grain
				*grain = unit->mGrains[--unit->mNumActive];
			}
		}
	}	
	
	unit->mNextGrain = nextGrain;
}

void Warp1_Ctor(Warp1 *unit)
{	
	SETCALC(Warp1_next);
	
	unit->mNumActive = 0;
	unit->mNextGrain = 1;
	
	ClearUnitOutputs(unit, 1);
}

////////////////////////////////////////////////////////////////////////////////////////////////////////
// LPC Ugens built by Richard Karpen, Don Craig and Josh Parmenter
////////////////////////////////////////////////////////////////////////////////////////////////////////

void LPCVals_Ctor(LPCVals *unit)
{
	if (INRATE(1) == calc_FullRate) {
	    SETCALC(LPCVals_next_a);
	    } else {
	    SETCALC(LPCVals_next_k);
	    }
	ClearUnitOutputs(unit, 1);
	unit->m_fbufnum = -1e9f;
}

void LPCVals_next_k(LPCVals *unit, int inNumSamples)
{
	float fbufnum  = ZIN0(0);
	
	if (fbufnum != unit->m_fbufnum) {
		uint32 bufnum = (int)fbufnum;
		World *world = unit->mWorld;
		if (bufnum >= world->mNumSndBufs) bufnum = 0;
		unit->m_fbufnum = fbufnum;
		unit->m_buf = world->mSndBufs + bufnum;
		}
		
	SndBuf *buf = unit->m_buf;
	float *bufData __attribute__((__unused__)) = buf->data; 
	float* fileFrames = bufData + 1;
	int numframes = (int)fileFrames[0];
	int rmsoStart = 3 + numframes;
	int nrmerrStart = 3 + (numframes * 2);
	int pchcpsStart = 3 + (numframes * 3);
	
	float *cpsout = OUT(0);
	float *rmsoout = OUT(1);
	float *errout = OUT(2);
	
	int intframes = ((int)numframes - 1);
	int frameround1, frameround2;
	float frame = unit->m_framestart * intframes;
	float frameend = IN0(1);

	float cps1, cps2, rmso1, rmso2, err1, err2, framepct;
	
	float frameinc = CALCSLOPE(frameend * intframes, frame);
		    	
	if (!bufData) { 
                unit->mDone = true; 
		ClearUnitOutputs(unit, inNumSamples); 
		return;
		};
	
	for (int i = 0; i < inNumSamples; ++i){

	    frameround1 = (int)frame; 
	    frameround2 = frameround1 + 1;
	    framepct = frame-(float)frameround1;
	    if (frameround2 > numframes) frameround2 = frameround1;

	    cps1 = *(bufData + (pchcpsStart + frameround1));
	    cps2 = *(bufData + (pchcpsStart + frameround2));
	    cpsout[i] = cps1 + ((cps2 - cps1) * framepct);
	
	    rmso1 = *(bufData + (rmsoStart + frameround1));
	    rmso2 = *(bufData + (rmsoStart + frameround2));
	    rmsoout[i] = rmso1 + ((rmso2 - rmso1) * framepct);
	
	    err1 = *(bufData + (nrmerrStart + frameround1));
	    err2 = *(bufData + (nrmerrStart + frameround2));
	    errout[i] = err1 + ((err2 - err1) * framepct);
	    
	    frame += frameinc;
	    }
	    
	    unit->m_framestart = frameend;
}


void LPCVals_next_a(LPCVals *unit, int inNumSamples)
{
	float fbufnum  = ZIN0(0);
	
	if (fbufnum != unit->m_fbufnum) {
		uint32 bufnum = (int)fbufnum;
		World *world = unit->mWorld;
		if (bufnum >= world->mNumSndBufs) bufnum = 0;
		unit->m_fbufnum = fbufnum;
		unit->m_buf = world->mSndBufs + bufnum;
		}
		
	SndBuf *buf = unit->m_buf;
	float *bufData __attribute__((__unused__)) = buf->data; 
	float* fileFrames = bufData + 1;
	int numframes = (int)fileFrames[0];
	int rmsoStart = 3 + numframes;
	int nrmerrStart = 3 + (numframes * 2);
	int pchcpsStart = 3 + (numframes * 3);
	
	float *cpsout = OUT(0);
	float *rmsoout = OUT(1);
	float *errout = OUT(2);
	
	int intframes = ((int)numframes - 1);
	int frameround1, frameround2;
	float frame; //= unit->m_framestart * intframes;
	float *framein = IN(1);

	float cps1, cps2, rmso1, rmso2, err1, err2, framepct;
	
		    	
	if (!bufData) { 
                unit->mDone = true; 
		ClearUnitOutputs(unit, inNumSamples); 
		return;
		};
	
	for (int i = 0; i < inNumSamples; ++i){
	    frame = framein[i] * (float)intframes;
	    frameround1 = (int)frame; 
	    frameround2 = frameround1 + 1;
	    framepct = frame-(float)frameround1;
	    if (frameround2 > numframes) frameround2 = frameround1;

	    cps1 = *(bufData + (pchcpsStart + frameround1));
	    cps2 = *(bufData + (pchcpsStart + frameround2));
	    cpsout[i] = cps1 + ((cps2 - cps1) * framepct);
	
	    rmso1 = *(bufData + (rmsoStart + frameround1));
	    rmso2 = *(bufData + (rmsoStart + frameround2));
	    rmsoout[i] = rmso1 + ((rmso2 - rmso1) * framepct);
	
	    err1 = *(bufData + (nrmerrStart + frameround1));
	    err2 = *(bufData + (nrmerrStart + frameround2));
	    errout[i] = err1 + ((err2 - err1) * framepct);
	    
	    }
	    
}

void LPCSynth_Ctor(LPCSynth *unit)
{
	if (INRATE(2) == calc_FullRate) {
	    SETCALC(LPCSynth_next_a);
	    } else {
	    SETCALC(LPCSynth_next_k);
	    }
	    
	ClearUnitOutputs(unit, 1);
	unit->m_fbufnum = -1e9f;
	unit->m_valindex = 0;
	unit->m_framestart = IN0(2);
	// init the past sample's values to 0.0
	for (int i = 0; i < 100; ++i){
	    unit->m_storevals[i] = 0.0f;
	    }

}

void LPCSynth_next_k(LPCSynth *unit, int inNumSamples)
{
	float fbufnum  = ZIN0(0);
	
	if (fbufnum != unit->m_fbufnum) {
		uint32 bufnum = (int)fbufnum;
		World *world = unit->mWorld;
		if (bufnum >= world->mNumSndBufs) bufnum = 0;
		unit->m_fbufnum = fbufnum;
		unit->m_buf = world->mSndBufs + bufnum;
	}
	
	SndBuf *buf = unit->m_buf;
	float *bufData __attribute__((__unused__)) = buf->data; 
	float start;
	float *out = OUT(0);
	float *signal = IN(1);
	float* numPoles = bufData;
	float* fileFrames = bufData + 1;
	int numframes = (int)fileFrames[0];
	int coefsStart = 3 + (numframes * 4);
	int npoles = (int)numPoles[0];

	float *polevals = new float[npoles];
	int valindex = unit->m_valindex;
	int intframes = ((int)numframes - 1);
	float frame = unit->m_framestart * intframes;
	float frameend = IN0(2);
	float frameendscale = frameend * intframes;
	float coef1, coef2;
	float framepct, frameinc;
	int frameround1, frameround2;
	
    	frameinc = CALCSLOPE(frameendscale, frame);
	
	if (!bufData) { 
                unit->mDone = true; 
		ClearUnitOutputs(unit, inNumSamples); 
		return;
		}; 

	//the all-pole filter	
	for (int i = 0; i < inNumSamples; ++i){
		// get values for interpolating the current frames coefs
		frameround1 = (int)frame; 
		frameround2 = frameround1 + 1;
		framepct = frame-(float)frameround1;
		if (frameround2 > numframes) frameround2 = frameround1;

	    //fill polevals with interpolated coefs
	    for (int j = 0; j < (int)npoles; j++){
		coef1 = *(bufData + (coefsStart + ((int)numframes * j) + frameround1));
		coef2 = *(bufData + (coefsStart + ((int)numframes * j) + frameround2));
		polevals[j] = (coef1 + ((coef2 - coef1) * framepct)); // lin interp for a pole val between frames
	    }   
	    start = signal[i];
    	    frame += frameinc;
	    for (int j = valindex, n=0; n < npoles; n++, j++)
		start +=  polevals[n] * unit->m_storevals[j];
		out[i] = unit->m_storevals[valindex] = unit->m_storevals[valindex + npoles] = zapgremlins(start);
		valindex = (valindex + 1) % npoles;
	    }
	
	unit->m_framestart = frameend;
	unit->m_valindex = valindex;
	delete [] polevals;
}

void LPCSynth_next_a(LPCSynth *unit, int inNumSamples)
{
	float fbufnum  = ZIN0(0);
	
	if (fbufnum != unit->m_fbufnum) {
		uint32 bufnum = (int)fbufnum;
		World *world = unit->mWorld;
		if (bufnum >= world->mNumSndBufs) bufnum = 0;
		unit->m_fbufnum = fbufnum;
		unit->m_buf = world->mSndBufs + bufnum;
	}
	
	SndBuf *buf = unit->m_buf;
	float *bufData __attribute__((__unused__)) = buf->data; 
	float start;
	float *out = OUT(0);
	float *signal = IN(1);
	float* numPoles = bufData;
	float* fileFrames = bufData + 1;
	int numframes = (int)fileFrames[0];
	int coefsStart = 3 + (numframes * 4);
	int npoles = (int)numPoles[0];

	float *polevals = new float[npoles];
	int valindex = unit->m_valindex;
	int intframes = ((int)numframes - 1);
	float *framein = IN(2);
	float coef1, coef2;
	float framepct, frame;
	int frameround1, frameround2;
		
	if (!bufData) { 
                unit->mDone = true; 
		ClearUnitOutputs(unit, inNumSamples); 
		delete [] polevals;
		return;
		}; 

//the all-pole filter	
	for (int i = 0; i < inNumSamples; ++i){
		// get values for interpolating the current frames coefs
		frame = framein[i] * intframes;
		frameround1 = (int)frame; 
		frameround2 = frameround1 + 1;
		framepct = frame-(float)frameround1;
		if (frameround2 > numframes) frameround2 = frameround1;

	    //fill polevals with interpolated coefs
	    for (int j = 0; j < (int)npoles; j++){
		coef1 = *(bufData + (coefsStart + ((int)numframes * j) + frameround1));
		coef2 = *(bufData + (coefsStart + ((int)numframes * j) + frameround2));
		polevals[j] = (coef1 + ((coef2 - coef1) * framepct)); // lin interp for a pole val between frames
	    }   
	    start = signal[i];

	    for (int j = valindex, n=0; n < npoles; n++, j++)
		start +=  polevals[n] * unit->m_storevals[j];
		out[i] = unit->m_storevals[valindex] = unit->m_storevals[valindex + npoles] = zapgremlins(start);
		valindex = (valindex + 1) % npoles;
	    }
	
	unit->m_valindex = valindex;
	delete [] polevals;
}

/////////////////////////////////////////////////////////////////////////////////////////////////

void MonoGrain_next(MonoGrain *unit, int inNumSamples)
{
	ClearUnitOutputs(unit, inNumSamples);
	float *out = OUT(0);
	float *in = IN(0);
	float winSize = IN0(1);
	float grainFreq = IN0(2);
	float randSize = IN0(3);
	for (int i=0; i < unit->mNumActive; ) {
		Grain *grain = unit->mGrains + i;
		double b1 = grain->b1;
		double y1 = grain->y1;
		double y2 = grain->y2;
		float *out1 = out; // + i;
		int nsmps = sc_min(grain->counter, inNumSamples);
		for (int j=0; j<nsmps; ++j) {
		    float amp = y1 * y1;
		    float outval = amp * in[j];
		    out1[j] += outval;
		    double y0 = b1 * y1 - y2;
		    y2 = y1;
		    y1 = y0;
		    }
		grain->y1 = y1;
		grain->y2 = y2;
		grain->counter -= nsmps;
		if (grain->counter <= 0) {
			// remove grain
			*grain = unit->mGrains[--unit->mNumActive];
		} else ++i;
	}
	int nextGrain = unit->m_NextGrain;
	for (int i=0; i<inNumSamples; ++i) {
		--nextGrain;
		if (nextGrain <= 0) {
			// start a grain
			if (unit->mNumActive+1 >= kMaxGrains) break;
    			Grain *grain = unit->mGrains + unit->mNumActive++;
			RGET
			s1 = trand(s1, s2, s3);
			s2 = trand(s2, s3, s1);
			s3 = trand(s3, s1, s2);
			if (grainFreq <= 0.0) {printf("GrainFreq must be greater than 0"); break;}
			double counter = winSize * SAMPLERATE;
			double winrandamt = frand2(s1, s2, s3) * randSize;
			counter = floor(counter + (counter * winrandamt));
			counter = sc_max(4., counter);
			grain->counter = (int)counter;
			nextGrain = (int)(SAMPLERATE / grainFreq);
			double w = pi / counter;
			double b1 = grain->b1 = 2. * cos(w);
			double y1 = sin(w);
			double y2 = 0.;
			float *out1 = out; // + i;
			int nsmps = sc_min(grain->counter, inNumSamples - i);
			for (int j=0; j<nsmps; ++j) {
			    float amp = y1 * y1;
			    float outval = amp * in[j];
			    out1[j] += outval;
			    double y0 = b1 * y1 - y2;
			    y2 = y1;
			    y1 = y0;
			}	
			grain->y1 = y1;
			grain->y2 = y2;
			// store random values
			RPUT
			// end change
			grain->counter -= nsmps;
			if (grain->counter <= 0) {
				// remove grain
				*grain = unit->mGrains[--unit->mNumActive];
			}
		}
	}	
	unit->m_NextGrain = nextGrain;
}

void MonoGrain_Ctor(MonoGrain *unit)
{	
	SETCALC(MonoGrain_next);
	unit->mNumActive = 0;
	unit->m_NextGrain = 1;
	ClearUnitOutputs(unit, 1);
}

// a BFormat version of the previous MonoGrain
void MonoGrainBF_next(MonoGrainBF *unit, int inNumSamples)
{
	ClearUnitOutputs(unit, inNumSamples);
	float *Wout = OUT(0);
	float *Xout = OUT(1);
	float *Yout = OUT(2);
	float *Zout = OUT(3);
	float *in = IN(0);
	float winSize = IN0(1);
	float grainFreq = IN0(2);
	float randSize = IN0(3);
	float loc = IN0(4); // same as azimuth
	float locrand = IN0(5);
	float elevation = IN0(6); 
	float elevationrand = IN0(7);
	float rho = IN0(8); // place at the speaker radius for now
	
	for (int i=0; i < unit->mNumActive; ) {
		Grainbf *grainbf = unit->mGrains + i;
		float sinint, cosint;
		double b1 = grainbf->b1;
		double y1 = grainbf->y1;
		double y2 = grainbf->y2;
		double thisloc = grainbf->loc;
		double thiselevation = grainbf->elevation;
		float W_amp = grainbf->m_W_amp;
		float X_amp = grainbf->m_X_amp;
		float Y_amp = grainbf->m_Y_amp;
		float Z_amp = grainbf->m_Z_amp;
		float *Wout1 = Wout; 
		float *Xout1 = Xout; 
		float *Yout1 = Yout; 
		float *Zout1 = Zout; 
    		float sina = sin(thisloc);
		float sinb = sin(thiselevation);		
		float cosa = cos(thisloc);
		float cosb = cos(thiselevation);
		if(rho >= 1) 
			{sinint = (rsqrt2 * (sin(0.78539816339745 * 1.0))) / pow(rho, (float)1.5);
			cosint =  (rsqrt2 * (cos(0.78539816339745 * 1.0))) / pow(rho, (float)1.5);}
		    else 
			{sinint = rsqrt2 * (sin(0.78539816339745 * rho));
			cosint = rsqrt2 * (cos(0.78539816339745 * rho));
			};
		float next_W_amp = rsqrt2 * cosint;
		float next_X_amp = cosa * cosb * sinint;
		float next_Y_amp = sina * cosb * sinint;
		float next_Z_amp = sinb * sinint;
		float W_slope = CALCSLOPE(next_W_amp, W_amp);
		float X_slope = CALCSLOPE(next_X_amp, X_amp);
		float Y_slope = CALCSLOPE(next_Y_amp, Y_amp);
		float Z_slope = CALCSLOPE(next_Z_amp, Z_amp);
		int nsmps = sc_min(grainbf->counter, inNumSamples);
		for (int j=0; j<nsmps; ++j) {
		    float amp = y1 * y1;
		    float outval = amp * in[j];
		    // do the BFormat encoding on the resultant signal
		    Wout1[j] += outval * W_amp;
		    Xout1[j] += outval * X_amp;
		    Yout1[j] += outval * Y_amp;
		    Zout1[j] += outval * Z_amp;
		    W_amp += W_slope;
		    X_amp += X_slope;
		    Y_amp += Y_slope;
		    Z_amp += Z_slope;
		    double y0 = b1 * y1 - y2;
		    y2 = y1;
		    y1 = y0;
		    }
		grainbf->y1 = y1;
		grainbf->y2 = y2;
		grainbf->counter -= nsmps;
		grainbf->m_W_amp = W_amp;
		grainbf->m_X_amp = X_amp;
		grainbf->m_Y_amp = Y_amp;
		grainbf->m_Z_amp = Z_amp;
		if (grainbf->counter <= 0) {
			// remove grain
			*grainbf = unit->mGrains[--unit->mNumActive];
		} else ++i;
	}
	int nextGrain = unit->m_NextGrain;
	for (int i=0; i<inNumSamples; ++i) {
		--nextGrain;
		if (nextGrain <= 0) {
			// start a grain
			if (unit->mNumActive+1 >= kMaxGrains) break;
			float sinint, cosint;
    			Grainbf *grainbf = unit->mGrains + unit->mNumActive++;
			RGET
			s1 = trand(s1, s2, s3);
			s2 = trand(s2, s3, s1);
			s3 = trand(s3, s1, s2);
			if (grainFreq <= 0.0) {printf("GrainFreq must be greater than 0"); break;}
			double counter = winSize * SAMPLERATE;
			double ranval = frand2(s1, s2, s3);
			double winrandamt = ranval * randSize;
			counter = floor(counter + (counter * winrandamt));
			counter = sc_max(4., counter);
			grainbf->counter = (int)counter;
			nextGrain = (int)(SAMPLERATE / grainFreq);
			double w = pi / counter;
			double b1 = grainbf->b1 = 2. * cos(w);
			double y1 = sin(w);
			double y2 = 0.;
			double thisloc = loc + (locrand * (ranval * pi));
			double thiselevation = elevation + (elevationrand * (ranval * 0.5 * pi));
			float *Wout1 = Wout;
			float *Xout1 = Xout; 
			float *Yout1 = Yout; 
			float *Zout1 = Zout; 
			float sina = sin(thisloc);
			float sinb = sin(thiselevation);		
			float cosa = cos(thisloc);
			float cosb = cos(thiselevation);
			if(rho >= 1) 
			    {sinint = (rsqrt2 * (sin(0.78539816339745 * 1.0))) / pow(rho, (float)1.5);
			    cosint =  (rsqrt2 * (cos(0.78539816339745 * 1.0))) / pow(rho, (float)1.5);}
			else 
			    {sinint = rsqrt2 * (sin(0.78539816339745 * rho));
			    cosint = rsqrt2 * (cos(0.78539816339745 * rho));};
			float W_amp = rsqrt2 * cosint;
			float X_amp = cosa * cosb * sinint;
			float Y_amp = sina * cosb * sinint;
			float Z_amp = sinb * sinint;
			float W_slope = CALCSLOPE(W_amp, W_amp);
			float X_slope = CALCSLOPE(X_amp, X_amp);
			float Y_slope = CALCSLOPE(Y_amp, Y_amp);
			float Z_slope = CALCSLOPE(Z_amp, Z_amp);
			int nsmps = sc_min(grainbf->counter, inNumSamples - i);
			for (int j=0; j<nsmps; ++j) {
			    float amp = y1 * y1;
			    float outval = amp * in[j];
			    // do the BFormat encoding on the resultant signal
			    Wout1[j] += outval * W_amp;
			    Xout1[j] += outval * X_amp;
			    Yout1[j] += outval * Y_amp;
			    Zout1[j] += outval * Z_amp;
			    W_amp += W_slope;
			    X_amp += X_slope;
			    Y_amp += Y_slope;
			    Z_amp += Z_slope;
			    double y0 = b1 * y1 - y2;
			    y2 = y1;
			    y1 = y0;
			}
			grainbf->y1 = y1;
			grainbf->y2 = y2;
			grainbf->loc = thisloc;
			grainbf->elevation = thiselevation;
			grainbf->m_W_amp = W_amp;
			grainbf->m_X_amp = X_amp;
			grainbf->m_Y_amp = Y_amp;
			grainbf->m_Z_amp = Z_amp;
			// store random values
			RPUT
			grainbf->counter -= nsmps;
			if (grainbf->counter <= 0) {
				// remove grain
				*grainbf = unit->mGrains[--unit->mNumActive];
			}
		}
	}	
	unit->m_NextGrain = nextGrain;
}

void MonoGrainBF_Ctor(MonoGrainBF *unit)
{	
	SETCALC(MonoGrainBF_next);
	unit->mNumActive = 0;
	unit->m_NextGrain = 1;
	ClearUnitOutputs(unit, 1);
}

void PVSynth_Ctor(AtsSynth* unit)
{
	SETCALC(AtsSynth_next);
	int tableSizeSin = ft->mSineSize;
	unit->m_numPartials = ZIN0(1);
	unit->m_partialStart = ZIN0(2);
	unit->m_partialSkip = ZIN0(3);
	int numPartials = (int)unit->m_numPartials;
	unit->m_lomask = (tableSizeSin - 1) << 3; 
	unit->m_radtoinc = tableSizeSin * (rtwopi * 65536.); 
	unit->m_cpstoinc = tableSizeSin * SAMPLEDUR * 65536.; 
	
	unit->m_phase = (int32*)RTAlloc(unit->mWorld, numPartials * sizeof(int32));
	
	unit->m_fbufnum = -1e9f;

	ClearUnitOutputs(unit, 1);

}

void PVSynth_Dtor(AtsSynth* unit)
{	
	RTFree(unit->mWorld, unit->m_phase);
}

void PVInfo_Ctor(AtsParInfo* unit)
{
	SETCALC(AtsParInfo_next);
	unit->m_partialNum = ZIN0(1);
	unit->m_framestart = ZIN0(2);
	unit->m_fbufnum = -1e9f;
	ClearUnitOutputs(unit, 1);

}


void Balance_Ctor(Balance* unit)
{
	unit->m_scaler = 0.f;
	if (INRATE(1) == calc_FullRate) {
	    SETCALC(Balance_next_za);
	    } else {
	    SETCALC(Balance_next_zk);
	    }
	unit->m_numPeriods = (int)IN0(2);
	unit->m_pastTestSig = 0.f;
	unit->m_pastTestSums = (float*)RTAlloc(unit->mWorld, unit->m_numPeriods * sizeof(float));
	unit->m_scaler = 0.f;
	// zero out the m_pastTestSums array
	for (int i = 0; i < unit->m_numPeriods; i++){
	    unit->m_pastTestSums[i] = 0.f;}
	ClearUnitOutputs(unit, 1);
	unit->m_numTurns = 0;
}

void Balance_Dtor(Balance* unit)
{	
	RTFree(unit->mWorld, unit->m_pastTestSums);
}


void Balance_next_za(Balance* unit, int inNumSamples)
{
	float* insig = IN(0);
	float* testsig = IN(1);
	float* out = OUT(0);
	float insum = 0.f;
	float testsum = 0.f;
	float scaler = unit->m_scaler;
	float newinrms, newtestrms, newscale, levelinc, sumholder;
	
	for(int i = 0; i < inNumSamples; i++){
	    insum += pow(insig[i], 2);
	    testsum += pow(testsig[i], 2);
	    }
	
	sumholder = testsum;
	
	// place the current value into 
	// work in the past samples
	for(int i = 0; i < unit->m_numTurns; i++){
	    testsum += unit->m_pastTestSums[i];
	    }
	    
	newinrms = sqrt(insum/inNumSamples);
	newtestrms = sqrt(testsum/(inNumSamples + (inNumSamples * unit->m_numTurns)));
	
	// protect a divide by zero... perhaps there is a better way to do this
	if(newinrms <= 0.f) {
	    newscale = 0.f; levelinc = 0.f;
	    } else {
	    newscale = newtestrms / newinrms;
	    levelinc = CALCSLOPE(newscale, scaler);
	    }
		
	for(int i = 0; i < inNumSamples; i++){
	    out[i] = insig[i] * scaler;
	    // overwrite old data with new data
	    scaler += levelinc;
	    }
	
	unit->m_scaler = newscale;
	
	unit->m_pastTestSums[unit->m_numTurns] = sumholder; 
	// numTurns will keep advancing the pointer into the buffer of past samples
	unit->m_numTurns = (unit->m_numTurns + 1);
	if(unit->m_numTurns >= (unit->m_numPeriods - 2)) SETCALC(Balance_next_a);
	
}
	
void Balance_next_a(Balance* unit, int inNumSamples)
{
	float* insig = IN(0);
	float* testsig = IN(1);
	float* out = OUT(0);
	float insum = 0.f;
	float testsum = 0.f;
	float scaler = unit->m_scaler;
	float newinrms, newtestrms, newscale, levelinc, sumholder;
	
	for(int i = 0; i < inNumSamples; i++){
	    insum += pow(insig[i], 2);
	    testsum += pow(testsig[i], 2);
	    }
	
	sumholder = testsum;
	
	// place the current value into 
	// work in the past samples
	for(int i = 0; i < (unit->m_numPeriods - 1); i++){
	    testsum += unit->m_pastTestSums[i];
	    }

	newinrms = sqrt(insum/inNumSamples);
	newtestrms = sqrt(testsum / (inNumSamples * unit->m_numPeriods));
		
	// protect a divide by zero... perhaps there is a better way to do this
	if(newinrms <= 0.f) {
	    newscale = 0.f; levelinc = 0.f;
	    } else {
	    newscale = newtestrms / newinrms;
	    levelinc = CALCSLOPE(newscale, scaler);
	    }
		
	for(int i = 0; i < inNumSamples; i++){
	    out[i] = insig[i] * scaler;
	    // overwrite old data with new data
	    scaler += levelinc;
	    }
	    
	unit->m_scaler = scaler;
	
//	Print("scaler : %6.6f\n", scaler);
	unit->m_pastTestSums[unit->m_numTurns] = sumholder; 
	// numTurns will keep advancing the pointer into the buffer of past samples
	unit->m_numTurns = (unit->m_numTurns + 1)%(unit->m_numPeriods);
	
}


void Balance_next_zk(Balance* unit, int inNumSamples)
{
	float* insig = IN(0);
	float testsig = IN0(1);
	float* out = OUT(0);
	float insum = 0.f;
	float testsum = 0.f;
	float scaler = unit->m_scaler;
	float newinrms, newtestrms, newscale, levelinc, sumholder;
	float lastTestSig = unit->m_pastTestSig;
	float testSigSlope = CALCSLOPE(testsig, lastTestSig);
	
	for(int i = 0; i < inNumSamples; i++){
		testsum += pow(lastTestSig, 2);
		insum += pow(insig[i], 2);
		lastTestSig += testSigSlope;
	    }
	
	sumholder = testsum;
	
	// place the current value into 
	// work in the past samples
	for(int i = 0; i < unit->m_numTurns; i++){
	    testsum += unit->m_pastTestSums[i];
	    }
	    
	newinrms = sqrt(insum/inNumSamples);
	newtestrms = sqrt(testsum/(inNumSamples + (inNumSamples * unit->m_numTurns)));
	
	// protect a divide by zero... perhaps there is a better way to do this
	if(newinrms <= 0.f) {
	    newscale = 0.f; levelinc = 0.f;
	    } else {
	    newscale = newtestrms / newinrms;
	    levelinc = CALCSLOPE(newscale, scaler);
	    }
		
	for(int i = 0; i < inNumSamples; i++){
	    out[i] = insig[i] * scaler;
	    // overwrite old data with new data
	    scaler += levelinc;
	    }
	
	unit->m_scaler = newscale;
	
	unit->m_pastTestSums[unit->m_numTurns] = sumholder; 
	// numTurns will keep advancing the pointer into the buffer of past samples
	unit->m_numTurns = (unit->m_numTurns + 1);
	if(unit->m_numTurns >= (unit->m_numPeriods - 2)) SETCALC(Balance_next_k);
	
}
	
void Balance_next_k(Balance* unit, int inNumSamples)
{
	float* insig = IN(0);
	float testsig = IN0(1);
	float* out = OUT(0);
	float insum = 0.f;
	float testsum = 0.f;
	float scaler = unit->m_scaler;
	float newinrms, newtestrms, newscale, levelinc, sumholder;
	
	float lastTestSig = unit->m_pastTestSig;
	float testSigSlope = CALCSLOPE(testsig, lastTestSig);
	
	for(int i = 0; i < inNumSamples; i++){
		testsum += pow(lastTestSig, 2);
		insum += pow(insig[i], 2);
		lastTestSig += testSigSlope;
	    }
	
	sumholder = testsum;
	
	// place the current value into 
	// work in the past samples
	for(int i = 0; i < (unit->m_numPeriods - 1); i++){
	    testsum += unit->m_pastTestSums[i];
	    }
	    
	newinrms = sqrt(insum/inNumSamples);
	newtestrms = sqrt(testsum/(inNumSamples * unit->m_numPeriods));
	
	// protect a divide by zero... perhaps there is a better way to do this
	if(newinrms <= 0.f) {
	    newscale = 0.f; levelinc = 0.f;
	    } else {
	    newscale = newtestrms / newinrms;
	    levelinc = CALCSLOPE(newscale, scaler);
	    }
		
	for(int i = 0; i < inNumSamples; i++){
	    out[i] = insig[i] * scaler;
	    // overwrite old data with new data
	    scaler += levelinc;
	    }
	unit->m_scaler = newscale;
	
//	Print("scaler : %6.6f\n", scaler);
	unit->m_pastTestSums[unit->m_numTurns] = sumholder; 
	// numTurns will keep advancing the pointer into the buffer of past samples
	unit->m_numTurns = (unit->m_numTurns + 1)%(unit->m_numPeriods - 1);
	
}



/*

moogvcf is a digital emulation of the Moog diode ladder filter configuration. This emulation is based loosely 
on the paper "Analyzing the Moog VCF with Considerations for Digital Implemnetation" by Stilson and Smith (CCRMA). 
This version was originally coded in Csound by Josep Comajuncosas. Some modifications and conversion to C were 
done by Hans Mikelson
SC version by Joshua Parmenter
Note: This filter requires that the input signal be normalized to one.
MoogVCF.ar(in, fco, res, mul, add); fco = filter cut-off, res = resonance (0 -> 1)

*/

void MoogVCF_Ctor(MoogVCF* unit) {
    if (INRATE(1) == calc_FullRate) {
	if (INRATE(2) == calc_FullRate) {
		SETCALC(MoogVCF_next_aa); // a-rate fco and res
		} else {
		SETCALC(MoogVCF_next_ak); // a-rate fco, k-rate res
		}
	    } else {
	if (INRATE(2) == calc_FullRate) {
		SETCALC(MoogVCF_next_ka); // k-rate fco, a-rate res
		} else {
		SETCALC(MoogVCF_next_kk);
		} // k-rate fco and res
	    }
	unit->m_fco = (IN0(1) * 2) / SAMPLERATE;
	unit->m_res = IN0(2);
	unit->m_xnm1 = unit->m_y1nm1 = unit->m_y2nm1 = unit->m_y3nm1 = 0.0f;
	unit->m_y1n = unit->m_y2n = unit->m_y3n = unit->m_y4n = 0.0f;
	ClearUnitOutputs(unit, 1);
}

void MoogVCF_next_kk(MoogVCF *unit, int inNumSamples){
    float* in = IN(0);
    float* out = OUT(0);
    float nextfco = IN0(1);
    float nextres = IN0(2);
    float fco = unit->m_fco; // already normalized
    float res = unit->m_res;
    float fcon = (nextfco * 2.0) / SAMPLERATE;   // filt freq, normalized to 0 to Nyquist

    float fcoslope = CALCSLOPE(fcon, fco);
    float resslope = CALCSLOPE(nextres, res);

    float xnm1 = unit->m_xnm1;
    float y1nm1 = unit->m_y1nm1;
    float y2nm1 = unit->m_y2nm1;
    float y3nm1 = unit->m_y3nm1;
    float y1n = unit->m_y1n;
    float y2n = unit->m_y2n;
    float y3n = unit->m_y3n;
    float y4n = unit->m_y4n;

    float kp, pp1d2, scale, xn, k;

    for (int i = 0; i < inNumSamples; i++) {
	kp = (3.6 * fco) - ((1.6 * fco) * fco) - 1.0;  /* Emperical tuning     */
	pp1d2 = (kp + 1.0) * 0.5;              /* Timesaver                  */
	scale = (float)exp((1.0 - (double)pp1d2 ) * 1.386249); /* Scaling factor  */
	k = res * scale;
	xn = in[i]; // make this similar to the CSound stuff for now... easier translation
	xn = xn - (k * y4n); /* Inverted feed back for corner peaking */
	
	/* Four cascaded onepole filters (bilinear transform) */
	y1n   = (xn  * pp1d2) + (xnm1  * pp1d2) - (kp * y1n);
	y2n   = (y1n * pp1d2) + (y1nm1 * pp1d2) - (kp * y2n);
	y3n   = (y2n * pp1d2) + (y2nm1 * pp1d2) - (kp * y3n);
	y4n   = (y3n * pp1d2) + (y3nm1 * pp1d2) - (kp * y4n);
	/* Clipper band limited sigmoid */
	y4n   = y4n - (((y4n * y4n) * y4n) / 6.0);
	xnm1  = xn;       /* Update Xn-1  */
	y1nm1 = y1n;      /* Update Y1n-1 */
	y2nm1 = y2n;      /* Update Y2n-1 */
	y3nm1 = y3n;      /* Update Y3n-1 */
	out[i] = y4n;
	fco += fcoslope;
	res += resslope;
    }

    unit->m_fco = fcon; // store the normalized frequency
    unit->m_res = nextres;
    unit->m_xnm1 = zapgremlins(xnm1);
    unit->m_y1nm1 = zapgremlins(y1nm1);
    unit->m_y2nm1 = zapgremlins(y2nm1);
    unit->m_y3nm1 = zapgremlins(y3nm1);
    unit->m_y1n = zapgremlins(y1n);
    unit->m_y2n = zapgremlins(y2n);
    unit->m_y3n = zapgremlins(y3n);
    unit->m_y4n = zapgremlins(y4n);
}

void MoogVCF_next_ka(MoogVCF *unit, int inNumSamples){
    float* in = IN(0);
    float* out = OUT(0);
    float nextfco = IN0(1);
    float* res = IN(2);
    float fco = unit->m_fco; // already normalized
    float fcon = (nextfco * 2.0) / SAMPLERATE;   // filt freq, normalized to 0 to Nyquist

    float fcoslope = CALCSLOPE(fcon, fco);

    float xnm1 = unit->m_xnm1;
    float y1nm1 = unit->m_y1nm1;
    float y2nm1 = unit->m_y2nm1;
    float y3nm1 = unit->m_y3nm1;
    float y1n = unit->m_y1n;
    float y2n = unit->m_y2n;
    float y3n = unit->m_y3n;
    float y4n = unit->m_y4n;

    float kp, pp1d2, scale, xn, k;

    for (int i = 0; i < inNumSamples; i++) {
	kp = (3.6 * fco) - ((1.6 * fco) * fco) - 1.0;  /* Emperical tuning     */
	pp1d2 = (kp + 1.0) * 0.5;              /* Timesaver                  */
	scale = (float)exp((1.0 - (double)pp1d2 ) * 1.386249); /* Scaling factor  */
	k = res[i] * scale;
	xn = in[i]; // make this similar to the CSound stuff for now... easier translation
	xn = xn - (k * y4n); /* Inverted feed back for corner peaking */
	
	/* Four cascaded onepole filters (bilinear transform) */
	y1n   = (xn  * pp1d2) + (xnm1  * pp1d2) - (kp * y1n);
	y2n   = (y1n * pp1d2) + (y1nm1 * pp1d2) - (kp * y2n);
	y3n   = (y2n * pp1d2) + (y2nm1 * pp1d2) - (kp * y3n);
	y4n   = (y3n * pp1d2) + (y3nm1 * pp1d2) - (kp * y4n);
	/* Clipper band limited sigmoid */
	y4n   = y4n - (((y4n * y4n) * y4n) / 6.0);
	xnm1  = xn;       /* Update Xn-1  */
	y1nm1 = y1n;      /* Update Y1n-1 */
	y2nm1 = y2n;      /* Update Y2n-1 */
	y3nm1 = y3n;      /* Update Y3n-1 */
	out[i] = y4n;
	fco += fcoslope;
    }

    unit->m_fco = fcon; // store the normalized frequency
    unit->m_xnm1 = zapgremlins(xnm1);
    unit->m_y1nm1 = zapgremlins(y1nm1);
    unit->m_y2nm1 = zapgremlins(y2nm1);
    unit->m_y3nm1 = zapgremlins(y3nm1);
    unit->m_y1n = zapgremlins(y1n);
    unit->m_y2n = zapgremlins(y2n);
    unit->m_y3n = zapgremlins(y3n);
    unit->m_y4n = zapgremlins(y4n);
}

void MoogVCF_next_ak(MoogVCF *unit, int inNumSamples){
    float* in = IN(0);
    float* out = OUT(0);
    float* fco = IN(1);
    float nextres = IN0(2);
    float res = unit->m_res;
    float fcon = 2.0 / SAMPLERATE;   // multiplier for filt freq, normalized to 0 to Nyquist

    float resslope = CALCSLOPE(nextres, res);

    float xnm1 = unit->m_xnm1;
    float y1nm1 = unit->m_y1nm1;
    float y2nm1 = unit->m_y2nm1;
    float y3nm1 = unit->m_y3nm1;
    float y1n = unit->m_y1n;
    float y2n = unit->m_y2n;
    float y3n = unit->m_y3n;
    float y4n = unit->m_y4n;

    float kp, pp1d2, scale, xn, k, thisfco;

    for (int i = 0; i < inNumSamples; i++) {
	thisfco = fco[i] * fcon;
	kp = (3.6 * thisfco) - ((1.6 * thisfco) * thisfco) - 1.0;  /* Emperical tuning     */
	pp1d2 = (kp + 1.0) * 0.5;              /* Timesaver                  */
	scale = (float)exp((1.0 - (double)pp1d2 ) * 1.386249); /* Scaling factor  */
	k = res * scale;
	xn = in[i]; // make this similar to the CSound stuff for now... easier translation
	xn = xn - (k * y4n); /* Inverted feed back for corner peaking */
	
	/* Four cascaded onepole filters (bilinear transform) */
	y1n   = (xn  * pp1d2) + (xnm1  * pp1d2) - (kp * y1n);
	y2n   = (y1n * pp1d2) + (y1nm1 * pp1d2) - (kp * y2n);
	y3n   = (y2n * pp1d2) + (y2nm1 * pp1d2) - (kp * y3n);
	y4n   = (y3n * pp1d2) + (y3nm1 * pp1d2) - (kp * y4n);
	/* Clipper band limited sigmoid */
	y4n   = y4n - (((y4n * y4n) * y4n) / 6.0);
	xnm1  = xn;       /* Update Xn-1  */
	y1nm1 = y1n;      /* Update Y1n-1 */
	y2nm1 = y2n;      /* Update Y2n-1 */
	y3nm1 = y3n;      /* Update Y3n-1 */
	out[i] = y4n;
	res += resslope;
    }

    unit->m_res = nextres;
    unit->m_xnm1 = zapgremlins(xnm1);
    unit->m_y1nm1 = zapgremlins(y1nm1);
    unit->m_y2nm1 = zapgremlins(y2nm1);
    unit->m_y3nm1 = zapgremlins(y3nm1);
    unit->m_y1n = zapgremlins(y1n);
    unit->m_y2n = zapgremlins(y2n);
    unit->m_y3n = zapgremlins(y3n);
    unit->m_y4n = zapgremlins(y4n);
}

void MoogVCF_next_aa(MoogVCF *unit, int inNumSamples){
    float* in = IN(0);
    float* out = OUT(0);
    float* fco = IN(1);
    float* res = IN(2);
    float fcon = 2.0 / SAMPLERATE;   // multiplier for filt freq, normalized to 0 to Nyquist

    float xnm1 = unit->m_xnm1;
    float y1nm1 = unit->m_y1nm1;
    float y2nm1 = unit->m_y2nm1;
    float y3nm1 = unit->m_y3nm1;
    float y1n = unit->m_y1n;
    float y2n = unit->m_y2n;
    float y3n = unit->m_y3n;
    float y4n = unit->m_y4n;

    float kp, pp1d2, scale, xn, k, thisfco;

    for (int i = 0; i < inNumSamples; i++) {
	thisfco = fco[i] * fcon;
	kp = (3.6 * thisfco) - ((1.6 * thisfco) * thisfco) - 1.0;  /* Emperical tuning     */
	pp1d2 = (kp + 1.0) * 0.5;              /* Timesaver                  */
	scale = (float)exp((1.0 - (double)pp1d2 ) * 1.386249); /* Scaling factor  */
	k = res[i] * scale;
	xn = in[i]; // make this similar to the CSound stuff for now... easier translation
	xn = xn - (k * y4n); /* Inverted feed back for corner peaking */
	
	/* Four cascaded onepole filters (bilinear transform) */
	y1n   = (xn  * pp1d2) + (xnm1  * pp1d2) - (kp * y1n);
	y2n   = (y1n * pp1d2) + (y1nm1 * pp1d2) - (kp * y2n);
	y3n   = (y2n * pp1d2) + (y2nm1 * pp1d2) - (kp * y3n);
	y4n   = (y3n * pp1d2) + (y3nm1 * pp1d2) - (kp * y4n);
	/* Clipper band limited sigmoid */
	y4n   = y4n - (((y4n * y4n) * y4n) / 6.0);
	xnm1  = xn;       /* Update Xn-1  */
	y1nm1 = y1n;      /* Update Y1n-1 */
	y2nm1 = y2n;      /* Update Y2n-1 */
	y3nm1 = y3n;      /* Update Y3n-1 */
	out[i] = y4n;
    }

    unit->m_xnm1 = zapgremlins(xnm1);
    unit->m_y1nm1 = zapgremlins(y1nm1);
    unit->m_y2nm1 = zapgremlins(y2nm1);
    unit->m_y3nm1 = zapgremlins(y3nm1);
    unit->m_y1n = zapgremlins(y1n);
    unit->m_y2n = zapgremlins(y2n);
    unit->m_y3n = zapgremlins(y3n);
    unit->m_y4n = zapgremlins(y4n);
}



struct IEnvGen : public Unit
{
    float m_level;
    float m_startpoint, m_numvals, m_pointin;
    float* m_envvals;
};


extern "C"
{
	void IEnvGen_next_a(IEnvGen *unit, int inNumSamples);
	void IEnvGen_next_k(IEnvGen *unit, int inNumSamples);
	void IEnvGen_Ctor(IEnvGen* unit);
	void IEnvGen_Dtor(IEnvGen* unit);
}

enum {
	shape_Step,
	shape_Linear,
	shape_Exponential,
	shape_Sine,
	shape_Welch,
	shape_Curve,
	shape_Squared,
	shape_Cubed
};

#define GET_ENV_VAL \
	switch (shape) \
	    { \
	    case shape_Step : \
		    level = unit->m_level = endLevel; \
		    break; \
	    case shape_Linear : \
	    default: \
		    level = unit->m_level = pos * (endLevel - begLevel) + begLevel; \
		    break; \
	    case shape_Exponential : \
		    level = unit->m_level = begLevel * pow(endLevel / begLevel, pos); \
		    break; \
	    case shape_Sine : \
		    level = unit->m_level = begLevel + (endLevel - begLevel) * (-cos(pi * pos) * 0.5 + 0.5); \
		    break; \
	    case shape_Welch : \
	    { \
		    if (begLevel < endLevel) \
			    level = unit->m_level = begLevel + (endLevel - begLevel) * sin(pi2 * pos); \
		    else \
			    level = unit->m_level = endLevel - (endLevel - begLevel) * sin(pi2 - pi2 * pos); \
		    break; \
	    } \
	    case shape_Curve : \
		    if (fabs((float)curve) < 0.0001) { \
			    level = unit->m_level = pos * (endLevel - begLevel) + begLevel; \
		    } else { \
			    double denom = 1. - exp((float)curve); \
			    double numer = 1. - exp((float)(pos * curve)); \
			    level = unit->m_level = begLevel + (endLevel - begLevel) * (numer/denom); \
		    } \
		    break; \
	    case shape_Squared : \
	    { \
		    double sqrtBegLevel = sqrt(begLevel); \
		    double sqrtEndLevel = sqrt(endLevel); \
		    double sqrtLevel = pos * (sqrtEndLevel - sqrtBegLevel) + sqrtBegLevel; \
		    level = unit->m_level = sqrtLevel * sqrtLevel; \
		    break; \
	    } \
	    case shape_Cubed : \
	    { \
		    double cbrtBegLevel = pow(begLevel, (float)0.3333333); \
		    double cbrtEndLevel = pow(endLevel, (float)0.3333333); \
		    double cbrtLevel = pos * (cbrtEndLevel - cbrtBegLevel) + cbrtBegLevel; \
		    level = unit->m_level = cbrtLevel * cbrtLevel * cbrtLevel; \
		    break; \
	    } \
	} \


void IEnvGen_Ctor(IEnvGen *unit)
{

	if (INRATE(0) == calc_FullRate) {
		SETCALC(IEnvGen_next_a);
		} else {
		SETCALC(IEnvGen_next_k); 
		}
		
	// pointer,
	// initlevel, numstages, totaldur,
	// [dur, shape, curve, level] * numvals
	int numstages = (int)IN0(2);
	int numvals = numstages * 4; // initlevel + (levels, dur, shape, curves) * stages
	float point = unit->m_pointin = IN0(0);
	unit->m_envvals = (float*)RTAlloc(unit->mWorld, (int)(numvals + 1.) * sizeof(float));
	
	unit->m_envvals[0] = IN0(1);
//	Print("%3.3f\n", unit->m_envvals[0]);
	// fill m_envvals with the values;
	for (int i = 1; i <= numvals; i++) {
	    unit->m_envvals[i] = IN0(3 + i);
//	    Print("%3.3f\n", unit->m_envvals[i]);
	}
	
//	float out = OUT0(0);
	int numStages = (int)IN0(2);
	float totalDur = IN0(3);
	float level = 0.f;
	float newtime = 0.f;
	int stage = 0;
	float seglen = 0.f;
	if (point >= totalDur) {
	    unit->m_level = level = unit->m_envvals[numStages * 4]; // grab the last value
	    } else {
	    if (point <= 0.0) {
		unit->m_level = level = unit->m_envvals[0];
		} else {
		float segpos = point;
		// determine which segment the current time pointer needs calculated
		for(int j = 0; point >= newtime; j++) {
			seglen = unit->m_envvals[(j * 4) + 1];
			newtime += seglen;
			segpos -= seglen;
			stage = j;
		    }
			
		segpos = segpos + seglen;
		float begLevel = unit->m_envvals[(stage * 4)];
		int shape = (int)unit->m_envvals[(stage * 4) + 2];
		int curve = (int)unit->m_envvals[(stage * 4) + 3];
		float endLevel = unit->m_envvals[(stage * 4) + 4];
		float pos = (segpos / seglen);
			
		GET_ENV_VAL
		}
	    }
	OUT0(0) = level;
}

void IEnvGen_Dtor(IEnvGen *unit)
{
	RTFree(unit->mWorld, unit->m_envvals);
}


void IEnvGen_next_a(IEnvGen *unit, int inNumSamples)
{
    float* out = OUT(0);
    float level = unit->m_level;
    float* pointin = IN(0);
    int numStages = (int)IN0(2);
    float point; // = unit->m_pointin;

    float totalDur = IN0(3);

    int stagemul;
    // pointer,
    // level0, numstages, totaldur,    
    // [initval, [dur, shape, curve, level] * N ]
	
    for( int i = 0; i < inNumSamples; i++) {
	    if (pointin[i] == unit->m_pointin){
		out[i] = level;
		} else {
		unit->m_pointin = point = pointin[i];
		float newtime = 0.f;
		int stage = 0;
		float seglen = 0.f;
		if (point >= totalDur) {
		    unit->m_level = level = unit->m_envvals[numStages * 4]; // grab the last value
		    } else {
		    if (point <= 0.0) {
			unit->m_level = level = unit->m_envvals[0];
			} else {
			float segpos = point;
			// determine which segment the current time pointer needs
			for(int j = 0; point >= newtime; j++) {
				seglen = unit->m_envvals[(j * 4) + 1];
				newtime += seglen;
				segpos -= seglen;
				stage = j;
			    }
			stagemul = stage * 4;
			segpos = segpos + seglen;
			float begLevel = unit->m_envvals[stagemul];
			int shape = (int)unit->m_envvals[stagemul + 2];
			int curve = (int)unit->m_envvals[stagemul + 3];
			float endLevel = unit->m_envvals[stagemul + 4];
			float pos = (segpos / seglen); 
			
			GET_ENV_VAL
		    }
		}
	    }
	out[i] = level;
	}
}


void IEnvGen_next_k(IEnvGen *unit, int inNumSamples)
{
    float* out = OUT(0);
    float level = unit->m_level;
    float pointin = IN0(0);
    int numStages = (int)IN0(2);

    float totalDur = IN0(3);
    // pointer,
    // level0, numstages, totaldur,
    float point = unit->m_pointin;

    // [initval, [dur, shape, curve, level] * N ]

    if (pointin == unit->m_pointin) {
	for(int i = 0; i < inNumSamples; i++){
	    out[i] = level;
	    }
	} else {
	float pointslope = CALCSLOPE(pointin, point);	
	for( int i = 0; i < inNumSamples; i++) {
		float newtime = 0.f;
		int stage = 0;
		float seglen = 0.f;
		if (point >= totalDur) {
		    unit->m_level = level = unit->m_envvals[numStages * 4]; // grab the last value
		    } else {
		    if (point <= 0.0) {
			unit->m_level = level = unit->m_envvals[0];
			} else {
			float segpos = point;
			// determine which segment the current time pointer needs calculated
			for(int j = 0; point >= newtime; j++) {
				seglen = unit->m_envvals[(j * 4) + 1];
				newtime += seglen;
				segpos -= seglen;
				stage = j;
			    }
			
			segpos = segpos + seglen;
			float begLevel = unit->m_envvals[(stage * 4)];
			int shape = (int)unit->m_envvals[(stage * 4) + 2];
			float curve = unit->m_envvals[(stage * 4) + 3];
			float endLevel = unit->m_envvals[(stage * 4) + 4];
			float pos = (segpos / seglen);
			
			GET_ENV_VAL
		    }
		}
	    out[i] = level;
	    point += pointslope;
	    }

   	unit->m_pointin = pointin;
    }
}
	
void PosRatio_Ctor(PosRatio* unit){
    float period = unit->period = IN0(1);
	unit->maxsamples = (int)(SAMPLERATE / period);
    float firstsample = unit->lastsample = IN0(0);
    SETCALC(PosRatio_next);
    if(firstsample > 0.f) unit->posvals = 1; else unit->posvals = 0;
	ClearUnitOutputs(unit, 1);
}



void PosRatio_next(PosRatio *unit, int inNumSamples){
    float *in = IN(0);
    float *out = OUT(0);
	float period = IN0(1); // period can only change AFTER the previous period is complete
	float thresh = IN0(2);
	float lastperiod = unit->lastperiod;
	int maxsamples = unit->maxsamples;
	 int posvals = unit->posvals;
	int numvals = unit->numvals;
	int testposvals = unit->testposvals;
	int testnumvals = unit->testnumvals; 
	float ratio = unit->ratio;
	float testratio = unit->testratio;
	float outval = unit->outval;
	float lastsample = unit->lastsample;
	float cursample;
	int i;

	
	// loop for inNumSamples
	// check for zero crossings...  and for a given period of pos values.
	// if there is a zero crossing, calc new ratio and IF it is within
	// given threshold, output that value
	// Forget interpolation for now between ratios
	
	for(i = 0; i < inNumSamples; i++){
		cursample = in[i];
		if(cursample > 0.) {posvals++; testposvals++;}
		numvals++;
		testnumvals++;
		// adjust the test values... if new ones are needed, also adjust period and maxsamples if needed
		if(testnumvals == maxsamples){
			testratio = testposvals / testnumvals;
			testposvals = testnumvals = 0;
			if(lastperiod != period){unit->maxsamples = maxsamples = (int)(SAMPLERATE/period); unit->lastperiod = period;}
			}
		// check for a zero crossing OR if maxvals has been reached
		if(((lastsample <= 0.) && (cursample > 0.f)) || (numvals == maxsamples)){
			ratio = (float)posvals / (float)numvals;
			posvals = numvals = 0; // zero out the bookkeeping
			if((testratio * (1. + thresh) > ratio) || ((testratio * (1. - thresh)) < ratio))
				outval = ratio;
				else
				outval = testratio;
			}

		lastsample = cursample;
		out[i] = outval;
	}
		unit->lastsample = lastsample;
		unit->posvals = posvals;
		unit->numvals = numvals;
		unit->testposvals = testposvals;
		unit->testnumvals = testnumvals;
		unit->ratio = ratio;
		unit->testratio = testratio;
		unit->outval = outval;
	}


void DelayUnit_AllocDelayLine(DelayUnit *unit)
{
	long delaybufsize = (long)ceil(unit->m_maxdelaytime * SAMPLERATE + 1.f);
	delaybufsize = delaybufsize + BUFLENGTH;
	delaybufsize = NEXTPOWEROFTWO(delaybufsize);  // round up to next power of two
	unit->m_fdelaylen = unit->m_idelaylen = delaybufsize;
	
	RTFree(unit->mWorld, unit->m_dlybuf);
	unit->m_dlybuf = (float*)RTAlloc(unit->mWorld, delaybufsize * sizeof(float));
	unit->m_mask = delaybufsize - 1;
}

float CalcDelay(DelayUnit *unit, float delaytime);
float CalcDelay(DelayUnit *unit, float delaytime)
{
	float next_dsamp = delaytime * SAMPLERATE;
	return sc_clip(next_dsamp, 1.f, unit->m_fdelaylen);
}

void DelayUnit_Reset(DelayUnit *unit)
{
	unit->m_maxdelaytime = ZIN0(2);
	unit->m_delaytime = ZIN0(3);
	unit->m_dlybuf = 0;
	
	DelayUnit_AllocDelayLine(unit);

	unit->m_dsamp = CalcDelay(unit, unit->m_delaytime);	
	
	unit->m_numoutput = 0;
	unit->m_iwrphase = 0;
}


void DelayUnit_Dtor(DelayUnit *unit)
{
	RTFree(unit->mWorld, unit->m_dlybuf);
}

inline float CalcFeedback(float delaytime, float decaytime)
{
	if (delaytime == 0.f) {
		return 0.f;
	} else if (decaytime > 0.f) {
		return exp(log001 * delaytime / decaytime);
	} else if (decaytime < 0.f) {
		return -exp(log001 * delaytime / -decaytime);
	} else {
		return 0.f;
	}
}

void FeedbackDelay_Reset(FeedbackDelay *unit)
{
	unit->m_decaytime = ZIN0(4);
	
	DelayUnit_Reset(unit);
	
	unit->m_feedbk = CalcFeedback(unit->m_delaytime, unit->m_decaytime);
}

////////////////////////////////////////////////////////////////////////////////////////////////////////
/*

*/
void Pluck_Ctor(Pluck *unit)
{
	FeedbackDelay_Reset(unit);
	if (INRATE(1) == calc_FullRate) {
	    if(INRATE(5) == calc_FullRate){
		SETCALC(Pluck_next_aa_z);
		} else {
		SETCALC(Pluck_next_ak_z); //ak
		}
	    } else {
	    if(INRATE(5) == calc_FullRate){ 
		SETCALC(Pluck_next_ka_z); //ka
		} else { 
		SETCALC(Pluck_next_kk_z); //kk
		}
	    }
	OUT0(0) = unit->m_lastsamp = 0.f;
	unit->m_prevtrig = 0.f;
	unit->m_inputsamps = 0;
	unit->m_coef = IN0(5);
}

void Pluck_next_aa(Pluck *unit, int inNumSamples)
{
	float *out = OUT(0);
	float *in = IN(0);
	float *trig = IN(1);
	float delaytime = IN0(3);
	float decaytime = IN0(4);
	float *coef = IN(5);
	float lastsamp = unit->m_lastsamp;
	unsigned long inputsamps = unit->m_inputsamps;
	
	float *dlybuf = unit->m_dlybuf;
	long iwrphase = unit->m_iwrphase;
	float dsamp = unit->m_dsamp;
	float feedbk = unit->m_feedbk;
	long mask = unit->m_mask;
	float thisin, curtrig;
	float prevtrig = unit->m_prevtrig;

	if (delaytime == unit->m_delaytime && decaytime == unit->m_decaytime) {
		long idsamp = (long)dsamp;
		float frac = dsamp - idsamp;
	for(int i = 0; i < inNumSamples; i++){
			curtrig = trig[i];
			if ((prevtrig <= 0.f) && (curtrig > 0.f)) {
			    inputsamps = (long)(delaytime * unit->mRate->mSampleRate + .5f);
			    }
			prevtrig = curtrig;
			long irdphase1 = iwrphase - idsamp;
			long irdphase2 = irdphase1 - 1;
			long irdphase3 = irdphase1 - 2;
			long irdphase0 = irdphase1 + 1;
			if (inputsamps > 0) {
			    thisin = in[i];
			    --inputsamps;
			    } else {
			    thisin = 0.f;
			    }
			float d0 = dlybuf[irdphase0 & mask];
			float d1 = dlybuf[irdphase1 & mask];
			float d2 = dlybuf[irdphase2 & mask];
			float d3 = dlybuf[irdphase3 & mask];
			float value = cubicinterp(frac, d0, d1, d2, d3);
			float thiscoef = coef[i];
			float onepole = ((1. - fabs(thiscoef)) * value) + (thiscoef * lastsamp);
			dlybuf[iwrphase & mask] = thisin + feedbk * onepole;				
			out[i] = lastsamp = onepole;
			iwrphase++;
		};
	} else {
	
		float next_dsamp = CalcDelay(unit, delaytime);
		float dsamp_slope = CALCSLOPE(next_dsamp, dsamp);
		
		float next_feedbk = CalcFeedback(delaytime, decaytime);
		float feedbk_slope = CALCSLOPE(next_feedbk, feedbk);

	    for(int i = 0; i < inNumSamples; i++){
			curtrig = trig[i];
			if ((prevtrig <= 0.f) && (curtrig > 0.f)) {
			    inputsamps = (long)(delaytime * unit->mRate->mSampleRate + .5f);
			    }
			prevtrig = curtrig;
			dsamp += dsamp_slope;
			long idsamp = (long)dsamp;
			float frac = dsamp - idsamp;
			long irdphase1 = iwrphase - idsamp;
			long irdphase2 = irdphase1 - 1;
			long irdphase3 = irdphase1 - 2;
			long irdphase0 = irdphase1 + 1;
			if (inputsamps > 0) {
			    thisin = in[i];
			    --inputsamps;
			    } else {
			    thisin = 0.f;
			    }
			float d0 = dlybuf[irdphase0 & mask];
			float d1 = dlybuf[irdphase1 & mask];
			float d2 = dlybuf[irdphase2 & mask];
			float d3 = dlybuf[irdphase3 & mask];
			float value = cubicinterp(frac, d0, d1, d2, d3);
			float thiscoef = coef[i];
			float onepole = ((1. - fabs(thiscoef)) * value) + (thiscoef * lastsamp);
			dlybuf[iwrphase & mask] = thisin + feedbk * onepole;				
			out[i] = lastsamp = onepole;
			feedbk += feedbk_slope;
			iwrphase++;
		};
		unit->m_feedbk = feedbk;
		unit->m_dsamp = dsamp;
		unit->m_delaytime = delaytime;
		unit->m_decaytime = decaytime;
	}
	
	unit->m_prevtrig = prevtrig;
	unit->m_inputsamps = inputsamps;
	unit->m_lastsamp = zapgremlins(lastsamp);
	unit->m_iwrphase = iwrphase;

}


void Pluck_next_aa_z(Pluck *unit, int inNumSamples)
{
	float *out = OUT(0);
	float *in = IN(0);
	float *trig = IN(1);
	float delaytime = IN0(3);
	float decaytime = IN0(4);
	float *coef = IN(5);
	float lastsamp = unit->m_lastsamp;
	
	float *dlybuf = unit->m_dlybuf;
	long iwrphase = unit->m_iwrphase;
	float dsamp = unit->m_dsamp;
	float feedbk = unit->m_feedbk;
	long mask = unit->m_mask;
	float d0, d1, d2, d3;
	float thisin, curtrig;
	unsigned long inputsamps = unit->m_inputsamps;
	float prevtrig = unit->m_prevtrig;
	
	if (delaytime == unit->m_delaytime && decaytime == unit->m_decaytime) {
		long idsamp = (long)dsamp;
		float frac = dsamp - idsamp;
		for(int i = 0; i < inNumSamples; i++){
			curtrig = trig[i];
			if ((prevtrig <= 0.f) && (curtrig > 0.f)) {
			    inputsamps = (long)(delaytime * unit->mRate->mSampleRate + .5f);
			    }
			prevtrig = curtrig;
			long irdphase1 = iwrphase - idsamp;
			long irdphase2 = irdphase1 - 1;
			long irdphase3 = irdphase1 - 2;
			long irdphase0 = irdphase1 + 1;
			if (inputsamps > 0) {
			    thisin = in[i];
			    --inputsamps;
			    } else {
			    thisin = 0.f;
			    }
			if (irdphase0 < 0) {
				dlybuf[iwrphase & mask] = thisin; 
				out[i] = 0.f;
			} else {
				if (irdphase1 < 0) {
					d1 = d2 = d3 = 0.f;
					d0 = dlybuf[irdphase0 & mask];
				} else if (irdphase2 < 0) {
					d1 = d2 = d3 = 0.f;
					d0 = dlybuf[irdphase0 & mask];
					d1 = dlybuf[irdphase1 & mask];
				} else if (irdphase3 < 0) {
					d3 = 0.f;
					d0 = dlybuf[irdphase0 & mask];
					d1 = dlybuf[irdphase1 & mask];
					d2 = dlybuf[irdphase2 & mask];
				} else {
					d0 = dlybuf[irdphase0 & mask];
					d1 = dlybuf[irdphase1 & mask];
					d2 = dlybuf[irdphase2 & mask];
					d3 = dlybuf[irdphase3 & mask];
				}
				float value = cubicinterp(frac, d0, d1, d2, d3);
			float thiscoef = coef[i];
			float onepole = ((1. - fabs(thiscoef)) * value) + (thiscoef * lastsamp);
			dlybuf[iwrphase & mask] = thisin + feedbk * onepole;				
			out[i] = lastsamp = onepole;
			}
			iwrphase++;
		};
	} else {
	
		float next_dsamp = CalcDelay(unit, delaytime);
		float dsamp_slope = CALCSLOPE(next_dsamp, dsamp);
		
		float next_feedbk = CalcFeedback(delaytime, decaytime);
		float feedbk_slope = CALCSLOPE(next_feedbk, feedbk);

		for(int i = 0; i < inNumSamples; i++) {
			curtrig = trig[i];
			if ((prevtrig <= 0.f) && (curtrig > 0.f)) {
			    inputsamps = (long)(delaytime * unit->mRate->mSampleRate + .5f);
			    }
			prevtrig = curtrig;
			dsamp += dsamp_slope;
			long idsamp = (long)dsamp;
			float frac = dsamp - idsamp;
			long irdphase1 = iwrphase - idsamp;
			long irdphase2 = irdphase1 - 1;
			long irdphase3 = irdphase1 - 2;
			long irdphase0 = irdphase1 + 1;
			if (inputsamps > 0) {
			    thisin = in[i];
			    --inputsamps;
			    } else {
			    thisin = 0.f;
			    }
			if (irdphase0 < 0) {
				dlybuf[iwrphase & mask] = thisin;
				out[i] = 0.f;
			    } else {
				if (irdphase1 < 0) {
					d1 = d2 = d3 = 0.f;
					d0 = dlybuf[irdphase0 & mask];
				} else if (irdphase2 < 0) {
					d1 = d2 = d3 = 0.f;
					d0 = dlybuf[irdphase0 & mask];
					d1 = dlybuf[irdphase1 & mask];
				} else if (irdphase3 < 0) {
					d3 = 0.f;
					d0 = dlybuf[irdphase0 & mask];
					d1 = dlybuf[irdphase1 & mask];
					d2 = dlybuf[irdphase2 & mask];
				} else {
					d0 = dlybuf[irdphase0 & mask];
					d1 = dlybuf[irdphase1 & mask];
					d2 = dlybuf[irdphase2 & mask];
					d3 = dlybuf[irdphase3 & mask];
				}
				float value = cubicinterp(frac, d0, d1, d2, d3);
			float thiscoef = coef[i];
			float onepole = ((1. - fabs(thiscoef)) * value) + (thiscoef * lastsamp);
			dlybuf[iwrphase & mask] = thisin + feedbk * onepole;				
			out[i] = lastsamp = onepole;
			}
			feedbk += feedbk_slope;
			iwrphase++;
		};
		unit->m_feedbk = feedbk;
		unit->m_dsamp = dsamp;
		unit->m_delaytime = delaytime;
		unit->m_decaytime = decaytime;
	}
	
	unit->m_inputsamps = inputsamps;
	unit->m_prevtrig = prevtrig;
	unit->m_lastsamp = zapgremlins(lastsamp);
	unit->m_iwrphase = iwrphase;

	unit->m_numoutput += inNumSamples;
	if (unit->m_numoutput >= unit->m_idelaylen) {
		SETCALC(Pluck_next_aa);
	}
}

void Pluck_next_kk(Pluck *unit, int inNumSamples)
{
	float *out = OUT(0);
	float *in = IN(0);
	float trig = IN0(1);
	float delaytime = IN0(3);
	float decaytime = IN0(4);
	float coef = IN0(5);
	float lastsamp = unit->m_lastsamp;
	unsigned long inputsamps = unit->m_inputsamps;
	
	float *dlybuf = unit->m_dlybuf;
	long iwrphase = unit->m_iwrphase;
	float dsamp = unit->m_dsamp;
	float feedbk = unit->m_feedbk;
	long mask = unit->m_mask;
	float thisin;
	
	if ((unit->m_prevtrig <= 0.f) && (trig > 0.f)) {
	    inputsamps = (long)(delaytime * unit->mRate->mSampleRate + .5f);
	    }
	unit->m_prevtrig = trig;

	if (delaytime == unit->m_delaytime && decaytime == unit->m_decaytime && coef == unit->m_coef) {
		long idsamp = (long)dsamp;
		float frac = dsamp - idsamp;

	    for(int i = 0; i < inNumSamples; i++){
			long irdphase1 = iwrphase - idsamp;
			long irdphase2 = irdphase1 - 1;
			long irdphase3 = irdphase1 - 2;
			long irdphase0 = irdphase1 + 1;
			if (inputsamps > 0) {
			    thisin = in[i];
			    --inputsamps;
			    } else {
			    thisin = 0.f;
			    }
			float d0 = dlybuf[irdphase0 & mask];
			float d1 = dlybuf[irdphase1 & mask];
			float d2 = dlybuf[irdphase2 & mask];
			float d3 = dlybuf[irdphase3 & mask];
			float value = cubicinterp(frac, d0, d1, d2, d3);
			float onepole = ((1. - fabs(coef)) * value) + (coef * lastsamp);
			dlybuf[iwrphase & mask] = thisin + (feedbk * onepole);				
			out[i] = lastsamp = onepole; //value;
			iwrphase++;
		};
	} else {
	
		float next_dsamp = CalcDelay(unit, delaytime);
		float dsamp_slope = CALCSLOPE(next_dsamp, dsamp);
		
		float next_feedbk = CalcFeedback(delaytime, decaytime);
		float feedbk_slope = CALCSLOPE(next_feedbk, feedbk);
		
		float curcoef = unit->m_coef;
		float coef_slope = CALCSLOPE(coef, curcoef);
		    
	    for(int i = 0; i < inNumSamples; i++){
			dsamp += dsamp_slope;
			long idsamp = (long)dsamp;
			float frac = dsamp - idsamp;
			long irdphase1 = iwrphase - idsamp;
			long irdphase2 = irdphase1 - 1;
			long irdphase3 = irdphase1 - 2;
			long irdphase0 = irdphase1 + 1;
			if (inputsamps > 0) {
			    thisin = in[i];
			    --inputsamps;
			    } else {
			    thisin = 0.f;
			    }
			float d0 = dlybuf[irdphase0 & mask];
			float d1 = dlybuf[irdphase1 & mask];
			float d2 = dlybuf[irdphase2 & mask];
			float d3 = dlybuf[irdphase3 & mask];
			float value = cubicinterp(frac, d0, d1, d2, d3);
			float onepole = ((1. - fabs(curcoef)) * value) + (curcoef * lastsamp);
			dlybuf[iwrphase & mask] = thisin + (feedbk * onepole);				
			out[i] = lastsamp = onepole; //value;
			feedbk += feedbk_slope;
			curcoef += coef_slope;
			iwrphase++;
		};
		unit->m_feedbk = feedbk;
		unit->m_coef = coef;
		unit->m_dsamp = dsamp;
		unit->m_delaytime = delaytime;
		unit->m_decaytime = decaytime;
	}
	
	unit->m_inputsamps = inputsamps;
	unit->m_lastsamp = zapgremlins(lastsamp);
	unit->m_iwrphase = iwrphase;

}


void Pluck_next_kk_z(Pluck *unit, int inNumSamples)
{
	float *out = OUT(0);
	float *in = IN(0);
	float trig = IN0(1);
	float delaytime = IN0(3);
	float decaytime = IN0(4);
	float coef = IN0(5);
	float lastsamp = unit->m_lastsamp;
	
	float *dlybuf = unit->m_dlybuf;
	long iwrphase = unit->m_iwrphase;
	float dsamp = unit->m_dsamp;
	float feedbk = unit->m_feedbk;
	long mask = unit->m_mask;
	float d0, d1, d2, d3;
	float thisin;
	unsigned long inputsamps = unit->m_inputsamps;

	if ((unit->m_prevtrig <= 0.f) && (trig > 0.f)) {
	    inputsamps = (long)(delaytime * unit->mRate->mSampleRate + .5f);
	    }
	unit->m_prevtrig = trig;
				
	if (delaytime == unit->m_delaytime && decaytime == unit->m_decaytime && coef == unit->m_coef) {
		long idsamp = (long)dsamp;
		float frac = dsamp - idsamp;
		for(int i = 0; i < inNumSamples; i++){

			long irdphase1 = iwrphase - idsamp;
			long irdphase2 = irdphase1 - 1;
			long irdphase3 = irdphase1 - 2;
			long irdphase0 = irdphase1 + 1;
			if (inputsamps > 0) {
			    thisin = in[i];
			    --inputsamps;
			    } else {
			    thisin = 0.f;
			    }
			if (irdphase0 < 0) {
				dlybuf[iwrphase & mask] = thisin; 
				out[i] = 0.f;
			} else {
				if (irdphase1 < 0) {
					d1 = d2 = d3 = 0.f;
					d0 = dlybuf[irdphase0 & mask];
				} else if (irdphase2 < 0) {
					d1 = d2 = d3 = 0.f;
					d0 = dlybuf[irdphase0 & mask];
					d1 = dlybuf[irdphase1 & mask];
				} else if (irdphase3 < 0) {
					d3 = 0.f;
					d0 = dlybuf[irdphase0 & mask];
					d1 = dlybuf[irdphase1 & mask];
					d2 = dlybuf[irdphase2 & mask];
				} else {
					d0 = dlybuf[irdphase0 & mask];
					d1 = dlybuf[irdphase1 & mask];
					d2 = dlybuf[irdphase2 & mask];
					d3 = dlybuf[irdphase3 & mask];
				}
				float value = cubicinterp(frac, d0, d1, d2, d3);
			float onepole = ((1. - fabs(coef)) * value) + (coef * lastsamp);
			dlybuf[iwrphase & mask] = thisin + (feedbk * onepole);				
			out[i] = lastsamp = onepole; //value;
			}
			iwrphase++;
		};
	} else {
	
		float next_dsamp = CalcDelay(unit, delaytime);
		float dsamp_slope = CALCSLOPE(next_dsamp, dsamp);
		
		float next_feedbk = CalcFeedback(delaytime, decaytime);
		float feedbk_slope = CALCSLOPE(next_feedbk, feedbk);
		
		float curcoef = unit->m_coef;
		float coef_slope = CALCSLOPE(coef, curcoef);

		for(int i = 0; i < inNumSamples; i++) {
			dsamp += dsamp_slope;
			long idsamp = (long)dsamp;
			float frac = dsamp - idsamp;
			long irdphase1 = iwrphase - idsamp;
			long irdphase2 = irdphase1 - 1;
			long irdphase3 = irdphase1 - 2;
			long irdphase0 = irdphase1 + 1;
			if (inputsamps > 0) {
			    thisin = in[i];
			    --inputsamps;
			    } else {
			    thisin = 0.f;
			    }
			if (irdphase0 < 0) {
				dlybuf[iwrphase & mask] = thisin;
				out[i] = 0.f;
			    } else {
				if (irdphase1 < 0) {
					d1 = d2 = d3 = 0.f;
					d0 = dlybuf[irdphase0 & mask];
				} else if (irdphase2 < 0) {
					d1 = d2 = d3 = 0.f;
					d0 = dlybuf[irdphase0 & mask];
					d1 = dlybuf[irdphase1 & mask];
				} else if (irdphase3 < 0) {
					d3 = 0.f;
					d0 = dlybuf[irdphase0 & mask];
					d1 = dlybuf[irdphase1 & mask];
					d2 = dlybuf[irdphase2 & mask];
				} else {
					d0 = dlybuf[irdphase0 & mask];
					d1 = dlybuf[irdphase1 & mask];
					d2 = dlybuf[irdphase2 & mask];
					d3 = dlybuf[irdphase3 & mask];
				}
			float value = cubicinterp(frac, d0, d1, d2, d3);
			float onepole = ((1. - fabs(curcoef)) * value) + (curcoef * lastsamp);
			dlybuf[iwrphase & mask] = thisin + (feedbk * onepole);				
			out[i] = lastsamp = onepole; //value;
			}
			feedbk += feedbk_slope;
			curcoef += coef_slope;
			iwrphase++;
		};
		unit->m_feedbk = feedbk;
		unit->m_dsamp = dsamp;
		unit->m_delaytime = delaytime;
		unit->m_decaytime = decaytime;
		unit->m_coef = coef;
	}
	
	unit->m_inputsamps = inputsamps;
	unit->m_lastsamp = zapgremlins(lastsamp);
	unit->m_iwrphase = iwrphase;

	unit->m_numoutput += inNumSamples;
	if (unit->m_numoutput >= unit->m_idelaylen) {
		SETCALC(Pluck_next_kk);
	}
}

void Pluck_next_ak(Pluck *unit, int inNumSamples)
{
	float *out = OUT(0);
	float *in = IN(0);
	float *trig = IN(1);
	float delaytime = IN0(3);
	float decaytime = IN0(4);
	float coef = IN0(5);
	float lastsamp = unit->m_lastsamp;
	unsigned long inputsamps = unit->m_inputsamps;
	
	float *dlybuf = unit->m_dlybuf;
	long iwrphase = unit->m_iwrphase;
	float dsamp = unit->m_dsamp;
	float feedbk = unit->m_feedbk;
	long mask = unit->m_mask;
	float thisin, curtrig;
	float prevtrig = unit->m_prevtrig;

	if (delaytime == unit->m_delaytime && decaytime == unit->m_decaytime) {
		long idsamp = (long)dsamp;
		float frac = dsamp - idsamp;
	for(int i = 0; i < inNumSamples; i++){
			curtrig = trig[i];
			if ((prevtrig <= 0.f) && (curtrig > 0.f)) {
			    inputsamps = (long)(delaytime * unit->mRate->mSampleRate + .5f);
			    }
			prevtrig = curtrig;
			long irdphase1 = iwrphase - idsamp;
			long irdphase2 = irdphase1 - 1;
			long irdphase3 = irdphase1 - 2;
			long irdphase0 = irdphase1 + 1;
			if (inputsamps > 0) {
			    thisin = in[i];
			    --inputsamps;
			    } else {
			    thisin = 0.f;
			    }
			float d0 = dlybuf[irdphase0 & mask];
			float d1 = dlybuf[irdphase1 & mask];
			float d2 = dlybuf[irdphase2 & mask];
			float d3 = dlybuf[irdphase3 & mask];
			float value = cubicinterp(frac, d0, d1, d2, d3);
			float onepole = ((1. - fabs(coef)) * value) + (coef * lastsamp);
			dlybuf[iwrphase & mask] = thisin + feedbk * onepole;				
			out[i] = lastsamp = onepole;
			iwrphase++;
		};
	} else {
	
		float next_dsamp = CalcDelay(unit, delaytime);
		float dsamp_slope = CALCSLOPE(next_dsamp, dsamp);
		
		float next_feedbk = CalcFeedback(delaytime, decaytime);
		float feedbk_slope = CALCSLOPE(next_feedbk, feedbk);

		float curcoef = unit->m_coef;
		float coef_slope = CALCSLOPE(coef, curcoef);

	    for(int i = 0; i < inNumSamples; i++){
			curtrig = trig[i];
			if ((prevtrig <= 0.f) && (curtrig > 0.f)) {
			    inputsamps = (long)(delaytime * unit->mRate->mSampleRate + .5f);
			    }
			prevtrig = curtrig;
			dsamp += dsamp_slope;
			long idsamp = (long)dsamp;
			float frac = dsamp - idsamp;
			long irdphase1 = iwrphase - idsamp;
			long irdphase2 = irdphase1 - 1;
			long irdphase3 = irdphase1 - 2;
			long irdphase0 = irdphase1 + 1;
			if (inputsamps > 0) {
			    thisin = in[i];
			    --inputsamps;
			    } else {
			    thisin = 0.f;
			    }
			float d0 = dlybuf[irdphase0 & mask];
			float d1 = dlybuf[irdphase1 & mask];
			float d2 = dlybuf[irdphase2 & mask];
			float d3 = dlybuf[irdphase3 & mask];
			float value = cubicinterp(frac, d0, d1, d2, d3);
			float onepole = ((1. - fabs(curcoef)) * value) + (curcoef * lastsamp);
			dlybuf[iwrphase & mask] = thisin + feedbk * onepole;				
			out[i] = lastsamp = onepole;
			feedbk += feedbk_slope;
			curcoef += coef_slope;
			iwrphase++;
		};
		unit->m_feedbk = feedbk;
		unit->m_dsamp = dsamp;
		unit->m_delaytime = delaytime;
		unit->m_decaytime = decaytime;
		unit->m_coef = coef;
	}
	
	unit->m_prevtrig = prevtrig;
	unit->m_inputsamps = inputsamps;
	unit->m_lastsamp = zapgremlins(lastsamp);
	unit->m_iwrphase = iwrphase;

}


void Pluck_next_ak_z(Pluck *unit, int inNumSamples)
{
	float *out = OUT(0);
	float *in = IN(0);
	float *trig = IN(1);
	float delaytime = IN0(3);
	float decaytime = IN0(4);
	float coef = IN0(5);
	float lastsamp = unit->m_lastsamp;
	
	float *dlybuf = unit->m_dlybuf;
	long iwrphase = unit->m_iwrphase;
	float dsamp = unit->m_dsamp;
	float feedbk = unit->m_feedbk;
	long mask = unit->m_mask;
	float d0, d1, d2, d3;
	float thisin, curtrig;
	unsigned long inputsamps = unit->m_inputsamps;
	float prevtrig = unit->m_prevtrig;
	
	if (delaytime == unit->m_delaytime && decaytime == unit->m_decaytime && coef == unit->m_coef) {
		long idsamp = (long)dsamp;
		float frac = dsamp - idsamp;
		for(int i = 0; i < inNumSamples; i++){
			curtrig = trig[i];
			if ((prevtrig <= 0.f) && (curtrig > 0.f)) {
			    inputsamps = (long)(delaytime * unit->mRate->mSampleRate + .5f);
			    }
			prevtrig = curtrig;
			long irdphase1 = iwrphase - idsamp;
			long irdphase2 = irdphase1 - 1;
			long irdphase3 = irdphase1 - 2;
			long irdphase0 = irdphase1 + 1;
			if (inputsamps > 0) {
			    thisin = in[i];
			    --inputsamps;
			    } else {
			    thisin = 0.f;
			    }
			if (irdphase0 < 0) {
				dlybuf[iwrphase & mask] = thisin; 
				out[i] = 0.f;
			} else {
				if (irdphase1 < 0) {
					d1 = d2 = d3 = 0.f;
					d0 = dlybuf[irdphase0 & mask];
				} else if (irdphase2 < 0) {
					d1 = d2 = d3 = 0.f;
					d0 = dlybuf[irdphase0 & mask];
					d1 = dlybuf[irdphase1 & mask];
				} else if (irdphase3 < 0) {
					d3 = 0.f;
					d0 = dlybuf[irdphase0 & mask];
					d1 = dlybuf[irdphase1 & mask];
					d2 = dlybuf[irdphase2 & mask];
				} else {
					d0 = dlybuf[irdphase0 & mask];
					d1 = dlybuf[irdphase1 & mask];
					d2 = dlybuf[irdphase2 & mask];
					d3 = dlybuf[irdphase3 & mask];
				}
				float value = cubicinterp(frac, d0, d1, d2, d3);
			float onepole = ((1. - fabs(coef)) * value) + (coef * lastsamp);
			dlybuf[iwrphase & mask] = thisin + feedbk * onepole;				
			out[i] = lastsamp = onepole;
			}
			iwrphase++;
		};
	} else {
	
		float next_dsamp = CalcDelay(unit, delaytime);
		float dsamp_slope = CALCSLOPE(next_dsamp, dsamp);
		
		float next_feedbk = CalcFeedback(delaytime, decaytime);
		float feedbk_slope = CALCSLOPE(next_feedbk, feedbk);
		
		float curcoef = unit->m_coef;
		float coef_slope = CALCSLOPE(coef, curcoef);
		
		for(int i = 0; i < inNumSamples; i++) {
			curtrig = trig[i];
			if ((prevtrig <= 0.f) && (curtrig > 0.f)) {
			    inputsamps = (long)(delaytime * unit->mRate->mSampleRate + .5f);
			    }
			prevtrig = curtrig;
			dsamp += dsamp_slope;
			long idsamp = (long)dsamp;
			float frac = dsamp - idsamp;
			long irdphase1 = iwrphase - idsamp;
			long irdphase2 = irdphase1 - 1;
			long irdphase3 = irdphase1 - 2;
			long irdphase0 = irdphase1 + 1;
			if (inputsamps > 0) {
			    thisin = in[i];
			    --inputsamps;
			    } else {
			    thisin = 0.f;
			    }
			if (irdphase0 < 0) {
				dlybuf[iwrphase & mask] = thisin;
				out[i] = 0.f;
			    } else {
				if (irdphase1 < 0) {
					d1 = d2 = d3 = 0.f;
					d0 = dlybuf[irdphase0 & mask];
				} else if (irdphase2 < 0) {
					d1 = d2 = d3 = 0.f;
					d0 = dlybuf[irdphase0 & mask];
					d1 = dlybuf[irdphase1 & mask];
				} else if (irdphase3 < 0) {
					d3 = 0.f;
					d0 = dlybuf[irdphase0 & mask];
					d1 = dlybuf[irdphase1 & mask];
					d2 = dlybuf[irdphase2 & mask];
				} else {
					d0 = dlybuf[irdphase0 & mask];
					d1 = dlybuf[irdphase1 & mask];
					d2 = dlybuf[irdphase2 & mask];
					d3 = dlybuf[irdphase3 & mask];
				}
				float value = cubicinterp(frac, d0, d1, d2, d3);
			float onepole = ((1. - fabs(curcoef)) * value) + (curcoef * lastsamp);
			dlybuf[iwrphase & mask] = thisin + feedbk * onepole;				
			out[i] = lastsamp = onepole;
			}
			feedbk += feedbk_slope;
			curcoef +=coef_slope;
			iwrphase++;
		};
		unit->m_feedbk = feedbk;
		unit->m_dsamp = dsamp;
		unit->m_delaytime = delaytime;
		unit->m_decaytime = decaytime;
		unit->m_coef = coef;
	}
	
	unit->m_inputsamps = inputsamps;
	unit->m_prevtrig = prevtrig;
	unit->m_lastsamp = zapgremlins(lastsamp);
	unit->m_iwrphase = iwrphase;

	unit->m_numoutput += inNumSamples;
	if (unit->m_numoutput >= unit->m_idelaylen) {
		SETCALC(Pluck_next_ak);
	}
}


void Pluck_next_ka(Pluck *unit, int inNumSamples)
{
	float *out = OUT(0);
	float *in = IN(0);
	float trig = IN0(1);
	float delaytime = IN0(3);
	float decaytime = IN0(4);
	float *coef = IN(5);
	float lastsamp = unit->m_lastsamp;
	unsigned long inputsamps = unit->m_inputsamps;
	
	float *dlybuf = unit->m_dlybuf;
	long iwrphase = unit->m_iwrphase;
	float dsamp = unit->m_dsamp;
	float feedbk = unit->m_feedbk;
	long mask = unit->m_mask;
	float thisin;

	if ((unit->m_prevtrig <= 0.f) && (trig > 0.f)) {
	    inputsamps = (long)(delaytime * unit->mRate->mSampleRate + .5f);
	    }
	unit->m_prevtrig = trig;
	
	if (delaytime == unit->m_delaytime && decaytime == unit->m_decaytime) {
		long idsamp = (long)dsamp;
		float frac = dsamp - idsamp;
	for(int i = 0; i < inNumSamples; i++){
			long irdphase1 = iwrphase - idsamp;
			long irdphase2 = irdphase1 - 1;
			long irdphase3 = irdphase1 - 2;
			long irdphase0 = irdphase1 + 1;
			if (inputsamps > 0) {
			    thisin = in[i];
			    --inputsamps;
			    } else {
			    thisin = 0.f;
			    }
			float d0 = dlybuf[irdphase0 & mask];
			float d1 = dlybuf[irdphase1 & mask];
			float d2 = dlybuf[irdphase2 & mask];
			float d3 = dlybuf[irdphase3 & mask];
			float value = cubicinterp(frac, d0, d1, d2, d3);
			float thiscoef = coef[i];
			float onepole = ((1. - fabs(thiscoef)) * value) + (thiscoef * lastsamp);
			dlybuf[iwrphase & mask] = thisin + feedbk * onepole;				
			out[i] = lastsamp = onepole;
			iwrphase++;
		};
	} else {
	
		float next_dsamp = CalcDelay(unit, delaytime);
		float dsamp_slope = CALCSLOPE(next_dsamp, dsamp);
		
		float next_feedbk = CalcFeedback(delaytime, decaytime);
		float feedbk_slope = CALCSLOPE(next_feedbk, feedbk);

	    for(int i = 0; i < inNumSamples; i++){
			dsamp += dsamp_slope;
			long idsamp = (long)dsamp;
			float frac = dsamp - idsamp;
			long irdphase1 = iwrphase - idsamp;
			long irdphase2 = irdphase1 - 1;
			long irdphase3 = irdphase1 - 2;
			long irdphase0 = irdphase1 + 1;
			if (inputsamps > 0) {
			    thisin = in[i];
			    --inputsamps;
			    } else {
			    thisin = 0.f;
			    }
			float d0 = dlybuf[irdphase0 & mask];
			float d1 = dlybuf[irdphase1 & mask];
			float d2 = dlybuf[irdphase2 & mask];
			float d3 = dlybuf[irdphase3 & mask];
			float value = cubicinterp(frac, d0, d1, d2, d3);
			float thiscoef = coef[i];
			float onepole = ((1. - fabs(thiscoef)) * value) + (thiscoef * lastsamp);
			dlybuf[iwrphase & mask] = thisin + feedbk * onepole;				
			out[i] = lastsamp = onepole;
			feedbk += feedbk_slope;
			iwrphase++;
		};
		unit->m_feedbk = feedbk;
		unit->m_dsamp = dsamp;
		unit->m_delaytime = delaytime;
		unit->m_decaytime = decaytime;
	}
	
	unit->m_inputsamps = inputsamps;
	unit->m_lastsamp = zapgremlins(lastsamp);
	unit->m_iwrphase = iwrphase;

}


void Pluck_next_ka_z(Pluck *unit, int inNumSamples)
{
	float *out = OUT(0);
	float *in = IN(0);
	float trig = IN0(1);
	float delaytime = IN0(3);
	float decaytime = IN0(4);
	float *coef = IN(5);
	float lastsamp = unit->m_lastsamp;
	
	float *dlybuf = unit->m_dlybuf;
	long iwrphase = unit->m_iwrphase;
	float dsamp = unit->m_dsamp;
	float feedbk = unit->m_feedbk;
	long mask = unit->m_mask;
	float d0, d1, d2, d3;
	float thisin;
	unsigned long inputsamps = unit->m_inputsamps;

	if ((unit->m_prevtrig <= 0.f) && (trig > 0.f)) {
	    inputsamps = (long)(delaytime * unit->mRate->mSampleRate + .5f);
	    }	
	
	unit->m_prevtrig = trig;
	
	if (delaytime == unit->m_delaytime && decaytime == unit->m_decaytime) {
		long idsamp = (long)dsamp;
		float frac = dsamp - idsamp;
		for(int i = 0; i < inNumSamples; i++){
			long irdphase1 = iwrphase - idsamp;
			long irdphase2 = irdphase1 - 1;
			long irdphase3 = irdphase1 - 2;
			long irdphase0 = irdphase1 + 1;
			if (inputsamps > 0) {
			    thisin = in[i];
			    --inputsamps;
			    } else {
			    thisin = 0.f;
			    }
			if (irdphase0 < 0) {
				dlybuf[iwrphase & mask] = thisin; 
				out[i] = 0.f;
			} else {
				if (irdphase1 < 0) {
					d1 = d2 = d3 = 0.f;
					d0 = dlybuf[irdphase0 & mask];
				} else if (irdphase2 < 0) {
					d1 = d2 = d3 = 0.f;
					d0 = dlybuf[irdphase0 & mask];
					d1 = dlybuf[irdphase1 & mask];
				} else if (irdphase3 < 0) {
					d3 = 0.f;
					d0 = dlybuf[irdphase0 & mask];
					d1 = dlybuf[irdphase1 & mask];
					d2 = dlybuf[irdphase2 & mask];
				} else {
					d0 = dlybuf[irdphase0 & mask];
					d1 = dlybuf[irdphase1 & mask];
					d2 = dlybuf[irdphase2 & mask];
					d3 = dlybuf[irdphase3 & mask];
				}
				float value = cubicinterp(frac, d0, d1, d2, d3);
			float thiscoef = coef[i];
			float onepole = ((1. - fabs(thiscoef)) * value) + (thiscoef * lastsamp);
			dlybuf[iwrphase & mask] = thisin + feedbk * onepole;				
			out[i] = lastsamp = onepole;
			}
			iwrphase++;
		};
	} else {
	
		float next_dsamp = CalcDelay(unit, delaytime);
		float dsamp_slope = CALCSLOPE(next_dsamp, dsamp);
		
		float next_feedbk = CalcFeedback(delaytime, decaytime);
		float feedbk_slope = CALCSLOPE(next_feedbk, feedbk);

		for(int i = 0; i < inNumSamples; i++) {
			dsamp += dsamp_slope;
			long idsamp = (long)dsamp;
			float frac = dsamp - idsamp;
			long irdphase1 = iwrphase - idsamp;
			long irdphase2 = irdphase1 - 1;
			long irdphase3 = irdphase1 - 2;
			long irdphase0 = irdphase1 + 1;
			if (inputsamps > 0) {
			    thisin = in[i];
			    --inputsamps;
			    } else {
			    thisin = 0.f;
			    }
			if (irdphase0 < 0) {
				dlybuf[iwrphase & mask] = thisin;
				out[i] = 0.f;
			    } else {
				if (irdphase1 < 0) {
					d1 = d2 = d3 = 0.f;
					d0 = dlybuf[irdphase0 & mask];
				} else if (irdphase2 < 0) {
					d1 = d2 = d3 = 0.f;
					d0 = dlybuf[irdphase0 & mask];
					d1 = dlybuf[irdphase1 & mask];
				} else if (irdphase3 < 0) {
					d3 = 0.f;
					d0 = dlybuf[irdphase0 & mask];
					d1 = dlybuf[irdphase1 & mask];
					d2 = dlybuf[irdphase2 & mask];
				} else {
					d0 = dlybuf[irdphase0 & mask];
					d1 = dlybuf[irdphase1 & mask];
					d2 = dlybuf[irdphase2 & mask];
					d3 = dlybuf[irdphase3 & mask];
				}
				float value = cubicinterp(frac, d0, d1, d2, d3);
			float thiscoef = coef[i];
			float onepole = ((1. - fabs(thiscoef)) * value) + (thiscoef * lastsamp);
			dlybuf[iwrphase & mask] = thisin + feedbk * onepole;				
			out[i] = lastsamp = onepole;
			}
			feedbk += feedbk_slope;
			iwrphase++;
		};
		unit->m_feedbk = feedbk;
		unit->m_dsamp = dsamp;
		unit->m_delaytime = delaytime;
		unit->m_decaytime = decaytime;
	}
	
	unit->m_inputsamps = inputsamps;
	unit->m_lastsamp = zapgremlins(lastsamp);
	unit->m_iwrphase = iwrphase;

	unit->m_numoutput += inNumSamples;
	if (unit->m_numoutput >= unit->m_idelaylen) {
		SETCALC(Pluck_next_ka);
	}
}

/* CombLP... gated input to delay line */

void CombLP_Ctor(CombLP *unit)
{
	FeedbackDelay_Reset(unit);
	if (INRATE(1) == calc_FullRate) {
	    if(INRATE(5) == calc_FullRate){
		SETCALC(CombLP_next_aa_z);
		} else {
		SETCALC(CombLP_next_ak_z); //ak
		}
	    } else {
	    if(INRATE(5) == calc_FullRate){ 
		SETCALC(CombLP_next_ka_z); //ka
		} else { 
		SETCALC(CombLP_next_kk_z); //kk
		}
	    }
	OUT0(0) = unit->m_lastsamp = 0.f;
	unit->m_prevtrig = 0.f;
	unit->m_coef = IN0(5);
}

void CombLP_next_aa(CombLP *unit, int inNumSamples)
{
	float *out = OUT(0);
	float *in = IN(0);
	float *gate = IN(1);
	float delaytime = IN0(3);
	float decaytime = IN0(4);
	float *coef = IN(5);
	float lastsamp = unit->m_lastsamp;
	
	float *dlybuf = unit->m_dlybuf;
	long iwrphase = unit->m_iwrphase;
	float dsamp = unit->m_dsamp;
	float feedbk = unit->m_feedbk;
	long mask = unit->m_mask;
	float thisin;
	
	if (delaytime == unit->m_delaytime && decaytime == unit->m_decaytime) {
		long idsamp = (long)dsamp;
		float frac = dsamp - idsamp;
	for(int i = 0; i < inNumSamples; i++){
			long irdphase1 = iwrphase - idsamp;
			long irdphase2 = irdphase1 - 1;
			long irdphase3 = irdphase1 - 2;
			long irdphase0 = irdphase1 + 1;
			thisin = in[i] * gate[i];
			float d0 = dlybuf[irdphase0 & mask];
			float d1 = dlybuf[irdphase1 & mask];
			float d2 = dlybuf[irdphase2 & mask];
			float d3 = dlybuf[irdphase3 & mask];
			float value = cubicinterp(frac, d0, d1, d2, d3);
			float thiscoef = coef[i];
			float onepole = ((1. - fabs(thiscoef)) * value) + (thiscoef * lastsamp);
			dlybuf[iwrphase & mask] = thisin + feedbk * onepole;				
			out[i] = lastsamp = onepole;
			iwrphase++;
		};
	} else {
	
		float next_dsamp = CalcDelay(unit, delaytime);
		float dsamp_slope = CALCSLOPE(next_dsamp, dsamp);
		
		float next_feedbk = CalcFeedback(delaytime, decaytime);
		float feedbk_slope = CALCSLOPE(next_feedbk, feedbk);

	    for(int i = 0; i < inNumSamples; i++){
			dsamp += dsamp_slope;
			long idsamp = (long)dsamp;
			float frac = dsamp - idsamp;
			long irdphase1 = iwrphase - idsamp;
			long irdphase2 = irdphase1 - 1;
			long irdphase3 = irdphase1 - 2;
			long irdphase0 = irdphase1 + 1;
			thisin = in[i] * gate[i];
			float d0 = dlybuf[irdphase0 & mask];
			float d1 = dlybuf[irdphase1 & mask];
			float d2 = dlybuf[irdphase2 & mask];
			float d3 = dlybuf[irdphase3 & mask];
			float value = cubicinterp(frac, d0, d1, d2, d3);
			float thiscoef = coef[i];
			float onepole = ((1. - fabs(thiscoef)) * value) + (thiscoef * lastsamp);
			dlybuf[iwrphase & mask] = thisin + feedbk * onepole;				
			out[i] = lastsamp = onepole;
			feedbk += feedbk_slope;
			iwrphase++;
		};
		unit->m_feedbk = feedbk;
		unit->m_dsamp = dsamp;
		unit->m_delaytime = delaytime;
		unit->m_decaytime = decaytime;
	}
	
	unit->m_lastsamp = zapgremlins(lastsamp);
	unit->m_iwrphase = iwrphase;

}


void CombLP_next_aa_z(CombLP *unit, int inNumSamples)
{
	float *out = OUT(0);
	float *in = IN(0);
	float *gate = IN(1);
	float delaytime = IN0(3);
	float decaytime = IN0(4);
	float *coef = IN(5);
	float lastsamp = unit->m_lastsamp;
	
	float *dlybuf = unit->m_dlybuf;
	long iwrphase = unit->m_iwrphase;
	float dsamp = unit->m_dsamp;
	float feedbk = unit->m_feedbk;
	long mask = unit->m_mask;
	float d0, d1, d2, d3;
	float thisin;
	
	if (delaytime == unit->m_delaytime && decaytime == unit->m_decaytime) {
		long idsamp = (long)dsamp;
		float frac = dsamp - idsamp;
		for(int i = 0; i < inNumSamples; i++) {
			long irdphase1 = iwrphase - idsamp;
			long irdphase2 = irdphase1 - 1;
			long irdphase3 = irdphase1 - 2;
			long irdphase0 = irdphase1 + 1;
			thisin = in[i] * gate[i];
			if (irdphase0 < 0) {
				dlybuf[iwrphase & mask] = thisin; 
				out[i] = 0.f;
			} else {
				if (irdphase1 < 0) {
					d1 = d2 = d3 = 0.f;
					d0 = dlybuf[irdphase0 & mask];
				} else if (irdphase2 < 0) {
					d1 = d2 = d3 = 0.f;
					d0 = dlybuf[irdphase0 & mask];
					d1 = dlybuf[irdphase1 & mask];
				} else if (irdphase3 < 0) {
					d3 = 0.f;
					d0 = dlybuf[irdphase0 & mask];
					d1 = dlybuf[irdphase1 & mask];
					d2 = dlybuf[irdphase2 & mask];
				} else {
					d0 = dlybuf[irdphase0 & mask];
					d1 = dlybuf[irdphase1 & mask];
					d2 = dlybuf[irdphase2 & mask];
					d3 = dlybuf[irdphase3 & mask];
				}
				float value = cubicinterp(frac, d0, d1, d2, d3);
			float thiscoef = coef[i];
			float onepole = ((1. - fabs(thiscoef)) * value) + (thiscoef * lastsamp);
			dlybuf[iwrphase & mask] = thisin + feedbk * onepole;				
			out[i] = lastsamp = onepole;
			}
			iwrphase++;
		};
	} else {
	
		float next_dsamp = CalcDelay(unit, delaytime);
		float dsamp_slope = CALCSLOPE(next_dsamp, dsamp);
		
		float next_feedbk = CalcFeedback(delaytime, decaytime);
		float feedbk_slope = CALCSLOPE(next_feedbk, feedbk);

		for(int i = 0; i < inNumSamples; i++) {
			dsamp += dsamp_slope;
			long idsamp = (long)dsamp;
			float frac = dsamp - idsamp;
			long irdphase1 = iwrphase - idsamp;
			long irdphase2 = irdphase1 - 1;
			long irdphase3 = irdphase1 - 2;
			long irdphase0 = irdphase1 + 1;
			thisin = in[i] * gate[i];
			if (irdphase0 < 0) {
				dlybuf[iwrphase & mask] = thisin;
				out[i] = 0.f;
			    } else {
				if (irdphase1 < 0) {
					d1 = d2 = d3 = 0.f;
					d0 = dlybuf[irdphase0 & mask];
				} else if (irdphase2 < 0) {
					d1 = d2 = d3 = 0.f;
					d0 = dlybuf[irdphase0 & mask];
					d1 = dlybuf[irdphase1 & mask];
				} else if (irdphase3 < 0) {
					d3 = 0.f;
					d0 = dlybuf[irdphase0 & mask];
					d1 = dlybuf[irdphase1 & mask];
					d2 = dlybuf[irdphase2 & mask];
				} else {
					d0 = dlybuf[irdphase0 & mask];
					d1 = dlybuf[irdphase1 & mask];
					d2 = dlybuf[irdphase2 & mask];
					d3 = dlybuf[irdphase3 & mask];
				}
				float value = cubicinterp(frac, d0, d1, d2, d3);
			float thiscoef = coef[i];
			float onepole = ((1. - fabs(thiscoef)) * value) + (thiscoef * lastsamp);
			dlybuf[iwrphase & mask] = thisin + feedbk * onepole;				
			out[i] = lastsamp = onepole;
			}
			feedbk += feedbk_slope;
			iwrphase++;
		};
		unit->m_feedbk = feedbk;
		unit->m_dsamp = dsamp;
		unit->m_delaytime = delaytime;
		unit->m_decaytime = decaytime;
	}
	
	unit->m_lastsamp = zapgremlins(lastsamp);
	unit->m_iwrphase = iwrphase;

	unit->m_numoutput += inNumSamples;
	if (unit->m_numoutput >= unit->m_idelaylen) {
		SETCALC(CombLP_next_aa);
	}
}

void CombLP_next_kk(CombLP *unit, int inNumSamples)
{
	float *out = OUT(0);
	float *in = IN(0);
	float gate = IN0(1);
	float delaytime = IN0(3);
	float decaytime = IN0(4);
	float coef = IN0(5);
	float lastsamp = unit->m_lastsamp;
	
	float *dlybuf = unit->m_dlybuf;
	long iwrphase = unit->m_iwrphase;
	float dsamp = unit->m_dsamp;
	float feedbk = unit->m_feedbk;
	long mask = unit->m_mask;
	float thisin;

	if (delaytime == unit->m_delaytime && decaytime == unit->m_decaytime && coef == unit->m_coef) {
		long idsamp = (long)dsamp;
		float frac = dsamp - idsamp;

	    for(int i = 0; i < inNumSamples; i++){
			long irdphase1 = iwrphase - idsamp;
			long irdphase2 = irdphase1 - 1;
			long irdphase3 = irdphase1 - 2;
			long irdphase0 = irdphase1 + 1;
			thisin = in[i] * gate;
			float d0 = dlybuf[irdphase0 & mask];
			float d1 = dlybuf[irdphase1 & mask];
			float d2 = dlybuf[irdphase2 & mask];
			float d3 = dlybuf[irdphase3 & mask];
			float value = cubicinterp(frac, d0, d1, d2, d3);
			float onepole = ((1. - fabs(coef)) * value) + (coef * lastsamp);
			dlybuf[iwrphase & mask] = thisin + (feedbk * onepole);				
			out[i] = lastsamp = onepole; //value;
			iwrphase++;
		};
	} else {
	
		float next_dsamp = CalcDelay(unit, delaytime);
		float dsamp_slope = CALCSLOPE(next_dsamp, dsamp);
		
		float next_feedbk = CalcFeedback(delaytime, decaytime);
		float feedbk_slope = CALCSLOPE(next_feedbk, feedbk);
		
		float curcoef = unit->m_coef;
		float coef_slope = CALCSLOPE(coef, curcoef);
		    
	    for(int i = 0; i < inNumSamples; i++){
			dsamp += dsamp_slope;
			long idsamp = (long)dsamp;
			float frac = dsamp - idsamp;
			long irdphase1 = iwrphase - idsamp;
			long irdphase2 = irdphase1 - 1;
			long irdphase3 = irdphase1 - 2;
			long irdphase0 = irdphase1 + 1;
			thisin = in[i] * gate;
			
			float d0 = dlybuf[irdphase0 & mask];
			float d1 = dlybuf[irdphase1 & mask];
			float d2 = dlybuf[irdphase2 & mask];
			float d3 = dlybuf[irdphase3 & mask];
			float value = cubicinterp(frac, d0, d1, d2, d3);
			float onepole = ((1. - fabs(curcoef)) * value) + (curcoef * lastsamp);
			dlybuf[iwrphase & mask] = thisin + (feedbk * onepole);				
			out[i] = lastsamp = onepole; //value;
			feedbk += feedbk_slope;
			curcoef += coef_slope;
			iwrphase++;
		};
		unit->m_feedbk = feedbk;
		unit->m_coef = coef;
		unit->m_dsamp = dsamp;
		unit->m_delaytime = delaytime;
		unit->m_decaytime = decaytime;
	}
	
	unit->m_lastsamp = zapgremlins(lastsamp);
	unit->m_iwrphase = iwrphase;

}


void CombLP_next_kk_z(CombLP *unit, int inNumSamples)
{
	float *out = OUT(0);
	float *in = IN(0);
	float gate = IN0(1);
	float delaytime = IN0(3);
	float decaytime = IN0(4);
	float coef = IN0(5);
	float lastsamp = unit->m_lastsamp;
	
	float *dlybuf = unit->m_dlybuf;
	long iwrphase = unit->m_iwrphase;
	float dsamp = unit->m_dsamp;
	float feedbk = unit->m_feedbk;
	long mask = unit->m_mask;
	float d0, d1, d2, d3;
	float thisin;
				
	if (delaytime == unit->m_delaytime && decaytime == unit->m_decaytime && coef == unit->m_coef) {
		long idsamp = (long)dsamp;
		float frac = dsamp - idsamp;
		for(int i = 0; i < inNumSamples; i++){

			long irdphase1 = iwrphase - idsamp;
			long irdphase2 = irdphase1 - 1;
			long irdphase3 = irdphase1 - 2;
			long irdphase0 = irdphase1 + 1;
			thisin = in[i] * gate;
			if (irdphase0 < 0) {
				dlybuf[iwrphase & mask] = thisin; 
				out[i] = 0.f;
			} else {
				if (irdphase1 < 0) {
					d1 = d2 = d3 = 0.f;
					d0 = dlybuf[irdphase0 & mask];
				} else if (irdphase2 < 0) {
					d1 = d2 = d3 = 0.f;
					d0 = dlybuf[irdphase0 & mask];
					d1 = dlybuf[irdphase1 & mask];
				} else if (irdphase3 < 0) {
					d3 = 0.f;
					d0 = dlybuf[irdphase0 & mask];
					d1 = dlybuf[irdphase1 & mask];
					d2 = dlybuf[irdphase2 & mask];
				} else {
					d0 = dlybuf[irdphase0 & mask];
					d1 = dlybuf[irdphase1 & mask];
					d2 = dlybuf[irdphase2 & mask];
					d3 = dlybuf[irdphase3 & mask];
				}
				float value = cubicinterp(frac, d0, d1, d2, d3);
			float onepole = ((1. - fabs(coef)) * value) + (coef * lastsamp);
			dlybuf[iwrphase & mask] = thisin + (feedbk * onepole);				
			out[i] = lastsamp = onepole; //value;
			}
			iwrphase++;
		};
	} else {
	
		float next_dsamp = CalcDelay(unit, delaytime);
		float dsamp_slope = CALCSLOPE(next_dsamp, dsamp);
		
		float next_feedbk = CalcFeedback(delaytime, decaytime);
		float feedbk_slope = CALCSLOPE(next_feedbk, feedbk);
		
		float curcoef = unit->m_coef;
		float coef_slope = CALCSLOPE(coef, curcoef);

		for(int i = 0; i < inNumSamples; i++) {
			dsamp += dsamp_slope;
			long idsamp = (long)dsamp;
			float frac = dsamp - idsamp;
			long irdphase1 = iwrphase - idsamp;
			long irdphase2 = irdphase1 - 1;
			long irdphase3 = irdphase1 - 2;
			long irdphase0 = irdphase1 + 1;
			thisin = in[i] * gate;
			if (irdphase0 < 0) {
				dlybuf[iwrphase & mask] = thisin;
				out[i] = 0.f;
			    } else {
				if (irdphase1 < 0) {
					d1 = d2 = d3 = 0.f;
					d0 = dlybuf[irdphase0 & mask];
				} else if (irdphase2 < 0) {
					d1 = d2 = d3 = 0.f;
					d0 = dlybuf[irdphase0 & mask];
					d1 = dlybuf[irdphase1 & mask];
				} else if (irdphase3 < 0) {
					d3 = 0.f;
					d0 = dlybuf[irdphase0 & mask];
					d1 = dlybuf[irdphase1 & mask];
					d2 = dlybuf[irdphase2 & mask];
				} else {
					d0 = dlybuf[irdphase0 & mask];
					d1 = dlybuf[irdphase1 & mask];
					d2 = dlybuf[irdphase2 & mask];
					d3 = dlybuf[irdphase3 & mask];
				}
			float value = cubicinterp(frac, d0, d1, d2, d3);
			float onepole = ((1. - fabs(curcoef)) * value) + (curcoef * lastsamp);
			dlybuf[iwrphase & mask] = thisin + (feedbk * onepole);				
			out[i] = lastsamp = onepole; //value;
			}
			feedbk += feedbk_slope;
			curcoef += coef_slope;
			iwrphase++;
		};
		unit->m_feedbk = feedbk;
		unit->m_dsamp = dsamp;
		unit->m_delaytime = delaytime;
		unit->m_decaytime = decaytime;
		unit->m_coef = coef;
	}
	
	unit->m_lastsamp = zapgremlins(lastsamp);
	unit->m_iwrphase = iwrphase;

	unit->m_numoutput += inNumSamples;
	if (unit->m_numoutput >= unit->m_idelaylen) {
		SETCALC(CombLP_next_kk);
	}
}

void CombLP_next_ak(CombLP *unit, int inNumSamples)
{
	float *out = OUT(0);
	float *in = IN(0);
	float *gate = IN(1);
	float delaytime = IN0(3);
	float decaytime = IN0(4);
	float coef = IN0(5);
	float lastsamp = unit->m_lastsamp;
	
	float *dlybuf = unit->m_dlybuf;
	long iwrphase = unit->m_iwrphase;
	float dsamp = unit->m_dsamp;
	float feedbk = unit->m_feedbk;
	long mask = unit->m_mask;
	float thisin;

	if (delaytime == unit->m_delaytime && decaytime == unit->m_decaytime) {
		long idsamp = (long)dsamp;
		float frac = dsamp - idsamp;
	for(int i = 0; i < inNumSamples; i++){
			long irdphase1 = iwrphase - idsamp;
			long irdphase2 = irdphase1 - 1;
			long irdphase3 = irdphase1 - 2;
			long irdphase0 = irdphase1 + 1;
			thisin = in[i] * gate[i];
			float d0 = dlybuf[irdphase0 & mask];
			float d1 = dlybuf[irdphase1 & mask];
			float d2 = dlybuf[irdphase2 & mask];
			float d3 = dlybuf[irdphase3 & mask];
			float value = cubicinterp(frac, d0, d1, d2, d3);
			float onepole = ((1. - fabs(coef)) * value) + (coef * lastsamp);
			dlybuf[iwrphase & mask] = thisin + feedbk * onepole;				
			out[i] = lastsamp = onepole;
			iwrphase++;
		};
	} else {
	
		float next_dsamp = CalcDelay(unit, delaytime);
		float dsamp_slope = CALCSLOPE(next_dsamp, dsamp);
		
		float next_feedbk = CalcFeedback(delaytime, decaytime);
		float feedbk_slope = CALCSLOPE(next_feedbk, feedbk);

		float curcoef = unit->m_coef;
		float coef_slope = CALCSLOPE(coef, curcoef);

	    for(int i = 0; i < inNumSamples; i++){
			dsamp += dsamp_slope;
			long idsamp = (long)dsamp;
			float frac = dsamp - idsamp;
			long irdphase1 = iwrphase - idsamp;
			long irdphase2 = irdphase1 - 1;
			long irdphase3 = irdphase1 - 2;
			long irdphase0 = irdphase1 + 1;
			thisin = in[i] * gate[i];
			float d0 = dlybuf[irdphase0 & mask];
			float d1 = dlybuf[irdphase1 & mask];
			float d2 = dlybuf[irdphase2 & mask];
			float d3 = dlybuf[irdphase3 & mask];
			float value = cubicinterp(frac, d0, d1, d2, d3);
			float onepole = ((1. - fabs(curcoef)) * value) + (curcoef * lastsamp);
			dlybuf[iwrphase & mask] = thisin + feedbk * onepole;				
			out[i] = lastsamp = onepole;
			feedbk += feedbk_slope;
			curcoef += coef_slope;
			iwrphase++;
		};
		unit->m_feedbk = feedbk;
		unit->m_dsamp = dsamp;
		unit->m_delaytime = delaytime;
		unit->m_decaytime = decaytime;
		unit->m_coef = coef;
	}
	
	unit->m_lastsamp = zapgremlins(lastsamp);
	unit->m_iwrphase = iwrphase;

}


void CombLP_next_ak_z(CombLP *unit, int inNumSamples)
{
	float *out = OUT(0);
	float *in = IN(0);
	float *gate = IN(1);
	float delaytime = IN0(3);
	float decaytime = IN0(4);
	float coef = IN0(5);
	float lastsamp = unit->m_lastsamp;
	
	float *dlybuf = unit->m_dlybuf;
	long iwrphase = unit->m_iwrphase;
	float dsamp = unit->m_dsamp;
	float feedbk = unit->m_feedbk;
	long mask = unit->m_mask;
	float d0, d1, d2, d3;
	float thisin;
	
	if (delaytime == unit->m_delaytime && decaytime == unit->m_decaytime && coef == unit->m_coef) {
		long idsamp = (long)dsamp;
		float frac = dsamp - idsamp;
		for(int i = 0; i < inNumSamples; i++){
			long irdphase1 = iwrphase - idsamp;
			long irdphase2 = irdphase1 - 1;
			long irdphase3 = irdphase1 - 2;
			long irdphase0 = irdphase1 + 1;
			thisin = in[i] * gate[i];
			if (irdphase0 < 0) {
				dlybuf[iwrphase & mask] = thisin; 
				out[i] = 0.f;
			} else {
				if (irdphase1 < 0) {
					d1 = d2 = d3 = 0.f;
					d0 = dlybuf[irdphase0 & mask];
				} else if (irdphase2 < 0) {
					d1 = d2 = d3 = 0.f;
					d0 = dlybuf[irdphase0 & mask];
					d1 = dlybuf[irdphase1 & mask];
				} else if (irdphase3 < 0) {
					d3 = 0.f;
					d0 = dlybuf[irdphase0 & mask];
					d1 = dlybuf[irdphase1 & mask];
					d2 = dlybuf[irdphase2 & mask];
				} else {
					d0 = dlybuf[irdphase0 & mask];
					d1 = dlybuf[irdphase1 & mask];
					d2 = dlybuf[irdphase2 & mask];
					d3 = dlybuf[irdphase3 & mask];
				}
				float value = cubicinterp(frac, d0, d1, d2, d3);
			float onepole = ((1. - fabs(coef)) * value) + (coef * lastsamp);
			dlybuf[iwrphase & mask] = thisin + feedbk * onepole;				
			out[i] = lastsamp = onepole;
			}
			iwrphase++;
		};
	} else {
	
		float next_dsamp = CalcDelay(unit, delaytime);
		float dsamp_slope = CALCSLOPE(next_dsamp, dsamp);
		
		float next_feedbk = CalcFeedback(delaytime, decaytime);
		float feedbk_slope = CALCSLOPE(next_feedbk, feedbk);
		
		float curcoef = unit->m_coef;
		float coef_slope = CALCSLOPE(coef, curcoef);
		
		for(int i = 0; i < inNumSamples; i++) {
			dsamp += dsamp_slope;
			long idsamp = (long)dsamp;
			float frac = dsamp - idsamp;
			long irdphase1 = iwrphase - idsamp;
			long irdphase2 = irdphase1 - 1;
			long irdphase3 = irdphase1 - 2;
			long irdphase0 = irdphase1 + 1;
			thisin = in[i] * gate[i];
			if (irdphase0 < 0) {
				dlybuf[iwrphase & mask] = thisin;
				out[i] = 0.f;
			    } else {
				if (irdphase1 < 0) {
					d1 = d2 = d3 = 0.f;
					d0 = dlybuf[irdphase0 & mask];
				} else if (irdphase2 < 0) {
					d1 = d2 = d3 = 0.f;
					d0 = dlybuf[irdphase0 & mask];
					d1 = dlybuf[irdphase1 & mask];
				} else if (irdphase3 < 0) {
					d3 = 0.f;
					d0 = dlybuf[irdphase0 & mask];
					d1 = dlybuf[irdphase1 & mask];
					d2 = dlybuf[irdphase2 & mask];
				} else {
					d0 = dlybuf[irdphase0 & mask];
					d1 = dlybuf[irdphase1 & mask];
					d2 = dlybuf[irdphase2 & mask];
					d3 = dlybuf[irdphase3 & mask];
				}
				float value = cubicinterp(frac, d0, d1, d2, d3);
			float onepole = ((1. - fabs(curcoef)) * value) + (curcoef * lastsamp);
			dlybuf[iwrphase & mask] = thisin + feedbk * onepole;				
			out[i] = lastsamp = onepole;
			}
			feedbk += feedbk_slope;
			curcoef +=coef_slope;
			iwrphase++;
		};
		unit->m_feedbk = feedbk;
		unit->m_dsamp = dsamp;
		unit->m_delaytime = delaytime;
		unit->m_decaytime = decaytime;
		unit->m_coef = coef;
	}
	
	unit->m_lastsamp = zapgremlins(lastsamp);
	unit->m_iwrphase = iwrphase;

	unit->m_numoutput += inNumSamples;
	if (unit->m_numoutput >= unit->m_idelaylen) {
		SETCALC(CombLP_next_ak);
	}
}


void CombLP_next_ka(CombLP *unit, int inNumSamples)
{
	float *out = OUT(0);
	float *in = IN(0);
	float gate = IN0(1);
	float delaytime = IN0(3);
	float decaytime = IN0(4);
	float *coef = IN(5);
	float lastsamp = unit->m_lastsamp;
	
	float *dlybuf = unit->m_dlybuf;
	long iwrphase = unit->m_iwrphase;
	float dsamp = unit->m_dsamp;
	float feedbk = unit->m_feedbk;
	long mask = unit->m_mask;
	float thisin;

	
	if (delaytime == unit->m_delaytime && decaytime == unit->m_decaytime) {
		long idsamp = (long)dsamp;
		float frac = dsamp - idsamp;
	for(int i = 0; i < inNumSamples; i++){
			long irdphase1 = iwrphase - idsamp;
			long irdphase2 = irdphase1 - 1;
			long irdphase3 = irdphase1 - 2;
			long irdphase0 = irdphase1 + 1;
			thisin = in[i] * gate;
			float d0 = dlybuf[irdphase0 & mask];
			float d1 = dlybuf[irdphase1 & mask];
			float d2 = dlybuf[irdphase2 & mask];
			float d3 = dlybuf[irdphase3 & mask];
			float value = cubicinterp(frac, d0, d1, d2, d3);
			float thiscoef = coef[i];
			float onepole = ((1. - fabs(thiscoef)) * value) + (thiscoef * lastsamp);
			dlybuf[iwrphase & mask] = thisin + feedbk * onepole;				
			out[i] = lastsamp = onepole;
			iwrphase++;
		};
	} else {
	
		float next_dsamp = CalcDelay(unit, delaytime);
		float dsamp_slope = CALCSLOPE(next_dsamp, dsamp);
		
		float next_feedbk = CalcFeedback(delaytime, decaytime);
		float feedbk_slope = CALCSLOPE(next_feedbk, feedbk);

	    for(int i = 0; i < inNumSamples; i++){
			dsamp += dsamp_slope;
			long idsamp = (long)dsamp;
			float frac = dsamp - idsamp;
			long irdphase1 = iwrphase - idsamp;
			long irdphase2 = irdphase1 - 1;
			long irdphase3 = irdphase1 - 2;
			long irdphase0 = irdphase1 + 1;
			thisin = in[i] * gate;
			float d0 = dlybuf[irdphase0 & mask];
			float d1 = dlybuf[irdphase1 & mask];
			float d2 = dlybuf[irdphase2 & mask];
			float d3 = dlybuf[irdphase3 & mask];
			float value = cubicinterp(frac, d0, d1, d2, d3);
			float thiscoef = coef[i];
			float onepole = ((1. - fabs(thiscoef)) * value) + (thiscoef * lastsamp);
			dlybuf[iwrphase & mask] = thisin + feedbk * onepole;				
			out[i] = lastsamp = onepole;
			feedbk += feedbk_slope;
			iwrphase++;
		};
		unit->m_feedbk = feedbk;
		unit->m_dsamp = dsamp;
		unit->m_delaytime = delaytime;
		unit->m_decaytime = decaytime;
	}
	
	unit->m_lastsamp = zapgremlins(lastsamp);
	unit->m_iwrphase = iwrphase;

}


void CombLP_next_ka_z(CombLP *unit, int inNumSamples)
{
	float *out = OUT(0);
	float *in = IN(0);
	float gate = IN0(1);
	float delaytime = IN0(3);
	float decaytime = IN0(4);
	float *coef = IN(5);
	float lastsamp = unit->m_lastsamp;
	
	float *dlybuf = unit->m_dlybuf;
	long iwrphase = unit->m_iwrphase;
	float dsamp = unit->m_dsamp;
	float feedbk = unit->m_feedbk;
	long mask = unit->m_mask;
	float d0, d1, d2, d3;
	float thisin;
	
	if (delaytime == unit->m_delaytime && decaytime == unit->m_decaytime) {
		long idsamp = (long)dsamp;
		float frac = dsamp - idsamp;
		for(int i = 0; i < inNumSamples; i++){
			long irdphase1 = iwrphase - idsamp;
			long irdphase2 = irdphase1 - 1;
			long irdphase3 = irdphase1 - 2;
			long irdphase0 = irdphase1 + 1;
			thisin = in[i] * gate;
			if (irdphase0 < 0) {
				dlybuf[iwrphase & mask] = thisin; 
				out[i] = 0.f;
			} else {
				if (irdphase1 < 0) {
					d1 = d2 = d3 = 0.f;
					d0 = dlybuf[irdphase0 & mask];
				} else if (irdphase2 < 0) {
					d1 = d2 = d3 = 0.f;
					d0 = dlybuf[irdphase0 & mask];
					d1 = dlybuf[irdphase1 & mask];
				} else if (irdphase3 < 0) {
					d3 = 0.f;
					d0 = dlybuf[irdphase0 & mask];
					d1 = dlybuf[irdphase1 & mask];
					d2 = dlybuf[irdphase2 & mask];
				} else {
					d0 = dlybuf[irdphase0 & mask];
					d1 = dlybuf[irdphase1 & mask];
					d2 = dlybuf[irdphase2 & mask];
					d3 = dlybuf[irdphase3 & mask];
				}
				float value = cubicinterp(frac, d0, d1, d2, d3);
			float thiscoef = coef[i];
			float onepole = ((1. - fabs(thiscoef)) * value) + (thiscoef * lastsamp);
			dlybuf[iwrphase & mask] = thisin + feedbk * onepole;				
			out[i] = lastsamp = onepole;
			}
			iwrphase++;
		};
	} else {
	
		float next_dsamp = CalcDelay(unit, delaytime);
		float dsamp_slope = CALCSLOPE(next_dsamp, dsamp);
		
		float next_feedbk = CalcFeedback(delaytime, decaytime);
		float feedbk_slope = CALCSLOPE(next_feedbk, feedbk);

		for(int i = 0; i < inNumSamples; i++) {
			dsamp += dsamp_slope;
			long idsamp = (long)dsamp;
			float frac = dsamp - idsamp;
			long irdphase1 = iwrphase - idsamp;
			long irdphase2 = irdphase1 - 1;
			long irdphase3 = irdphase1 - 2;
			long irdphase0 = irdphase1 + 1;
			thisin = in[i] * gate;
			if (irdphase0 < 0) {
				dlybuf[iwrphase & mask] = thisin;
				out[i] = 0.f;
			    } else {
				if (irdphase1 < 0) {
					d1 = d2 = d3 = 0.f;
					d0 = dlybuf[irdphase0 & mask];
				} else if (irdphase2 < 0) {
					d1 = d2 = d3 = 0.f;
					d0 = dlybuf[irdphase0 & mask];
					d1 = dlybuf[irdphase1 & mask];
				} else if (irdphase3 < 0) {
					d3 = 0.f;
					d0 = dlybuf[irdphase0 & mask];
					d1 = dlybuf[irdphase1 & mask];
					d2 = dlybuf[irdphase2 & mask];
				} else {
					d0 = dlybuf[irdphase0 & mask];
					d1 = dlybuf[irdphase1 & mask];
					d2 = dlybuf[irdphase2 & mask];
					d3 = dlybuf[irdphase3 & mask];
				}
				float value = cubicinterp(frac, d0, d1, d2, d3);
			float thiscoef = coef[i];
			float onepole = ((1. - fabs(thiscoef)) * value) + (thiscoef * lastsamp);
			dlybuf[iwrphase & mask] = thisin + feedbk * onepole;				
			out[i] = lastsamp = onepole;
			}
			feedbk += feedbk_slope;
			iwrphase++;
		};
		unit->m_feedbk = feedbk;
		unit->m_dsamp = dsamp;
		unit->m_delaytime = delaytime;
		unit->m_decaytime = decaytime;
	}
	
	unit->m_lastsamp = zapgremlins(lastsamp);
	unit->m_iwrphase = iwrphase;

	unit->m_numoutput += inNumSamples;
	if (unit->m_numoutput >= unit->m_idelaylen) {
		SETCALC(CombLP_next_ka);
	}
}
////////////////////////////////////////////////////////////////////////////////////////////////////////

void load(InterfaceTable *inTable)
{
	ft = inTable;
	
	DefineDtorUnit(AtsSynth);
	DefineDtorUnit(AtsNoiSynth);
	DefineSimpleUnit(AtsFreq);
	DefineSimpleUnit(AtsAmp);
	DefineSimpleUnit(AtsNoise);
	DefineSimpleUnit(AtsParInfo);
	DefineSimpleUnit(SinTone);
	DefineSimpleUnit(BFEncode1);
	DefineSimpleUnit(BFEncode2);
//	DefineSimpleUnit(BFEnocdeSter);
	DefineSimpleCantAliasUnit(BFDecode1);
	DefineSimpleUnit(BFManipulate);
	DefineSimpleUnit(B2Ster);
	DefineSimpleCantAliasUnit(Maxamp);
	DefineSimpleCantAliasUnit(Warp1);
	DefineSimpleCantAliasUnit(MonoGrain);
	DefineSimpleCantAliasUnit(MonoGrainBF);
	DefineSimpleUnit(LPCSynth);
	DefineSimpleUnit(LPCVals);
	DefineSimpleUnit(AudioMSG);
	DefineSimpleUnit(A2B);
	DefineSimpleUnit(B2A);
	DefineSimpleUnit(UHJ2B);
	DefineSimpleUnit(B2UHJ);
	DefineDtorUnit(PVSynth);
	DefineSimpleUnit(PVInfo);
	DefineDtorCantAliasUnit(Balance);
	DefineSimpleUnit(MoogVCF);
	DefineDtorUnit(IEnvGen);
	DefineSimpleUnit(PosRatio);
    
#define DefineDelayUnit(name) \
    (*ft->fDefineUnit)(#name, sizeof(name), (UnitCtorFunc)&name##_Ctor, \
    (UnitDtorFunc)&DelayUnit_Dtor, 0);
    
	DefineDelayUnit(Pluck);
	DefineDelayUnit(CombLP);
}


	