// Multiband Frequency Shifter
// written by Matthias Wagner March 2015, for MAT 240B

// Takes an audio signal and splits it into 4 bands via Linwitz-Riley filters
// (implementation from Robert Bristow-Johnson @ http://www.musicdsp.org/files/Audio-EQ-Cookbook.txt),
// and applies a different amount of frequency shift to each band
// (Hilbert filter implementation from Sean M. Costello @
// https://github.com/csound/csound/blob/95925a2af8973a0594d826a259ce603d3c651df4/Opcodes/ugsc.c).
// Big thanks to Joseph Tilbian and Andrés Cabrera for help and advice!!!

#include <stdio.h>
#include <portaudio.h>
#include <math.h>
#include <sndfile.h>

#include "Hilbert.h"
#include "Linkwitz_Riley.h"

#define SR 44100
#define TABLE_SIZE 4096
#define FRAMES_PER_BUFFER 256
#define DUR 10.0

#ifndef M_PI
#define M_PI 3.14159265
#endif

typedef struct
{
    float sine[TABLE_SIZE];
    float cosine[TABLE_SIZE];
    float sinPhase[4];
    float sinPhaseIncr[4];
    float cosPhase[4];
    float cosPhaseIncr[4];

    // split signal in LP and HP
    LR_coefficients coef_LP_pre, coef_HP_pre;
    LR_state LP_pre_s1, LP_pre_s2, HP_pre_s1, HP_pre_s2;
    // split HP from original signal again into HP and LP
    LR_coefficients coef_LP_hi_mid, coef_HP_hi;
    LR_state LP_hi_mid_s1, LP_hi_mid_s2, HP_hi_s1, HP_hi_s2;
    // split LP from original signal again into HP and LP
    LR_coefficients coef_LP_lo, coef_HP_lo_mid;
    LR_state LP_lo_s1, LP_lo_s2, HP_lo_mid_s1, HP_lo_mid_s2;

    // Hilbert coefficients
    H_coefficients coef_Hil1, coef_Hil2, coef_Hil3, coef_Hil4;

    double * bufferPointer;
    SNDFILE *sndFile;
    SF_INFO sfInfo;
    int position;
} paData;

static int paCallback( const void *inputBuffer,
                       void *outputBuffer,
                           unsigned long framesPerBuffer,
                           const PaStreamCallbackTimeInfo* timeInfo,
                           PaStreamCallbackFlags statusFlags,
                           void *userData ) // only points to starting point of memory
{
    paData *data = (paData*)userData; // recast for memory mapping
    float *in = (float*)inputBuffer;
    float *out = (float*)outputBuffer;
    float pre_lo, pre_hi, lo, lo_mid, hi_mid, hi;
    float lo_sine, lo_mid_sine, hi_mid_sine, hi_sine;
    float lo_cosine, lo_mid_cosine, hi_mid_cosine, hi_cosine;
    float lo_out, lo_mid_out, hi_mid_out, hi_out;
    float randMax[4];
    float func_out;
    unsigned int i, j;
    float input;
    (void) inputBuffer; /* Prevent unused variable warning. */

    for (i=0; i<framesPerBuffer; i++) {

        // increment modulators
        for (j=0; j<4; j++) {
            data->sinPhase[j] += data->sinPhaseIncr[j];
            while (data->sinPhase[j] >= TABLE_SIZE) { data->sinPhase[j] -= TABLE_SIZE; }
            while (data->sinPhase[j] < 0) { data->sinPhase[j] += TABLE_SIZE; }
            data->cosPhase[j] += data->cosPhaseIncr[j];
            while (data->cosPhase[j] >= TABLE_SIZE) { data->cosPhase[j] -= TABLE_SIZE; }
            while (data->cosPhase[j] < 0) { data->cosPhase[j] += TABLE_SIZE; }
//        data->sinPhaseIncr[j] += 0.00005;
//        data->cosPhaseIncr[j] += 0.00005;
        }

//        input = data->sine[(int)data->sinPhase[0]];
        input = *in++;

       // Linkwitz_Riley Test
//     // split signal into 2 parts
//     pre_lo = LP_LRFilter(input, &data->coef_LP_pre, &data->LP_pre_s1, &data->LP_pre_s2);
//     pre_hi = HP_LRFilter(input, &data->coef_HP_pre, &data->HP_pre_s1, &data->HP_pre_s2);
//     // split signal further into lo
//     lo = LP_LRFilter(pre_lo, &data->coef_LP_lo, &data->LP_lo_s1, &data->LP_lo_s2);
//     lo_sine = lo * data->sine[(int)data->cosPhase[0]];
//     lo_cosine = lo * data->sine[(int)data->sinPhase[0]];
//     lo_out = lo_sine + lo_cosine;
//     // split signal further into lo_mid
//     lo_mid =  HP_LRFilter(pre_lo, &data->coef_HP_lo_mid, &data->HP_lo_mid_s1, &data->HP_lo_mid_s2);
//     lo_mid_sine = lo_mid * data->sine[(int)data->cosPhase[1]];
//     lo_mid_cosine = lo_mid * data->sine[(int)data->sinPhase[1]];
//     lo_mid_out = lo_mid_sine + lo_mid_cosine;
//     // split signal further into hi_mid
//     hi_mid = LP_LRFilter(pre_hi, &data->coef_LP_hi_mid, &data->LP_hi_mid_s1, &data->LP_hi_mid_s2);
//     hi_mid_sine = hi_mid * data->sine[(int)data->cosPhase[2]];
//     hi_mid_cosine = hi_mid * data->sine[(int)data->sinPhase[2]];
//     hi_mid_out = hi_mid_sine + hi_mid_cosine;
//     // split signal further into hi
//     hi = HP_LRFilter(pre_hi, &data->coef_HP_hi, &data->HP_hi_s1, &data->HP_hi_s2);
//     hi_sine = hi * data->sine[(int)data->cosPhase[3]];
//     hi_cosine = hi * data->sine[(int)data->sinPhase[3]];
//     hi_out = hi_sine + hi_cosine;

        // split signal into 2 parts
        pre_lo = LP_LRFilter(input, &data->coef_LP_pre, &data->LP_pre_s1, &data->LP_pre_s2);
        pre_hi = HP_LRFilter(input, &data->coef_HP_pre, &data->HP_pre_s1, &data->HP_pre_s2);
        // split signal further into lo
        lo = LP_LRFilter(pre_lo, &data->coef_LP_lo, &data->LP_lo_s1, &data->LP_lo_s2);
        lo_sine = Hil_real(lo, &data->coef_Hil1) * data->cosine[(int)data->cosPhase[0]];
        lo_cosine = Hil_imaginary(lo, &data->coef_Hil1) * data->sine[(int)data->sinPhase[0]];
        lo_out = lo_sine + lo_cosine;
        // split signal further into lo_mid
        lo_mid = HP_LRFilter(pre_lo, &data->coef_HP_lo_mid, &data->HP_lo_mid_s1, &data->HP_lo_mid_s2);
        lo_mid_sine = Hil_real(lo_mid, &data->coef_Hil2) * data->cosine[(int)data->cosPhase[1]];
        lo_mid_cosine = Hil_imaginary(lo_mid, &data->coef_Hil2) * data->sine[(int)data->sinPhase[1]];
        lo_mid_out = lo_mid_sine + lo_mid_cosine;
        // split signal further into hi_mid
        hi_mid = LP_LRFilter(pre_hi, &data->coef_LP_hi_mid, &data->LP_hi_mid_s1, &data->LP_hi_mid_s2);
        hi_mid_sine = Hil_real(hi_mid, &data->coef_Hil3) * data->cosine[(int)data->cosPhase[2]];
        hi_mid_cosine = Hil_imaginary(hi_mid, &data->coef_Hil3) * data->sine[(int)data->sinPhase[2]];
        hi_mid_out = hi_mid_sine + hi_mid_cosine;
        // split signal further into hi
        hi = HP_LRFilter(pre_hi, &data->coef_HP_hi, &data->HP_hi_s1, &data->HP_hi_s2);
        hi_sine = Hil_real(hi, &data->coef_Hil4) * data->cosine[(int)data->cosPhase[3]];
        hi_cosine = Hil_imaginary(hi, &data->coef_Hil4) * data->sine[(int)data->sinPhase[3]];
        hi_out = hi_sine + hi_cosine;

       *out++ = lo_out + lo_mid_out + hi_mid_out + hi_out;
//      func_out = Hil_real(data->sine[(int)data->sinPhase[0]], &data->coef_Hil1);
//      *out++ = lo + lo_mid + hi_mid + hi;
//      *out++ = hi;
//      *out++ = data->sine[(int)data->sin1Phase];

        // write output to buffer for libsndfile
        *(data->bufferPointer++) = lo_out + lo_mid_out + hi_mid_out + hi_out;
    }
    return paNoError;
}


int main(void)
{
    PaStream *stream;
    PaError err;
    paData data;

    int totalSamples = (int)(SR*DUR);
    double buffer[totalSamples];

    double shift_lo = 30.0;
    double shift_lo_mid = 50.0;
    double shift_hi_mid = 70.0;
    double shift_hi = 110.0;

    double first_split_freq = 1000.0;
    double second_lo_split_freq = 400.0;
    double second_hi_split_freq = 2000.0;

    int i;
    // set LP and HP filter samples to 0
    for (i=0; i<2; i++) {
      data.LP_pre_s1.x[i] = 0;
      data.LP_pre_s1.y[i] = 0;
      data.HP_pre_s1.x[i] = 0;
      data.HP_pre_s1.y[i] = 0;

      data.LP_lo_s1.x[i] = 0;
      data.LP_lo_s1.y[i] = 0;
      data.HP_lo_mid_s1.x[i] = 0;
      data.HP_lo_mid_s1.y[i] = 0;

      data.LP_hi_mid_s1.x[i] = 0;
      data.LP_hi_mid_s1.y[i] = 0;
      data.HP_hi_s1.x[i] = 0;
      data.HP_hi_s1.y[i] = 0;
    }
    // set LP and HP filter coefficients to 0
    for (i=0; i<3; i++) {
        data.coef_LP_pre.a[i] = 0;
        data.coef_LP_pre.b[i] = 0;
        data.coef_HP_pre.a[i] = 0;
        data.coef_HP_pre.b[i] = 0;

        data.coef_LP_lo.a[i] = 0;
        data.coef_LP_lo.b[i] = 0;
        data.coef_HP_lo_mid.a[i] = 0;
        data.coef_HP_lo_mid.b[i] = 0;

        data.coef_LP_hi_mid.a[i] = 0;
        data.coef_LP_hi_mid.b[i] = 0;
        data.coef_HP_hi.a[i] = 0;
        data.coef_HP_hi.b[i] = 0;
    }
    // set Hilbert filter samples to 0;
    for (i=0; i<12; i++) {
       data.coef_Hil1.xnm1[i] = 0;
       data.coef_Hil1.ynm1[i] = 0;
       data.coef_Hil2.xnm1[i] = 0;
       data.coef_Hil2.ynm1[i] = 0;
       data.coef_Hil3.xnm1[i] = 0;
       data.coef_Hil3.ynm1[i] = 0;
       data.coef_Hil4.xnm1[i] = 0;
       data.coef_Hil4.ynm1[i] = 0;
    }
    // calculate coefficients for LP and HP filters
    LP_coef(SR, first_split_freq, &data.coef_LP_pre);
    HP_coef(SR, first_split_freq, &data.coef_HP_pre);
    LP_coef(SR, second_hi_split_freq, &data.coef_LP_hi_mid);
    HP_coef(SR, second_hi_split_freq, &data.coef_HP_hi);
    LP_coef(SR, second_lo_split_freq, &data.coef_LP_lo);
    HP_coef(SR, second_lo_split_freq, &data.coef_HP_lo_mid);

    // calculate coefficients for Hilbert filters
    Hil_coef(SR, &data.coef_Hil1);
    Hil_coef(SR, &data.coef_Hil2);
    Hil_coef(SR, &data.coef_Hil3);
    Hil_coef(SR, &data.coef_Hil4);

    // set modulator frequencies and phase
    float modFreq[4] = {shift_lo, shift_lo_mid, shift_hi_mid, shift_hi};
    for (i=0; i<4; i++) {
        data.sinPhaseIncr[i] = data.cosPhaseIncr[i] = TABLE_SIZE * modFreq[i]/SR;
        data.sinPhase[i] = data.cosPhase[i] = 0.0;
    }
    // create sine and cosine tables for modulators
    for (i=0; i<TABLE_SIZE; i++) {
        data.sine[i] = (float)sin( ((double)i/(double)TABLE_SIZE)*M_PI*2.0);
        data.cosine[i] = (float)cos( ((double)i/(double)TABLE_SIZE)*M_PI*2.0);
    }

    // libsndfile buffer
    for (i=0; i<totalSamples; i++) {
       buffer[i] = 0;
    }
    data.bufferPointer = buffer;

    SF_INFO sfinfo;
    sfinfo.channels = 1;
    sfinfo.format = SF_FORMAT_WAV | SF_FORMAT_PCM_16;
    sfinfo.samplerate = SR;

    SNDFILE *sf = sf_open("out.wav", SFM_WRITE, &sfinfo);
    if (!sf) {
        printf("Error opening file!\n");
        return -1;
    }

    err = Pa_Initialize();
    if( err != paNoError )
        printf( "PortAudio error: %s\n", Pa_GetErrorText( err ) );

    err = Pa_OpenDefaultStream( &stream,
                                1, /* no input channels */
                                1, /* stereo output */
                                paFloat32, /* 32 bit floating point output */
                                SR,
                                FRAMES_PER_BUFFER, /* frames per buffer */
                                paCallback, /* this is your callback function */
                                &data ); /*This is a pointer that will be passed to your callback*/
    if( err != paNoError )
        printf( "PortAudio error: %s\n", Pa_GetErrorText( err ) );

    err = Pa_StartStream( stream );
    if( err != paNoError )
        printf( "PortAudio error: %s\n", Pa_GetErrorText( err ) );

    Pa_Sleep(DUR*1000);

    err = Pa_StopStream( stream );
    if( err != paNoError )
        printf( "PortAudio error: %s\n", Pa_GetErrorText( err ) );

    err= Pa_CloseStream(stream);
    if( err != paNoError )
        printf( "PortAudio error: %s\n", Pa_GetErrorText( err ) );

    err = Pa_Terminate();
    if( err != paNoError )
        printf( "PortAudio error: %s\n", Pa_GetErrorText( err ) );
    printf("\nHello World!\n");

    // write buffer to file
    sf_writef_double(sf, buffer, totalSamples);

    return 0;
}
