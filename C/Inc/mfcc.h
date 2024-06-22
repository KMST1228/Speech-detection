#ifndef __MFCC_H__
#define __MFCC_H__

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <math.h>

#include "preproc.h"

extern const int hop_length;     // 帧移
extern const int frame_length; // ֡帧长
extern const int n_fft;         // fft点数
extern const int n_mfcc;         // mfcc特征值维数，最后二阶差分后得到一个3*n_mfcc行的矩阵
extern const int n_mels;         // mel滤波器个数
extern const int sample_rate;  // 采样率

void Discrete_Cosine_Transform(int direction, int length, double *X);

void DCT(int direction, int length, double *X);

void Fast_Fourier_Transform(int direction, int length, double *Xr, double *Xi);

void FFT(int direction, int length, double *Xr, double *Xi);

double Mel_Scale(int direction, double x);

void MFCC(int length_frame, int length_DFT, int number_coefficients, int number_filterbanks, int sample_rate, double *frame, double *feature_vector);

double Get_Buffer(int index, int buffer_length, short *buffer16);

void getMFCC(const uint8_t *pcm_data, int pcm_length, int buffer_length, int number_feature_vectors, double **feature_vector);

double distance(double *x1, double *x2, int len);

double dtw(double **M1, int M1_len, double **M2, int M2_len, int vec_len);

#endif
