 #include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "pcm_data.h"
#include "preproc.h"
#include <stdint.h>

const int hop_length = 80;  // ����
extern const int frame_length; // ֡��
const int n_fft = 128; // ����Ҷ����
const int n_mfcc = 13; // ��ɢ�任ά�ȣ����յõ�3*number_coefficientsά����������
const int n_mels = 20; // ����������
const int sample_rate = 8000; // ����Ƶ��(ÿ��������), ��ʾÿ��ͨ���Ĳ����ٶ�

const double pi = 3.14159265358979323846;

char str1[] = "pcm_data1.txt";
char str2[] = "mfcc_data1.txt";

// DCT实现
void Discrete_Cosine_Transform(int direction, int length, double* X)
{
    double* x = (double*)malloc(length * sizeof(double));

    for (int i = 0; i < length; i++)
    {
        x[i] = X[i];
    }
    for (int k = 0; k < length; k++)
    {
        double sum = 0;

        if (direction == 1)
        {
            for (int n = 0; n < length; n++)
            {
                sum += ((k == 0) ? (sqrt(0.5)) : (1)) * x[n] * cos(pi * (n + 0.5) * k / length);
            }
        }
        else if (direction == -1)
        {
            for (int n = 0; n < length; n++)
            {
                sum += ((n == 0) ? (sqrt(0.5)) : (1)) * x[n] * cos(pi * n * (k + 0.5) / length);
            }
        }
        X[k] = sum * sqrt(2.0 / length);
    }
    free(x);
}

void DCT(int direction, int length, double* X)
{
    if (direction == 1 || direction == -1)
    {
        Discrete_Cosine_Transform(direction, length, X);
        return;
    }
}

// FFT�ľ���ʵ��
void Fast_Fourier_Transform(int direction, int length, double* Xr, double* Xi)
{
    int log_length = (int)(log((double)length) / log(2.0));
    double pi = 3.14159265358979323846;

    for (int i = 0, j = 0; i < length; i++, j = 0)
    {
        for (int k = 0; k < log_length; k++)
        {
            j = (j << 1) | (1 & (i >> k));
        }
        if (j < i)
        {
            double t;

            t = Xr[i];
            Xr[i] = Xr[j];
            Xr[j] = t;

            t = Xi[i];
            Xi[i] = Xi[j];
            Xi[j] = t;
        }
    }
    for (int i = 0; i < log_length; i++)
    {
        int L = (int)pow(2.0, i);

        for (int j = 0; j < length - 1; j += 2 * L)
        {
            for (int k = 0; k < L; k++)
            {
                double argument = direction * -pi * k / L;

                double xr = Xr[j + k + L] * cos(argument) - Xi[j + k + L] * sin(argument);
                double xi = Xr[j + k + L] * sin(argument) + Xi[j + k + L] * cos(argument);

                Xr[j + k + L] = Xr[j + k] - xr;
                Xi[j + k + L] = Xi[j + k] - xi;
                Xr[j + k] = Xr[j + k] + xr;
                Xi[j + k] = Xi[j + k] + xi;
            }
        }
    }
    if (direction == -1)
    {
        for (int k = 0; k < length; k++)
        {
            Xr[k] /= length;
            Xi[k] /= length;
        }
    }
}

void FFT(int direction, int length, double* Xr, double* Xi)
{
    int log_length = log((double)length) / log(2.0);

    if (direction != 1 && direction != -1)
    {
        return;
    }
    if (1 << log_length != length)
    {
        return;
    }
    Fast_Fourier_Transform(direction, length, Xr, Xi);
}

// ÷��Ƶ�ʷ�Χ
double Mel_Scale(int direction, double x)
{
    switch (direction)
    {
    case -1:
        return 700.0 * (exp(x / 1125.0) - 1);
    case 1:
        return 1125.0 * log(1 + x / 700.0);
    }
    return 0;
}

// ��ȡ������MFCC����(������ֵ)������FFT��ȡ����ֵ��Mel�˲���ȡ������DCT, ��󷵻�feature_vector[]һ֡����������
void MFCC(int length_frame, int length_DFT, int number_coefficients, int number_filterbanks, int sample_rate, double* frame, double* feature_vector)
{
    double max_Mels_frequency = Mel_Scale(1, sample_rate / 2);
    double min_Mels_frequency = Mel_Scale(1, 300);
    double interval = (max_Mels_frequency - min_Mels_frequency) / (number_filterbanks + 1);

    double* filterbank = (double*)malloc(number_filterbanks * sizeof(double));
    memset(filterbank, 0, (number_filterbanks) * sizeof(double));
    double* Xr = (double*)malloc(length_DFT * sizeof(double));
    double* Xi = (double*)malloc(length_DFT * sizeof(double));

    for (int i = 0; i < length_DFT; i++)
    {
        Xr[i] = (i < length_frame) ? (frame[i]) : (0);
        Xi[i] = 0;
    }

    // FFT
    FFT(1, length_DFT, Xr, Xi);

    for (int i = 0; i < length_DFT / 2 + 1; i++)
    {
        double frequency = (sample_rate / 2) * i / (length_DFT / 2);
        double Mel_frequency = Mel_Scale(1, frequency);
        double power = (Xr[i] * Xr[i] + Xi[i] * Xi[i]) / length_frame;

        // ÷���˲�
        for (int j = 0; j < number_filterbanks; j++)
        {
            double frequency_boundary[] = { min_Mels_frequency + interval * (j + 0), min_Mels_frequency + interval * (j + 1), min_Mels_frequency + interval * (j + 2) };

            if (frequency_boundary[0] <= Mel_frequency && Mel_frequency <= frequency_boundary[1])
            {
                double lower_frequency = Mel_Scale(-1, frequency_boundary[0]);
                double upper_frequency = Mel_Scale(-1, frequency_boundary[1]);

                filterbank[j] += power * (frequency - lower_frequency) / (upper_frequency - lower_frequency);
            }
            else if (frequency_boundary[1] <= Mel_frequency && Mel_frequency <= frequency_boundary[2])
            {
                double lower_frequency = Mel_Scale(-1, frequency_boundary[1]);
                double upper_frequency = Mel_Scale(-1, frequency_boundary[2]);

                filterbank[j] += power * (upper_frequency - frequency) / (upper_frequency - lower_frequency);
            }
        }
    }

    // ȡ����
    for (int i = 0; i < number_filterbanks; i++)
    {
        if (filterbank[i] <= 0) {
            filterbank[i] = 1e-10; // ��������ֵȡ����
        }
        filterbank[i] = log(filterbank[i]);
    }

    // DCT
    DCT(1, number_filterbanks, filterbank);

    // 提取MFCC特征
    for (int i = 1; i < number_coefficients; i++)
    {
        //i - 1的目的为去除第一列数据
        feature_vector[i - 1] = filterbank[i];
    }

    free(filterbank);
    free(Xr);
    free(Xi);
}

double Get_Buffer(int index, int buffer_length, short* buffer16)
{
    if (0 <= index && index < buffer_length)
    {
        return (buffer16[index] + 0.5) / 32767.5;
    }
    return 0;
}

void getMFCC(const uint8_t* pcm_data, int pcm_length, int buffer_length, int number_feature_vectors, double **feature_vector)
{
    //获取双字节幅度数据
    short* buffer16 = (short*)malloc(buffer_length * sizeof(short));
    for (int i = 0; i < buffer_length; i++)
    {
        buffer16[i] = (short)((pcm_data[2 * i + 1] << 8) | pcm_data[2 * i]);
    }

    printf("frame count: %d\n", number_feature_vectors);


    // MFCC
    for (int i = 0; i <= buffer_length - frame_length; i += hop_length)
    {
        double* frame = (double*)malloc(frame_length * sizeof(double));
        memset(frame, 0, frame_length * sizeof(double));

        // pre-emphasis预加重
        for (int j = 0; j < frame_length; j++)
        {
            if (i + j < buffer_length)
            {
                frame[j] = Get_Buffer(i + j, buffer_length, buffer16) - 0.95 * Get_Buffer(i + j - 1, buffer_length, buffer16);
            }
            else
            {
                frame[j] = 0;
            }
        }

        // windowing加汉明窗
        for (int j = 0; j < frame_length; j++)
        {
            frame[j] *= 0.54 - 0.46 * cos(2 * pi * j / (frame_length - 1));
        }

        MFCC(frame_length, n_fft, n_mfcc, n_mels, sample_rate, frame, feature_vector[i / hop_length]);

        free(frame);
    }

    // deltas
    for (int i = 0; i < number_feature_vectors; i++)
    {
        int prev = (i == 0) ? (0) : (i - 1);
        int next = (i == number_feature_vectors - 1) ? (number_feature_vectors - 1) : (i + 1);

        for (int j = 0; j < n_mfcc; j++)
        {
            feature_vector[i][n_mfcc + j] = (feature_vector[next][j] - feature_vector[prev][j]) / 2;
        }
    }

    // delta-deltas
    for (int i = 0; i < number_feature_vectors; i++)
    {
        int prev = (i == 0) ? (0) : (i - 1);
        int next = (i == number_feature_vectors - 1) ? (number_feature_vectors - 1) : (i + 1);

        for (int j = n_mfcc; j < 2 * n_mfcc; j++)
        {
            feature_vector[i][n_mfcc + j] = (feature_vector[next][j] - feature_vector[prev][j]) / 2;
        }
    }


    // FILE* file_pcm = fopen("pcm_data1.txt", "wt");
    // FILE* file_mfcc = fopen("mfcc_data1.txt", "wt");
    FILE* file_pcm = fopen(str1, "wt");
    FILE* file_mfcc = fopen(str2, "wt");
    str1[8]++;
    str2[9]++;

    // ��.wav��MFCC����д�뵽�ļ��У�ÿ֡һ�С�ÿ��39ά���ݡ�
    for (int i = 0; i < number_feature_vectors; i++)
    {
        for (int j = 0; j < 3 * n_mfcc; j++)
        {
            fprintf(file_mfcc, "%lf ", feature_vector[i][j]);
        }
        fprintf(file_mfcc, "\n");
    }
    for (int i = 0; i < buffer_length; i++) {
        fprintf(file_pcm, "%d\n", buffer16[i]);
    }

    fclose(file_pcm);
    fclose(file_mfcc);

    free(buffer16);
}

// 两个维数相等的向量之间的距离
double distance(double* x1, double* x2, int len) {
    double sum = 0;
    for (int i = 0; i < len; i++) {
        sum += fabs(x1[i] - x2[i]);
    }
    return sum;
}

// DTW 算法...
double dtw(double** M1, int M1_len, double** M2, int M2_len, int vec_len) {
    // 初始化 cost 数组
    double** cost = (double**)malloc(M1_len * sizeof(double*));
    for (int i = 0; i < M1_len; i++) {
        cost[i] = (double*)malloc(M2_len * sizeof(double));
    }

    // 初始化 dis 数组
    double** dis = (double**)malloc(M1_len * sizeof(double*));
    for (int i = 0; i < M1_len; i++) {
        dis[i] = (double*)malloc(M2_len * sizeof(double));
        for (int j = 0; j < M2_len; j++) {
            dis[i][j] = distance(M1[i], M2[j], vec_len);
        }
    }

    // 初始化 cost 的第 0 行和第 0 列
    cost[0][0] = dis[0][0];
    for (int i = 1; i < M1_len; i++) {
        cost[i][0] = cost[i - 1][0] + dis[i][0];
    }
    for (int j = 1; j < M2_len; j++) {
        cost[0][j] = cost[0][j - 1] + dis[0][j];
    }

    // 计算 cost 矩阵
    for (int i = 1; i < M1_len; i++) {
        for (int j = 1; j < M2_len; j++) {
            double min_cost = fmin(cost[i - 1][j] + dis[i][j] * 1,
                                   fmin(cost[i - 1][j - 1] + dis[i][j] * 2,
                                        cost[i][j - 1] + dis[i][j] * 1));
            cost[i][j] = min_cost;
        }
    }

    double result = cost[M1_len - 1][M2_len - 1];

    // 释放内存
    for (int i = 0; i < M1_len; i++) {
        free(cost[i]);
        free(dis[i]);
    }
    free(cost);
    free(dis);

    return result;
}

int main () {
    uint8_t cut_wav_data1[16000] = { 0 };
    uint8_t cut_wav_data2[16000] = { 0 };
    int pcm_length1 = preproc(pcm_data_cut_1, cut_wav_data1);
    int pcm_length2 = preproc(pcm_data_cut_2, cut_wav_data2);

    uint8_t* wavedata1 = (uint8_t*)malloc(pcm_length1 * sizeof(uint8_t));
    uint8_t* wavedata2 = (uint8_t*)malloc(pcm_length2 * sizeof(uint8_t));

    for (int i = 0; i < pcm_length1 / 2; i++) {
        wavedata1[i] = cut_wav_data1[i];
    }
    for (int i = 0; i < pcm_length2 / 2; i++) {
        wavedata2[i] = cut_wav_data2[i];
    }


    printf("len1:%d\nlen2:%d\n", pcm_length1, pcm_length2);

    // int pcm_length = 32000;
    int buffer_length1 = pcm_length1 / 2;
    int buffer_length2 = pcm_length2 / 2;
    int number_feature_vectors1 = (buffer_length1 - frame_length) / hop_length + 1; // ��.wav�ж���֡
    int number_feature_vectors2 = (buffer_length2 - frame_length) / hop_length + 1; // ��.wav�ж���֡

    double** mfcc1 = (double**)malloc(number_feature_vectors1 * sizeof(double*));
    for (int i = 0; i < number_feature_vectors1; i++)
    {
        mfcc1[i] = (double*)malloc(3 * n_mfcc * sizeof(double));
        memset(mfcc1[i], 0, 3 * n_mfcc * sizeof(double));
    }
    double** mfcc2 = (double**)malloc(number_feature_vectors2 * sizeof(double*));
    for (int i = 0; i < number_feature_vectors2; i++)
    {
        mfcc2[i] = (double*)malloc(3 * n_mfcc * sizeof(double));
        memset(mfcc2[i], 0, 3 * n_mfcc * sizeof(double));
    }

    getMFCC(wavedata1, pcm_length1, buffer_length1, number_feature_vectors1, mfcc1);
    getMFCC(wavedata2, pcm_length2, buffer_length2, number_feature_vectors2, mfcc2);

    double dis = dtw(mfcc1, number_feature_vectors1, mfcc2, number_feature_vectors2, 3 * n_mfcc);
    printf("distance: %lf\n", dis);

    for (int i = 0; i < number_feature_vectors2; i++)
    {
        free(mfcc2[i]);
    }
    free(mfcc2);
    for (int i = 0; i < number_feature_vectors1; i++)
    {
        free(mfcc1[i]);
    }
    free(mfcc1);
    free(wavedata1);
    free(wavedata2);

    //free(cut_wav_data1);
    //free(cut_wav_data2);
    return 0;
}

