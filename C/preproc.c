#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <math.h>
#include "preproc.h"
//#include "pcm_data.h"

int frame_length = 200;
int buffer_length = 16000;

int32_t energy[80] = {0};
double zeroCrossingRate[80] = {0};
int result[10] = {0};

uint8_t sgn(int16_t data)
{
    return (data > 0) ? 1 : 0;
}

void calEnergy(int16_t *wave_data, uint32_t len, int32_t *energy)
{
    int i = 0, j = 0;
    int32_t sum = 0;

    for (i = 0; i < len; i++)
    {
        sum = sum + (wave_data[i] * wave_data[i]);
        if ((i + 1) % frame_length == 0)
        {
            energy[j++] = sum;
            sum = 0;
        }
        else if (i == len - 1)
        {
            energy[j++] = sum;
        }
    }
}

void calZeroCrossingRate(int16_t *wave_data, uint16_t len, double *zeroCrossingRate)
{
    int i = 0, j = 0;
    double sum = 0;
    for (i = 0; i < len; i++)
    {
        if ((i % frame_length) == 0)
        {
            continue;
        }
        sum = sum + fabs(sgn(wave_data[i]) - sgn(wave_data[i - 1]));
        if (((i + 1) % frame_length) == 0)
        {
            zeroCrossingRate[j++] = (double)sum / (frame_length - 1);
            sum = 0;
        }
        else if (i == len - 1)
        {
            zeroCrossingRate[j++] = (double)sum / (frame_length - 1);
        }
    }
}

int endPointDetect(short *wave_data, int len, int32_t *energy, double *zeroCrossingRate, int *result)
{
    double sumEnergy = 0;
    double energyAverge = 0;
    double sumCzr = 0;
    int i = 0, j = 0;
    double threholdL;
    double threholdH;
    double zs = 0;
    int lenA = 0, lenB = 0;
    int edge = 0, left = 0;

    for (i = 0; i < frameCnt; i++)
    {
        sumEnergy += energy[i];
    }
    energyAverge = sumEnergy / frameCnt;

    sumEnergy = 0;
    for (i = 0; i < 5; i++)
    {
        sumEnergy = sumEnergy + energy[i];
    }
    threholdL = sumEnergy / 5.0;
    threholdH = energyAverge / 4.0;
    threholdL = (threholdL + threholdH) / 4.0;
    sumEnergy = 0;

    for (i = 0; i < 5; i++)
    {
        sumCzr += zeroCrossingRate[i];
    }
    zs = sumCzr / 5.0;

    int *arrayA = (int *)malloc(frameCnt * sizeof(int));
    int *arrayB = (int *)malloc(frameCnt * sizeof(int));

    memset(arrayA, 0, frameCnt);
    memset(arrayB, 0, frameCnt);
    arrayA[0] = -1;

    // higher threshold
    int flag = 0;
    for (i = 0, j = 0; i < frameCnt; i++)
    {
        if (arrayA[0] == -1 && flag == 0 && energy[i] > threholdH)
        {
            arrayA[j++] = i;
            flag = 1;
        }
        else if (flag == 0 && energy[i] > threholdH && (i - 21) > arrayA[j - 1])
        {
            arrayA[j++] = i;
            flag = 1;
        }
        else if (flag == 0 && energy[i] > threholdH && (i - 21) < arrayA[j - 1])
        {
            arrayA[j - 1] == 0;
            j--;
            flag = 1;
        }

        if (flag == 1 && energy[i] < threholdH)
        {
            arrayA[j++] = i;
            flag = 0;
        }
        if (flag == 1 && energy[i] > threholdH && i == frameCnt - 1)
        {
            arrayA[j++] = i;
            flag = 0;
        }
    }
    lenA = j;

    // lower threshold
    for (i = 0, j = 0; i < lenA; i++)
    {
        edge = arrayA[i];
        if ((i % 2) == 1)
        {
            while (edge < frameCnt && energy[edge] > threholdL)
            {
                
                edge++;
            }
            arrayB[j++] = edge;
        }
        else
        {
            while (edge > 0 && energy[edge] > threholdL)
            {
                edge--;
            }
            arrayB[j++] = edge;
        }
    }
    lenB = j;
    
    for (i = 0, j = 0; i < lenB; i++)
    {
        edge = arrayB[i];
        if ((i % 2) == 1)
        {
            while (edge < frameCnt && zeroCrossingRate[edge] >= 3 * zs)
            {
                edge++;
            }
            result[j++] = edge;
        }
        else
        {
            while (edge > 0 && zeroCrossingRate[edge] >= 3 * zs)
            {
                edge--;
            }
            result[j++] = edge;
        }
    }

    free(arrayA);
    free(arrayB);
    
    return j;
}

//通过result中记录的语音活动帧，将音频数据中的语音活动部分剪辑
void cutAudio(const uint8_t* wav_input_data, uint8_t** keyframe_data, int* result, int result_len, int frame_length, int audio_len) 
{
    //(*keyframe_data) = (uint8_t*)malloc(audio_len * 2 * sizeof(uint8_t));
    int cnt = 0;
    for (int j = 0; j < result_len; j += 2) {
        int start_index = result[j] * frame_length * 2;
        int end_index = (result[j + 1]) * frame_length * 2;

        printf("[%d, %d]\n", start_index, end_index);
        for (int i = start_index; i < end_index; i++) {
            (*keyframe_data)[cnt++] = wav_input_data[i];
        }
    }
}

//封装
int preproc(const uint8_t* wav_data, uint8_t* cut_wav_data)
{
    int i, j;
    short *buffer16 = (short *)malloc(buffer_length * sizeof(short));
    FILE* fp = fopen("wav01.txt", "w+");
    for (i = 0; i < buffer_length; i++)
    {
        buffer16[i] = (short)((wav_data[2 * i + 1] << 8) | wav_data[2 * i]);
        fprintf(fp, "%d\n", buffer16[i]);
    }
    fclose(fp);

    calEnergy(buffer16, buffer_length, energy);
    
    calZeroCrossingRate(buffer16, buffer_length, zeroCrossingRate);

    int result_len = endPointDetect(buffer16, buffer_length, energy, zeroCrossingRate, result);
    int audio_len = 0;
    for (int i = 0; i < result_len; i += 2) {
        audio_len += (result[i+1] - result[i]) * frame_length;
    }

    printf("audio_len: %d\n", audio_len);
    cutAudio(wav_data, &cut_wav_data, result, result_len,frame_length, audio_len);
    
    FILE* fp1 = fopen("wav01_cut.txt", "w+");
    for (i = 0; i < audio_len; i++)
    {
        buffer16[i] = (short)((cut_wav_data[2 * i + 1] << 8) | cut_wav_data[2 * i]);
        fprintf(fp1, "%d\n", buffer16[i]);
    }
    fclose(fp1);

    free(buffer16);
    return audio_len * 2;
}


