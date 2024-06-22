#include "mfcc.h"
#include "preproc.h"
#include "pcm_data.h"

int main()
{
    //从pcm_data中提取的语音数据长度为pcm_length*
    uint8_t cut_wav_data1[16000] = {0};
    uint8_t cut_wav_data2[16000] = {0};
    int pcm_length1 = preproc(pcm_data_cut_1, cut_wav_data1);
    int pcm_length2 = preproc(pcm_data_cut_2, cut_wav_data2);

    uint8_t *wavedata1 = (uint8_t *)malloc(pcm_length1 * sizeof(uint8_t));
    uint8_t *wavedata2 = (uint8_t *)malloc(pcm_length2 * sizeof(uint8_t));

    for (int i = 0; i < pcm_length1 / 2; i++)
    {
        wavedata1[i] = cut_wav_data1[i];
    }
    for (int i = 0; i < pcm_length2 / 2; i++)
    {
        wavedata2[i] = cut_wav_data2[i];
    }

    printf("len1:%d\nlen2:%d\n", pcm_length1, pcm_length2);

    int buffer_length1 = pcm_length1 / 2;
    int buffer_length2 = pcm_length2 / 2;
    int number_feature_vectors1 = (buffer_length1 - frame_length) / hop_length + 1;
    int number_feature_vectors2 = (buffer_length2 - frame_length) / hop_length + 1;

    double **mfcc1 = (double **)malloc(number_feature_vectors1 * sizeof(double *));
    for (int i = 0; i < number_feature_vectors1; i++)
    {
        mfcc1[i] = (double *)malloc(3 * n_mfcc * sizeof(double));
        memset(mfcc1[i], 0, 3 * n_mfcc * sizeof(double));
    }
    double **mfcc2 = (double **)malloc(number_feature_vectors2 * sizeof(double *));
    for (int i = 0; i < number_feature_vectors2; i++)
    {
        mfcc2[i] = (double *)malloc(3 * n_mfcc * sizeof(double));
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

    // free(cut_wav_data1);
    // free(cut_wav_data2);
    return 0;
}
