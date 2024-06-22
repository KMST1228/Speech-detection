#ifndef __PREPROC_H__
#define __PREPROC_H__

#define frameCnt 80

uint8_t sgn(int16_t data);
void calEnergy(int16_t *wave_data, uint32_t len, int32_t *energy);
void calZeroCrossingRate(int16_t *wave_data, uint16_t len, double *zeroCrossingRate);
int endPointDetect(short *wave_data, int len, int32_t *energy, double *zeroCrossingRate, int *result);
void cutAudio(const uint8_t* wav_input_data, uint8_t** keyframe_data, int* result, int result_len, int frame_length, int audio_len);
int preproc(const uint8_t* wav_data, uint8_t* cut_wav_data);

#endif
