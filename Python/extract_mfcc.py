import numpy as np
import librosa
import os
import re

def extract_MFCC(wav_file):
    sample_rate = 8000
    frame_t = 25
    hop_length_t = 10
    frame_length = 200
    hop_length = 80

    win_length = int(frame_t*sample_rate/1000)
    hop_length = int(hop_length_t*sample_rate/1000)
    n_fft = int(2**np.ceil(np.log2(win_length)))

    n_mels= 40
    n_mfcc= 20

    data, fs = librosa.load(wav_file, sr=sample_rate)
    data_trim, index = librosa.effects.trim(data, top_db=20, frame_length=frame_length, hop_length=hop_length)

    mfcc = librosa.feature.mfcc(y=data_trim, 
                                sr=sample_rate, 
                                n_mfcc=n_mfcc,
                                n_mels=n_mels,
                                n_fft=n_fft, 
                                win_length=win_length, 
                                hop_length=hop_length)
    mfcc = np.delete(mfcc, obj=0, axis=0)

    mfcc_delta = librosa.feature.delta(mfcc)
    mfcc_delta2 = librosa.feature.delta(mfcc, order=2)
    mfcc_d1_d2 = np.concatenate([mfcc, mfcc_delta, mfcc_delta2], axis=0)

    mfcc_d1_d2 = np.transpose(mfcc_d1_d2)
    return mfcc_d1_d2

if __name__ == "__main__":
    wav_path = "wav"
    for i in range(1, 7):
        file_path = os.path.join(wav_path, str(i))
        dirs = os.listdir(file_path)
        for file in dirs:
            if re.match(r'^0.*\.wav$', file):
                wav_file = os.path.join(file_path, file)
                mfcc = extract_MFCC(wav_file)
                name = f'./mfcc/{file.split(".")[0]}.txt'
                np.savetxt(name, mfcc, fmt="%s")
