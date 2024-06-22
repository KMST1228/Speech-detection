import numpy as np
import librosa
import os
import pyaudio
import wave

CHUNK = 1024
FORMAT = pyaudio.paInt16
CHANNELS = 2
RATE = 48000
RECORD_SECONDS = 3

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

def dtw(M1, M2):
    M1_len = len(M1)
    M2_len = len(M2)

    if M1.shape[1] != M2.shape[1]:
        min_len = min(M1.shape[1], M2.shape[1])
        M1 = M1[:, :min_len]
        M2 = M2[:, :min_len]

    cost = np.zeros((M1_len, M2_len))
    dis = np.zeros((M1_len, M2_len))

    for i in range(M1_len):
        for j in range(M2_len):
            dis[i, j] = np.linalg.norm(M1[i] - M2[j])

    cost[0, 0] = dis[0, 0]
    for i in range(1, M1_len):
        cost[i, 0] = cost[i - 1, 0] + dis[i, 0]
    for j in range(1, M2_len):
        cost[0, j] = cost[0, j - 1] + dis[0, j]

    for i in range(1, M1_len):
        for j in range(1, M2_len):
            cost[i, j] = dis[i, j] + min(cost[i - 1, j], cost[i, j - 1], cost[i - 1, j - 1])

    return cost[M1_len - 1, M2_len - 1]

def distance(x1, x2):
    return np.linalg.norm(x1 - x2)

if __name__ == "__main__":
    totalCnt = 0
    rightCnt = 0
    falseCnt = 0

    wav_test_path = "wav_test"
    mfcc_path = "mfcc"
    files = os.listdir(mfcc_path)

    for file_name in os.listdir(wav_test_path):
        if file_name.endswith(".wav"):
            disArray = np.full(10, np.inf)

            file_path = os.path.join(wav_test_path, file_name)
            mfcc = extract_MFCC(file_path)
            var_name = file_name.split(".")[0]
            for file in files:
                mfcc_file = os.path.join(mfcc_path, file)
                mfcc1 = np.loadtxt(mfcc_file)
                dis = dtw(mfcc, mfcc1)

                cnt = int(file.split("-")[1].split(".")[0])
                disArray[cnt] = dis

            result = np.argmin(disArray)
            totalCnt += 1
            if result == int(var_name.split("-")[1]):
                rightCnt += 1
                flag = "Yes"
            else:
                falseCnt += 1
                flag = "x"
            print(f"wav:{var_name};result={result} {flag}")
    correctRate = rightCnt / totalCnt
    print(f"CORRECT_RATE: {correctRate:.4f}")
