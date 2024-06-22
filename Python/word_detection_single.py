import numpy as np
import librosa
import os
import pyaudio
import wave
import tkinter as tk
from tkinter import messagebox

CHUNK = 1024
FORMAT = pyaudio.paInt16
CHANNELS = 2  # 单声道
RATE = 48000  # 采样率设置为48000
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

    n_mels= 40  # 增加梅尔滤波器的数量
    n_mfcc= 20  # 增加MFCC特征的数量

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
    
    # 检查形状是否一致，如果不一致则进行截断或补零处理
    if M1.shape[1] != M2.shape[1]:
        min_len = min(M1.shape[1], M2.shape[1])
        M1 = M1[:, :min_len]
        M2 = M2[:, :min_len]

    cost = np.zeros((M1_len, M2_len))
    dis = np.zeros((M1_len, M2_len))

    for i in range(M1_len):
        for j in range(M2_len):
            dis[i, j] = np.linalg.norm(M1[i] - M2[j])  # 使用欧几里得距离

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

    audio_dir = "wav_test_single"  # 文件夹名称

    WAVE_OUTPUT_FILENAME = os.path.join(audio_dir, "output.wav")
    audio = pyaudio.PyAudio()
    stream = audio.open(format=FORMAT, channels=CHANNELS,
                rate=RATE, input=True,
                frames_per_buffer=CHUNK)
    print("/*************开始录音*************/")
    frames = []
    for i in range(0, int(RATE / CHUNK * RECORD_SECONDS)):
        data = stream.read(CHUNK)
        frames.append(data)
    print("/*************结束录音*************/")

    stream.stop_stream()
    stream.close()
    audio.terminate()
    waveFile = wave.open(WAVE_OUTPUT_FILENAME, 'wb')
    waveFile.setnchannels(CHANNELS)
    waveFile.setsampwidth(audio.get_sample_size(FORMAT))
    waveFile.setframerate(RATE)
    waveFile.writeframes(b''.join(frames))
    waveFile.close()

    wav_single_test_path = "wav_test_single"
    for file_name in os.listdir(wav_single_test_path):
        if file_name.endswith(".wav"):
            disArray = np.full(10, np.inf)

            file_path = os.path.join(wav_single_test_path, file_name)
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
            print(f"wav:{var_name};result={result}")

            if result == 6:# Hi 芯原
                root = tk.Tk()
                root.withdraw()
                print("\r\nHi 芯原\r\n")
                messagebox.showinfo("识别结果", "Hi 芯原",)
                root.destroy()
            elif result == 7:# 测体温
                root = tk.Tk()
                root.withdraw()
                print("\r\n测体温\r\n")
                messagebox.showinfo("识别结果", "测体温")
                root.destroy()
            elif result == 8:# 测血压
                root = tk.Tk()
                root.withdraw()
                print("\r\n测血压\r\n")
                messagebox.showinfo("识别结果", "测血压")  
                root.destroy()              
            elif result == 9:# 测血糖
                root = tk.Tk()
                root.withdraw()
                print("\r\n测血糖\r\n")
                messagebox.showinfo("识别结果", "测血糖")
                root.destroy()
            else:
                root = tk.Tk()
                root.withdraw()
                print("\r\n无效命令\r\n")
                messagebox.showinfo("识别结果", "无效命令")
                root.destroy()




