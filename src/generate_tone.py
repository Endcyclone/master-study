# coding=utf-8

import pyaudio
import numpy as np
import threading
from time import time, sleep

def init_global():
    global TT, phase, freq, newfreq, level, level_flag
    TT = time()
    freq = 440
    newfreq = 440
    phase = 0
    level = 0
    level_flag = False
    return

def set_freq(freq):
    global newfreq
    newfreq = freq

def set_level(newlevel):
    global level, level_flag, level_filter
    ft_len = 1024
    r = np.arange(0, ft_len) * np.pi / ft_len
    diff = level - newlevel
    level_filter = (0.5 + np.sign(diff) * 0.5 * np.cos(r)) ** 0.5 * abs(diff) + min(level, newlevel)
    level_flag = True
    level = newlevel
    return

def do_level_filter(left):
    global level_flag, level_filter, level
    if level_flag:
        left *= level_filter
        level_flag = False
    else:
        left *= level
    return left

def callback(in_data, frame_count, time_info, status):
    global TT, phase, freq, newfreq
    if newfreq != freq:
        phase = 2*np.pi*TT*(freq-newfreq)+phase
        freq = newfreq
    left = do_level_filter(np.sin(phase+2*np.pi*freq*(TT+np.arange(frame_count)/44100)))
    data = np.zeros((left.shape[0]*2,), np.float32)
    data[::2] = left
    data[1::2] = left
    TT += frame_count/44100
    return (data, pyaudio.paContinue)

class generate_tone(threading.Thread):
    def __init__(self):
        threading.Thread.__init__(self)
        #self.setDaemon(True)
        init_global()
        self.p = pyaudio.PyAudio()
        self.stream = self.p.open(format=pyaudio.paFloat32, channels=2, rate=44100, output=True, stream_callback=callback)
        self.running = True
        return
    def run(self):
        self.stream.start_stream()
        while self.running:
            sleep(.1)
        self.stream.stop_stream()
        self.stream.close()
        self.p.terminate()