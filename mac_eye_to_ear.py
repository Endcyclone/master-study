import numpy as np
import pandas as pd
import sys
import matplotlib
import matplotlib.animation
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import json
import pyaudio
import os
from pynput.keyboard import Key, Listener
import threading
from time import sleep

#プロット
def pltplot(data):
    #print(len(data))
    #data = data[0:1000]
    plt.figure()
    x = range(len(data))
    plt.plot(x, data)
    plt.show()

#連続サイン波を生成する
def sines(data, rate):
    freqs = []
    for adata in data:
        tmp = np.ones(1000)
        freqs.append((600+adata*10) * tmp)
    freqs = np.concatenate(freqs)
    phazes_diff = 2 * np.pi * freqs / rate
    phazes = np.cumsum(phazes_diff)
    return np.sin(phazes)

#視覚情報を聴覚情報に
class eye_to_ear:
    def __init__(self, mapname):
        f = open(mapname, 'r')
        self.mapdata = json.load(f)
        f.close()
        return
    def show_map(self):
        print(self.mapdata)
        return
    def is_collision(self, pos):
        return self.in_objects(pos) or self.is_outside(pos)
    def in_objects(self, pos):
        obj = self.mapdata['objects']
        for i in obj:
            count = 0
            if (pos[0] > obj[i]['x'][0]) & (pos[0] < obj[i]['x'][1]): count += 1
            if (pos[1] > obj[i]['y'][0]) & (pos[1] < obj[i]['y'][1]): count += 1
            if (pos[2] > obj[i]['z'][0]) & (pos[2] < obj[i]['z'][1]): count += 1
            if count == 3: return True
        return False
    def is_outside(self, pos):
        if (pos[0] > self.mapdata['mapsize']['x_max']) | (pos[0] < 0): return True
        if (pos[1] > self.mapdata['mapsize']['y_max']) | (pos[1] < 0): return True
        if (pos[2] > self.mapdata['mapsize']['z_max']) | (pos[2] < 0): return True
        return False
    def calc_height(self):
        pos = self.pos.copy()
        obj = self.mapdata['objects']
        while(self.in_objects(pos) == False):
            pos[2] -= 0.01
            if pos[2] < 0:
                return (self.pos[2])
        return self.pos[2] - pos[2]
    def set_pos(self, x, y, z):
        pos = np.array([x, y, z], dtype="float64")
        if self.is_collision(pos):
            print("error: wrong position")
            return
        self.pos = pos
        print("now_pos =", self.pos)
        return
    def get_pos(self):
        return self.pos

    def set_deg(self, deg):
        self.deg = (deg + 360) % 360
        print("now_deg =", self.deg)
    def set_ele(self, ele):
        if ele < 0 or 180 <= ele:
            print("error: elevation must be [0,180)")
            return
        self.ele = ele
        print("now_ele =", ele)
    def set_angle(self, deg, ele):
        self.set_deg(deg)
        self.set_ele(ele)
        print("deg, ele =", deg, ele)
    def get_deg(self):
        return self.deg
    def get_ele(self):
        return self.ele
    def scan(self, deg, ele_max=120, ele_min=30, ele_step=1):
        deg = deg % 360
        disdata = []
        try:
            pos = self.pos.copy()
            acc = np.array([np.sin(np.pi*deg/180), -1*np.cos(np.pi*deg/180)])
        except ValueError:
            print("error: position or angle is not defined")
            return

        for ele in range(ele_min, ele_max, ele_step):
            pos_3d = pos.copy()
            acc_3d = acc * np.sin(np.pi*ele/180)
            acc_3d = np.append(acc_3d, -1*np.cos(np.pi*ele/180)) * 0.1
            #print(acc_3d, np.linalg.norm(acc_3d))
            dis = 0
            while(self.is_collision(pos_3d) == False):
                dis += 0.1
                pos_3d += acc_3d
            #print(ele, pos_3d)
            disdata.append(dis)
        #pltplot(disdata)
        return disdata
    #90度回転させながらscan_mapに格納
    def rotateData(self, scan_map, scandata, d):
        for i, data in enumerate(scandata):
            coo_x = int(i * np.cos(np.pi*d/180))
            coo_y = int(i * np.sin(np.pi*d/180))
            if scan_map[coo_x, coo_y] == 0:
                scan_map[coo_x, coo_y] = data
            else:
                scan_map[coo_x, coo_y] = (scan_map[coo_x, coo_y] + data) / 2

    def scan90(self, deg):
        scan_map = np.zeros([90,90])
        for d in range(90):
            scandata = self.scan(deg+d, 120, 30)
            scan_map[d] = np.array(scandata)
            #self.rotateData(scan_map, scandata, d)
        return scan_map
    def scan_square(self, size=45, step=1):
        reso = int(size/step)
        scan_map = np.zeros([reso,reso])
        ele_min = self.ele - int(size / 2)
        ele_max = ele_min + size
        for d in range(reso):
            deg = ((self.deg + 360) + (d - int(reso / 2)) * step) % 360
            scandata = self.scan(deg, ele_max=ele_max, ele_min=ele_min, ele_step=step)
            scan_map[d] = np.log(np.array(scandata))
        return scan_map.T[::-1,:]

#キー入力時の動作
def on_press(key):
    print("pressed")
    try:
        if key.char == "w" : ka.move(0)
        if key.char == "a" : ka.move(-90)
        if key.char == "s" : ka.move(180)
        if key.char == "d" : ka.move(90)
        if key.char == "c" : ka.savecsv()
    except AttributeError:
        if key == Key.up: ka.lookUp()
        if key == Key.down: ka.lookDown()
        if key == Key.left: ka.lookLeft()
        if key == Key.right: ka.lookRight()
        if key == Key.esc:
            print('キー受け付けを終了しました')
            ka.drawing = False
            ka.update_event.set()
            ps.playing = False
            return False

def on_release(key):
    pass

class key_animation:
    def __init__(self, ete, step):
        self.ete = ete
        self.step = step
        self.fig, self.ax = plt.subplots(1,1)
        self.update_event = threading.Event()
        self.updatefig()
        self.drawing = True
        #self.updating = False

    def updatefig(self):
        #self.updating = True
        self.figdata = self.ete.scan_square(step=self.step)
        ps.level = self.figdata[int(self.figdata.shape[0] / 2), int(self.figdata.shape[0] / 2)]
        #print(ps.level)
        #self.updating = False

        self.update_event.set()

    def move(self, exdir):
        deg = self.ete.get_deg() + exdir
        div = np.array([int(1*np.sin(np.pi*deg/180)), int(-1*np.cos(np.pi*deg/180)), 0])
        pos = self.ete.get_pos() + div
        if not self.ete.is_collision(pos):
            self.ete.set_pos(pos[0], pos[1], pos[2])
        self.updatefig()
    def lookUp(self):
        self.ete.set_ele(self.ete.get_ele()+15)
        self.updatefig()
    def lookDown(self):
        self.ete.set_ele(self.ete.get_ele()-15)
        self.updatefig()
    def lookRight(self):
        self.ete.set_deg(self.ete.get_deg()+15)
        self.updatefig()
    def lookLeft(self):
        self.ete.set_deg(self.ete.get_deg()-15)
        self.updatefig()
    def animate(self):
        prev = 0
        plt.imshow(self.figdata, cmap="magma_r")
        while self.drawing:
            plt.pause(.1)
            self.update_event.wait()
            self.update_event.clear()
            self.ax.clear()
            self.ax.imshow(self.figdata, cmap="magma_r")
    def savecsv(self):
        df = pd.DataFrame(self.figdata)
        filename = str(self.step)+str(self.ete.get_pos())+str(self.ete.get_deg())+str(self.ete.get_ele())+".csv"
        df.to_csv(filename)
        print("saved.")
        return

# オーディオプレーヤー
class play_sound:
        def __init__(self):
            self.p = pyaudio.PyAudio()
            self.stream = self.p.open(format=pyaudio.paFloat32, channels=1, rate=44100, output=1)
            self.level = -1
            self.playing = True

        def generate_tone(self, frequency, length, rate):
            length = int(length * rate)
            factor = float(frequency) * (np.pi * 2) / rate
            #print(length, factor, length * factor / (np.pi))
            return np.sin(np.arange(length) * factor)

        def generate_tone_ver2(self, rate = 44100):
            freqs = (600 - self.level * 100) * np.ones(44100)
            #for adata in data:
            #    tmp = np.ones(1000)
            #    freqs.append((600+adata*10) * tmp)
            #freqs = np.concatenate(freqs)
            phazes_diff = 2 * np.pi * freqs / rate
            phazes = np.cumsum(phazes_diff)
            return np.sin(phazes)

        def play(self):
            while self.playing:
                chunks = []
                if self.level < 0:
                    sleep(1)
                    continue
                freq = 600 - self.level * 100
                chunks.append(self.generate_tone(freq, 0.5, 44100))
                #chunks.append(self.generate_tone_ver2())
                chunk = np.concatenate(chunks) * 0.1
                self.stream.write(chunk.astype(np.float32).tostring())
                sleep(0.5)

            self.stream.close()
            self.p.terminate()

############### メインスレッド ###############

ps = play_sound()
ete = eye_to_ear("table in box.json")
ete.set_pos(1,1,2)
ete.set_angle(90,90)
ete.calc_height()
#scandata = ete.scan_square()
ka = key_animation(ete, step=3)

def key_listener():
    try:
        with Listener(on_press=on_press, on_release=on_release) as listener:
            listener.join()
    except:
        pass

key_thread = threading.Thread(target=key_listener)
key_thread.start()

ps_thread = threading.Thread(target=ps.play)
ps_thread.start()

#ka = key_animation(ete)
ka.animate()
