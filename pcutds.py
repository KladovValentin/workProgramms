import sys
from xml.dom import minicompat

import torch
import torch.nn as nn
import torch.optim as optim
from sklearn.ensemble import GradientBoostingRegressor
from sklearn.multioutput import MultiOutputRegressor
from torch.utils.data import Dataset, SubsetRandomSampler
from sklearn.tree import DecisionTreeRegressor
import pickle
import random

import numpy as np
import io

device = "cuda:0" if torch.cuda.is_available() else "cpu"
print(torch.cuda.get_device_name(0))
print(f"Using {device} device")
torch.cuda.empty_cache()
print(torch.cuda.memory_allocated(0))
print(torch.cuda.memory_reserved(0))
cuda = torch.device('cuda:0')


x_coord = 0
y_coord = 1
z_coord = 0
initial_static = 0
size_x = 100
size_y = 100


def cutDataset(fname):
    x = []
    y = []
    local_x = []
    local_y = 0
    f = io.open(fname, 'r')
    for i, line in enumerate(f):
        local_x = []
        local_y = 0
        tokens = line.strip().split('	')
        for j in range(11):
            local_x.append(float(tokens[j]))
        local_y = int(tokens[11])
        if local_x != []:
            x.append(local_x)
            y.append(local_y)
    f.close()
    print("data arrays are filled")
    
    print("asassfa")
    x = np.array(x)
    y = np.array(y)
    
    #making datasets the same size
    f = io.open("trainSetCorr.txt", 'w')
    typecount = [sum(yt == 0 for yt in y), sum(yt == 1 for yt in y), sum(yt == 2 for yt in y), sum(yt == 3 for yt in y)]
    print(str(typecount[0]) + "  " + str(typecount[1]) + "    " + str(typecount[2]) + "    " + str(typecount[3]))
    minumumCount = 10000000
    for j in range(4):
        if typecount[j] < minumumCount:
            minumumCount = typecount[j]
    minumumCount = typecount[0]
    for i in range(len(y)):
        writeYN = True
        for j in range(4):
            if (typecount[j] > minumumCount):
                if (y[i] == j and random.random() >= float(minumumCount)/float(typecount[j])):
                #if (y[i] == j and random.random() >= float(10000)/float(typecount[j])):
                    writeYN = False
        if writeYN:
            for j in range(11):
                f.write(str(float("{:.4f}".format(x[i][j]))) + " ")
            f.write(str(y[i]) + "\n")
    f.close()
    print("i am here")


cutDataset("trainSet.txt")

