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
import pandas

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


def writeEnPointDataset(datF,index):
    lastName = list(datF.columns)[-1]
    x = [len(datF[(datF[lastName]==i)]) for i in range(4)]
    minimumCount = np.amin(np.array(x))
    frames = [datF.loc[datF[lastName] == i].sample(minimumCount) for i in range(4)]
    resultTable = pandas.concat(frames).reset_index(drop=True)
    if (len(resultTable) > 0):
        resultTable.to_csv("trainSetCorr" + str(index) + ".csv")
        return index+1
    return index

def cutDataset(fname):
    df = pandas.read_table(fname,sep='	',header=None)
    dfs = dict(tuple(df.groupby(list(df.columns)[10])))
    listdf = [dfs[x] for x in dfs]
    inn = 0
    for dft in listdf:
        inn = writeEnPointDataset(dft,inn)



cutDataset("trainSet.txt")

