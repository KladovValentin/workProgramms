
import sys
import os
import subprocess
os.environ["PATH"] += os.pathsep + 'D:/Program Files/Graphviz/bin/'
from sklearn import tree
import numpy as np
import random
from matplotlib import pyplot as plt
import pickle
from dtreeviz.trees import dtreeviz
from joblib import dump, load
import torch
from torch import nn
from torch.utils.data import DataLoader
from torchvision import datasets, transforms

#device = "cuda:0" if torch.cuda.is_available() else "cpu"
#print(torch.cuda.get_device_name(0))
#print(f"Using {device} device")
#cuda = torch.device('cuda:0')

#pcutds = subprocess.run("python pcutds.py")
#ptrain = subprocess.run("python trainML.py NN")
ptest = subprocess.run("python predict.py NN")




"""
# test  
mod = "mod" #for train or test (mod or exp)

namei = "i"
nameo = "o"
if mod == "train":
  data = data2 = ""
  with open("testSetMod.txt") as fp:
    data = fp.read()
  with open("testSetBack.txt") as fp:
    data2 = fp.read()
  data += data2
  with open ("trainSet2.txt", 'w') as fp:
    fp.write(data)
  namei = "trainSet2.txt"
  nameo = "train2.txt"
elif mod == "mod":
  namei = "testSetMod.txt"
  nameo = "testMarksMod.txt"
else:
  namei = "testSetExp.txt"
  nameo = "testMarksExp.txt"

file11 = open(namei,"r")
file23 = open(nameo,"w")
XE = []
YE = []

# test 
file2 = open("testMarks" + mod + ".txt","w")

file2.close()

# x axis values
x = [1,2,3]
# corresponding y axis values
y = [2,4,1]
 """
# plotting the points
#plt.plot(x, y)
#plt.show()

#boolArr = (ZTrain == 4)
#XV = XTrain[boolArr]
#YV = YTrain[boolArr]
#ZV = ZTrain[boolArr]
#viz = dtreeviz(clf[4],XV,YV, target_name="target", feature_names=["p","dEx","amp"], class_names=["kaon","pion"])
#viz.view()