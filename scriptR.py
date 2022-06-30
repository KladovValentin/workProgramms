
import sys
import os
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

# Я оставляю так, в любом случае беру только те, у которых хотя бы в одной модели вероятности быть каоном/пионом обе > 0.5 (или домноженное на хи 2 отношение), так как 
# (вообще нужно долю таких, где и то и то маленькое посчитать для событий и фона) плохие по вероятности модели провалят отбор дальше по крайней мере. Можно разве что сместить порог итоговый пониже смотря по резам анализа эфф
# В принципе можно выбор модели тоже запихать в небольшую сетку, но анализ в любом случае поможет хотя бы для описания причин и за счет чего мл сеть может работать в данном случае. + понимание. 
# Дальше конечно нет уверенности, что хи2 меньше у той модели где лучше вероятности для событий фона (если нет, то хорошо разделять должно), но при этом есть уверенность, что хи2 больше у неправильной модели (по вероятностям) ккпи.
# Если вероятности и там и там > 0.5 (редко, <~1%), то надо брать их и выбирать чтобы одна из моделей была хороша и по хи и по вер, так как для ккпи оба должны быть в согласии, а для фонов не факт, 
# Брать только те в которых > 0.7 условно в ту или другую сторону после домножения на отношение хи 2 для каждой из вероятностей. Это в принципе можно применить ко всем, не только тем у которых и то и то > 0.5.
# В конце сделать маленькую нейронку на датабазе ккпи, 2к2пи (мб) и 4пи, на параметрах 2х вероятностей, хи2 ккпи и хи2 4пи. Оно будет норм потому что у ккпи будет обычно большое произведение вероятностей на хи2к/хи2пи, 
# А для допустим 4пи будут в основном меньше и вероятности и отношение хи квадратов.


# train
trainYN = 0
file1 = open("trainSet.txt","r")
X = []
Y = []
Z = []
for line in file1:
  features = line.strip().split('	')
  momenta = float(features[0])
  i = int(momenta*10)//1
  i = i-int((i)//9 * ((i+1)%9))
  X.append([float(features[0]),float(features[1]),float(features[2])])
  Y.append(int(features[3]))
  Z.append(i) 
XTrain = []
YTrain = []
ZTrain = []
XTest = []
YTest = []
ZTest = []
for i in range(len(Y)):
  randind = random.random()
  if randind < 0.95:
    XTrain.append(X[i])
    YTrain.append(Y[i])
    ZTrain.append(Z[i])
  elif randind >= 0.95:
    XTest.append(X[i])
    YTest.append(Y[i])
    ZTest.append(Z[i])
XTrain = np.array(XTrain)
YTrain = np.array(YTrain)
ZTrain = np.array(ZTrain)
XTest = np.array(XTest)
YTest = np.array(YTest)
ZTest = np.array(ZTest)
clft = []
for x in range(9):
  clft.append(tree.DecisionTreeClassifier(max_depth=7))
for k in range(9):
#  boolArr = (ZTrain == k+1)
  boolArr = (ZTrain == k)
  XV = XTrain[boolArr]
  YV = YTrain[boolArr]
  ZV = ZTrain[boolArr]
  print(XV)
  clft[k] = clft[k].fit(XV, YV)
file1.close()
XV = XTrain
YV = YTrain
for k in range(9):
  if trainYN == 1:
    dump(clft[k], 'clf' + str(k) + '.joblib') 


#____________________________________________
clf = []
for k in range(9):
  clf.append(load('clf' + str(k) + '.joblib'))
file1C = open("trainSet.txt","r")
acc=0
all=0
#for line in file1C:
for j in range(len(YTest)):
  #features = line.strip().split('	')
  #x=[float(features[0]),float(features[1]),float(features[2])]
  #y=int(features[3])
  #momenta = float(features[0])
  #i = int(momenta*10)//1
  #i = i-int((i)//9 * ((i+1)%9))
  i = ZTest[j]
  #if y == 0 or y == 1:
  if YTest[j] == 0 or YTest[j] == 1:
    all+=1.0
  #if clf[i].predict([x])[0] == y:
  #if clf[i-1].predict([XTest[j]])[0] == YTest[j]:
  if clf[i].predict([XTest[j]])[0] == YTest[j]:
    acc+=1.0
acc = acc/all
print(clf[2].predict_proba([[0.2,3.5,0.0]])[0])
print(clf[2].predict_proba([[0.5,0.5,0.2]])[0])
print(acc)
file1C.close()


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
Indexes = []
Processes = []
Chisquares = []
bad = 0
bad2pi = 0
bad2ka = 0
badnone = 0
badhz = 0
good1 = 0
good2 = 0
goodBoth = 0
all = 0
percK1 = 0
percP1 = 0
percK2 = 0
percP2 = 0
nonenone = 0
nonekk = 0
nonepp = 0
for line in file11:
  features = line.strip().split('	')
  xp1=[float(features[0]),float(features[1]),float(features[2])]
  xk1=[float(features[3]),float(features[4]),float(features[5])]
  x21=float(features[6])
  xp2=[float(features[7]),float(features[8]),float(features[9])]
  xk2=[float(features[10]),float(features[11]),float(features[12])]
  x22=float(features[13])
  x24pi=float(features[14])
  process=int(features[15]) #0-background,1-kkp,2-exp
  Processes.append(process)
  all += 1
  momentap1 = float(features[0])
  ip1 = int(momentap1*10)//1
  momentak1 = float(features[3])
  ik1 = int(momentak1*10)//1
  momentap2 = float(features[7])
  ip2 = int(momentap2*10)//1
  momentak2 = float(features[10])
  ik2 = int(momentak2*10)//1
  ip1 = ip1-int((ip1)//9 * ((ip1+1)%9))
  ik1 = ik1-int((ik1)//9 * ((ik1+1)%9))
  ip2 = ip2-int((ip2)//9 * ((ip2+1)%9))
  ik2 = ik2-int((ik2)//9 * ((ik2+1)%9))
  percP1 = 1
  percK1 = 1
  percP2 = 1
  percK2 = 1
  if ip1 < 9:
    percP1 = clf[ip1].predict_proba([xp1])[0][1]
  else:
    percP1 = 0
    print("bad")
  if ik1 < 9:
    percK1 = clf[ik1].predict_proba([xk1])[0][0]
  else:
    percK1 = 1
    print("bad")
  if ip2 < 9:
    percP2 = clf[ip2].predict_proba([xp2])[0][1]
  else:
    precP2 = 0
    print("bad")
  if ik2 < 9:
    percK2 = clf[ik2].predict_proba([xk2])[0][0]
  else:
    percK2 = 1
    print("bad")
  
  if (percK1 >= 0.5 and percP1 >= 0.5) and (percK2 < 0.5 or percP2 < 0.5):
    good1+=1
    XE.append([float(percP1),float(percK1),float(x21),float(x24pi)])
    file23.write(str(percP1) + " " + str(percK1) + " " + str(x21) + " " + str(x24pi))
    Indexes.append(int(1))
    Chisquares.append(x21)
    #file2.write("1  1 " + str(x21) + "\n")
    if percK2 < 0.5 and percP2 < 0.5 and x22 <90:
      nonenone+=1
    elif percK2 > 0.5 and x22 < 90:
      nonekk+=1
    elif percP2 > 0.5 and x22 < 90:
      nonepp+=1
  elif (percK2 >= 0.5 and percP2 >= 0.5) and (percK1 < 0.5 or percP1 < 0.5):
    good2+=1
    XE.append([float(percP2),float(percK2),float(x22),float(x24pi)])
    file23.write(str(percP2) + " " + str(percK2) + " " + str(x22) + " " + str(x24pi))
    Indexes.append(int(2))
    Chisquares.append(x22)
    #file2.write("1  2 " + str(x22) + "\n")
  elif percK2 >= 0.5 and percP2 >= 0.5 and percK1 >= 0.5 and percP1 >= 0.5:
    goodBoth+=1
    #print(str(percP1*percK1/(percP2*percK2)*x22/x21) + " " + str(x21) + " " + str(x22))
    if(percP1*percK1/(percP2*percK2)*x22/x21>=1.5):
      XE.append([float(percP1),float(percK1),float(x21),float(x24pi)])
      file23.write(str(percP1) + " " + str(percK1) + " " + str(x21) + " " + str(x24pi))
      Indexes.append(int(1))
      Chisquares.append(x21)
      #file2.write("1  1 " + str(x21) + "\n")
    elif(percP2*percK2/(percP1*percK1)*x21/x22>=1.5):
      XE.append([float(percP2),float(percK2),float(x22),float(x24pi)])
      file23.write(str(percP2) + " " + str(percK2) + " " + str(x22) + " " + str(x24pi))
      Indexes.append(int(2))
      Chisquares.append(x22)
      #file2.write("1  2 " + str(x22) + "\n")
    else:
      #print("not good")
      XE.append([float(percP1),float(percK1),float(x21),float(x24pi)])
      file23.write(str(percP1) + " " + str(percK1) + " " + str(x21) + " " + str(x24pi))
      Indexes.append(int(0))
      Chisquares.append(100)
      #file2.write("0  0 100.0\n")
  else:
    bad+=1
    if (percK1 < 0.5 and percP1 >= 0.5) and (percK2 < 0.5 and percP2 > 0.5):
      bad2pi+=1
      Indexes.append(int(4))
      Chisquares.append(x21)
    elif (percK1 < 0.5 and percP1 > 0.5) and (percK2 < 0.5 and percP2 >= 0.5):
      bad2pi+=1
      Indexes.append(int(5))
      Chisquares.append(x22)
    elif (percK1 >= 0.5 and percP1 < 0.5) or (percK2 >= 0.5 and percP2 < 0.5):
      bad2ka+=1
      Indexes.append(int(0))
      Chisquares.append(100)
    elif (percK1 < 0.5 and percP1 < 0.5) and (percK2 < 0.5 and percP2 < 0.5):
      badnone+=1
      Indexes.append(int(3))
      Chisquares.append(min(x22,x21))
    else:
      badhz+=1
      Indexes.append(int(0))
      Chisquares.append(100)
    XE.append([float(percP1),float(percK1),float(x21),float(x24pi)])
    file23.write(str(percP1) + " " + str(percK1) + " " + str(x21) + " " + str(x24pi))
    #file2.write("0  0 100.0\n")
  
  if mod == "train":
     file23.write(" " + str(process) + " " + str(int(Indexes[-1])))
  file23.write("\n")
file11.close()
#file2.close()
file23.close()
print(str(all) + " " + str(good1) + " " +  str(good2) + " " +  str(goodBoth) + " " + str(bad)+ " " + str(bad2pi)+ " " + str(bad2ka)+ " " + str(badnone)+ " " + str(badhz))
print(str(nonenone) + " " + str(nonekk) + " " + str(nonepp))
XE=np.array(XE) 
YE=np.array(YE)
Indexes=np.array(Indexes)

# test 
file2 = open("testMarks" + mod + ".txt","w")
allkkp = 0
selectedkkpT = 0
selectedkkpF = 0
all4pi = 0
selected4piT = 0
selected4piF = 0
all = 0
for i in range(Indexes.size):
  if Processes[i] == 1:
    allkkp += 1
  if Processes[i] == 0:
    all4pi += 1
  #if clfFinal.predict_proba([XE[i]])[0][1]>0.5 and Indexes[i]>0.5: #y = 1 for kkpi, indexes>0.5 to determine what kinfit out of 2 to use
  #if (XE[i][3] > XE[i][2]-15 or XE[i][3]>15) and Indexes[i]>0.5:
  if Indexes[i]==1 or Indexes[i]==2 or Indexes[i]==4 or Indexes[i]==5:
    file2.write("1  " + str(Indexes[i]) + " " + str(Chisquares[i]) + "\n")
    if Processes[i] == 1:
      selectedkkpT += 1
    if Processes[i] == 0:
      selected4piF += 1
  else:
    file2.write("0  " + str(Indexes[i]) + " " + str(Chisquares[i]) + "\n")
    if Processes[i] == 1:
      selectedkkpF += 1
    if Processes[i] == 0:
      selected4piT += 1
print("\n")
if mod == "train":
  print("process selection efficiency: " + str(selectedkkpT/allkkp))
  print("process selection inefficiency: " + str(selectedkkpF/allkkp))
  print("background selection efficiency: " + str(selected4piF/all4pi))
  print("background supression: " + str(selected4piT/all4pi))
file2.close()

# x axis values
x = [1,2,3]
# corresponding y axis values
y = [2,4,1]
 
# plotting the points
#plt.plot(x, y)
#plt.show()

#boolArr = (ZTrain == 4)
#XV = XTrain[boolArr]
#YV = YTrain[boolArr]
#ZV = ZTrain[boolArr]
#viz = dtreeviz(clf[4],XV,YV, target_name="target", feature_names=["p","dEx","amp"], class_names=["kaon","pion"])
#viz.view()