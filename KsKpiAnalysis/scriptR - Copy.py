
import sys
import os
os.environ["PATH"] += os.pathsep + 'D:/Program Files/Graphviz/bin/'
from sklearn import tree
import numpy as np
from matplotlib import pyplot as plt
import pickle
from dtreeviz.trees import dtreeviz

# Я оставляю так, в любом случае беру только те, у которых хотя бы в одной модели вероятности быть каоном/пионом обе > 0.5 (или домноженное на хи 2 отношение), так как 
# (вообще нужно долю таких, где и то и то маленькое посчитать для событий и фона) плохие по вероятности модели провалят отбор дальше по крайней мере. Можно разве что сместить порог итоговый пониже смотря по резам анализа эфф
# В принципе можно выбор модели тоже запихать в небольшую сетку, но анализ в любом случае поможет хотя бы для описания причин и за счет чего мл сеть может работать в данном случае. + понимание. 
# Дальше конечно нет уверенности, что хи2 меньше у той модели где лучше вероятности для событий фона (если нет, то хорошо разделять должно), но при этом есть уверенность, что хи2 больше у неправильной модели (по вероятностям) ккпи.
# Если вероятности и там и там > 0.5 (редко, <~1%), то надо брать их и выбирать чтобы одна из моделей была хороша и по хи и по вер, так как для ккпи оба должны быть в согласии, а для фонов не факт, 
# Брать только те в которых > 0.7 условно в ту или другую сторону после домножения на отношение хи 2 для каждой из вероятностей. Это в принципе можно применить ко всем, не только тем у которых и то и то > 0.5.
# В конце сделать маленькую нейронку на датабазе ккпи, 2к2пи (мб) и 4пи, на параметрах 2х вероятностей, хи2 ккпи и хи2 4пи. Оно будет норм потому что у ккпи будет обычно большое произведение вероятностей на хи2к/хи2пи, 
# А для допустим 4пи будут в основном меньше и вероятности и отношение хи квадратов.


# train
file1 = open("trainSet.txt","r")
X = []
Y = []
for line in file1:
  features = line.strip().split('	')
  X.append([float(features[0]),float(features[1]),float(features[2])])
  Y.append(int(features[3]))
X=np.array(X)
Y=np.array(Y)
clf = tree.DecisionTreeClassifier(max_depth=7)
clf = clf.fit(X, Y)
file1.close()
XV = X
YV = Y


file1C = open("trainSet.txt","r")
acc=0
all=0
for line in file1C:
  features = line.strip().split('	')
  x=[float(features[0]),float(features[1]),float(features[2])]
  y=int(features[3])
  if y == 0 or y == 1:
    all+=1.0
  if clf.predict([x])[0] == y:
    acc+=1.0
acc = acc/all
print(acc)
file1C.close()


# test 
mod = "train" #for train or test (mod or exp)

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
good1 = 0
good2 = 0
goodBoth = 0
all = 0
percK1 = 0
percP1 = 0
percK2 = 0
percP2 = 0
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
  percP1 = clf.predict_proba([xp1])[0][1]
  percK1 = clf.predict_proba([xk1])[0][0]
  percP2 = clf.predict_proba([xp2])[0][1]
  percK2 = clf.predict_proba([xk2])[0][0]
  if percK1 >= 0.5 and percP1 >= 0.5 and (percK2 < 0.5 or percP2 < 0.5):
    good1+=1
    XE.append([float(percP1),float(percK1),float(x21),float(x24pi)])
    file23.write(str(percP1) + " " + str(percK1) + " " + str(x21) + " " + str(x24pi))
    Indexes.append(int(1))
    Chisquares.append(x21)
    #file2.write("1  1 " + str(x21) + "\n")
  elif percK2 >= 0.5 and percP2 >= 0.5 and (percK1 < 0.5 or percP1 < 0.5):
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
    XE.append([float(percP1),float(percK1),float(x21),float(x24pi)])
    file23.write(str(percP1) + " " + str(percK1) + " " + str(x21) + " " + str(x24pi))
    Indexes.append(int(0))
    Chisquares.append(100)
    #file2.write("0  0 100.0\n")
  
  if mod == "train":
     file23.write(" " + str(process) + " " + str(int(Indexes[-1])))
  file23.write("\n")
file11.close()
#file2.close()
file23.close()
print(str(all) + " " + str(good1) + " " +  str(good2) + " " +  str(goodBoth) + " " + str(bad))
XE=np.array(XE) 
YE=np.array(YE)
Indexes=np.array(Indexes)


#train
file1 = open("train2.txt","r")
X = []
Y = []
for line in file1:
  #print(line)
  features = line.strip().split(' ')
  #print(features[4])
  #print(features[5])
  if (int(float(features[4])) == 0 or int(float(features[4]) == 1)) and int(float(features[5])) > 0.5:
    X.append([float(features[0]),float(features[1]),float(features[2]),float(features[3])])
    Y.append([int(features[4])])
X=np.array(X)
Y=np.array(Y)
clfFinal = tree.DecisionTreeClassifier(max_depth=10)
clfFinal = clfFinal.fit(X, Y)
file1.close()

file1C = open("train2.txt","r")
acc=0
all=0
allkkp = 0
all4pi = 0
selectedkkp = 0
selected4pi = 0
for line in file1C:
  features = line.strip().split(' ')
  if (int(float(features[4])) == 0 or int(float(features[4]) == 1)) and int(float(features[5])) > 0.5:
    x=[float(features[0]),float(features[1]),float(features[2]),float(features[3])]
    y=int(float(features[4]))
    all+=1.0
    if y == 1:
      allkkp+=1
    if y == 0:
      all4pi+=1
    if clfFinal.predict([x])[0] == y:
      acc+=1.0
      if y == 1:
        selectedkkp+=1
      else:
        selected4pi+=1
acc = acc/all
print("2nd ml accuracy: " + str(acc))
print("2nd ml process selection efficiency: " + str(selectedkkp/allkkp))
print("2nd ml 4pi selection efficiency: " + str(selected4pi/all4pi))
file1C.close()


# test 
file2 = open("testMarksMod.txt","w")
allkkp = 0
selectedkkpT = 0
selectedkkpF = 0
all4pi = 0
selected4piT = 0
selected4piF = 0
all = 0
for i in range(Indexes.size-1):
  if Processes[i] == 1:
    allkkp += 1
  if Processes[i] == 0:
    all4pi += 1
  #if clfFinal.predict_proba([XE[i]])[0][1]>0.5 and Indexes[i]>0.5: #y = 1 for kkpi, indexes>0.5 to determine what kinfit out of 2 to use
  if (XE[i][3] > XE[i][2]-15 or XE[i][3]>15) and Indexes[i]>0.5:
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
plt.plot(x, y)
plt.show()

viz = dtreeviz(clf,XV,YV, target_name="target", feature_names=["p","dEx","amp"], class_names=["kaon","pion"])
viz.view()