import sys

import torch
import torch.nn as nn
from torch.utils.data import Dataset, SubsetRandomSampler
from matplotlib import cm
import pickle
import random
import matplotlib.pyplot as plt

import numpy as np
import io
#import trainML

x_coord = 0
y_coord = 1
z_coord = 0
initial_static = 0
size_x = 100
size_y = 100


def get_x_test_and_y_weights(mod):
    f = io.open("testSet" + mod + ".txt", 'r')
    x_test = []
    local_x = []
    for line in f:
        local_x = []
        tokens = line.strip().split('	')
        for j in range(8):
            local_x.append(float(tokens[j]))
        local_x.append(float(tokens[8])/20.0)
        local_x.append(float(tokens[8])/20.0-float(tokens[9])/20.0)
        local_x.append(float(tokens[10]))
        if local_x != [] and float(tokens[10])>=50./750 and float(tokens[10])<1825./750.:
            x_test.append(local_x)
    f.close()

    #print(x_test)

    return list(np.array(x_test))


def predict_nn(mod):
    x_test = get_x_test_and_y_weights(mod)

    nn_model = nn.Sequential(
        nn.Linear(x_test[0].shape[0], 256, bias=True),
        nn.ReLU(inplace=True),
        nn.Linear(256, 1024),
        #nn.ReLU(inplace=True),
        #nn.Linear(1024, 1024),
        nn.ReLU(inplace=True),
        nn.Linear(1024, 4),
    )
    nn_model.type(torch.FloatTensor)
    nn_model.load_state_dict(torch.load("nnMode750.pt"))
    nn_model1 = nn.Sequential(
        nn.Linear(x_test[0].shape[0], 256, bias=True),
        nn.ReLU(inplace=True),
        nn.Linear(256, 1024),
        #nn.ReLU(inplace=True),
        #nn.Linear(1024, 1024),
        nn.ReLU(inplace=True),
        nn.Linear(1024, 4),
    )
    nn_model1.type(torch.FloatTensor)
    nn_model1.load_state_dict(torch.load("nnMode825.pt"))
    nn_model2 = nn.Sequential(
        nn.Linear(x_test[0].shape[0], 256, bias=True),
        nn.ReLU(inplace=True),
        nn.Linear(256, 1024),
        #nn.ReLU(inplace=True),
        #nn.Linear(1024, 1024),
        nn.ReLU(inplace=True),
        nn.Linear(1024, 4),
    )
    nn_model2.type(torch.FloatTensor)
    nn_model2.load_state_dict(torch.load("nnMode900.pt"))
    nn_model3 = nn.Sequential(
        nn.Linear(x_test[0].shape[0], 256, bias=True),
        nn.ReLU(inplace=True),
        nn.Linear(256, 1024),
        #nn.ReLU(inplace=True),
        #nn.Linear(1024, 1024),
        nn.ReLU(inplace=True),
        nn.Linear(1024, 4),
    )
    nn_model3.type(torch.FloatTensor)
    nn_model3.load_state_dict(torch.load("nnMode1100.pt"))

    #cross check
    x_test_cc = []
    for j in range(int(len(x_test)/2)):
        if (random.random()>0):
            x_test_cc.append(x_test[2*j])
            x_test_cc.append(x_test[2*j+1])
    x_test_cc = np.array(x_test_cc)

    #predict answer
    #pred = nn_model(torch.FloatTensor(x_test_cc)).softmax(dim=1)

    #y will be a vector of lines with the lenght of the number of classes
    y = []
    #print(pred[10, 0].item())
    counttempj=0
    for j in range(len(x_test_cc)):
        pred = 0
        counttempj+=1
        if x_test_cc[j,10]>=150./750 and x_test_cc[j,10]<750./750.:
            pred = nn_model(torch.FloatTensor(x_test_cc[j])).softmax(dim=0)
            print(str(counttempj) + " / " + str(len(x_test_cc)))
        elif x_test_cc[j,10]>=750./750 and x_test_cc[j,10]<825./750.:
            pred = nn_model1(torch.FloatTensor(x_test_cc[j])).softmax(dim=0)
            print(str(counttempj) + " / " + str(len(x_test_cc)))
        elif x_test_cc[j,10]>=825./750 and x_test_cc[j,10]<900./750.:
            pred = nn_model2(torch.FloatTensor(x_test_cc[j])).softmax(dim=0)
            print(str(counttempj) + " / " + str(len(x_test_cc)))
        elif x_test_cc[j,10]>=900./750 and x_test_cc[j,10]<1100./750.:
            pred = nn_model3(torch.FloatTensor(x_test_cc[j])).softmax(dim=0)
            print(str(counttempj) + " / " + str(len(x_test_cc)))
        tempy = [pred[0].item(), pred[1].item(), pred[2].item(), pred[3].item()]
        #tempy = [pred[j, 0].item(), pred[j, 1].item(), pred[j, 2].item(), pred[j, 3].item()]
        y.append(tempy)
    print(y[0])
    y = list(np.array(y))
    info = ""
    f = open("predictOutput.txt", "w")

    #for i in range(len(pred[0])):
    #    info += str(y[i]) + "\n"

    #calculate_efficiencies(y,x_test_cc)
    write_output(y,x_test_cc,mod)
    draw_probabilities_spread(y)
    #f.write(info)
    f.close()


def draw_NN():
    x_test = get_x_test_and_y_weights()

    nn_model = nn.Sequential(
        nn.Linear(x_test[0].shape[0], 256, bias=True),
        nn.ReLU(inplace=True),
        nn.Linear(256, 1024),
        nn.ReLU(inplace=True),
        nn.Linear(1024, 4),
    )
    nn_model.type(torch.FloatTensor)
    nn_model.load_state_dict(torch.load("nnModel.pt"))
    draw_train_results("NN", nn_model)


def calculate_efficiencies(outputs, testdata):
    efficiency = 0
    for i in range(len(outputs)/2):
        max = 0
        predictedClass = 0
        for j in range(4):
            if outputs[2*i][j]>max:
                max = outputs[2*i][j]
                predictedClass = j
        max2 = 0
        predictedClass2 = 0
        for j in range(4):
            if outputs[2*i+1][j]>max2:
                max2 = outputs[2*i+1][j]
                predictedClass2 = j
        
        if predictedClass == 0 and (predictedClass2 !=0 or testdata[2*i][8]):
            efficiency += 1
    print("efficiency = " + str(100*float(efficiency)/float(len(outputs))))


def draw_probabilities_spread(outputs):
    class_hist = []
    class_hist_selected = []
    for i in range(len(outputs)):
        max = 0
        predictedClass = 0
        current_probas = []
        for j in range(4):
            current_probas.append(outputs[i][j])
            if outputs[i][j]>max:
                max = outputs[i][j]
                predictedClass = j
        class_hist.append(current_probas)
        if predictedClass == 0:
            class_hist_selected.append(current_probas)
    class_hist = np.array(class_hist)
    class_hist_selected = np.array(class_hist_selected)

    bins = np.linspace(0, 1, 20)

    plt.rc('axes', titlesize=22)
    plt.rc('axes', labelsize=22)
    plt.rc('xtick', labelsize=20)
    plt.rc('legend',fontsize=22)

    plt.hist(class_hist_selected[:,2], bins, color='#0504aa',
                            alpha=0.7, rwidth=0.95, label = '$\pi-\pi$')
    plt.hist(class_hist_selected[:,0], bins, color='#228B22',
                            alpha=0.7, rwidth=0.95, label = '$K-\pi$')
    plt.legend(loc=[0.6,0.8])
    plt.grid(axis='y', alpha=0.75)
    plt.xlabel('probability')
    plt.show()


def write_output(outputs, testdata,mod):
    efficiency = 0
    allel = 0
    bad100 = 0
    info = ""
    f = open("testMarks" + mod + ".txt", "w")
    for i in range(int(len(outputs)/2)):
        max1 = 0
        allel+=1
        predictedClass1 = 0
        for j in range(4):
            if outputs[2*i][j]>max1:
                max1 = outputs[2*i][j]
                predictedClass1 = j
        max2 = 0
        predictedClass2 = 0
        for j in range(4):
            if outputs[2*i+1][j]>max2:
                max2 = outputs[2*i+1][j]
                predictedClass2 = j
        if (predictedClass1 == 0 and predictedClass2 != 0 and outputs[2*i][0]-outputs[2*i][1]>0.05 and outputs[2*i][0]-outputs[2*i][2]>0.07 and outputs[2*i][0]-outputs[2*i][3]>0.05 and outputs[2*i][0]>0.08):
            info += "1  " + "1  " + str(testdata[2*i][8]*20) + "\n"
            efficiency+=1
        elif (predictedClass2 == 0 and predictedClass1 != 0 and outputs[2*i+1][0]-outputs[2*i+1][1]>0.05 and outputs[2*i+1][0]-outputs[2*i+1][2]>0.07 and outputs[2*i+1][0]-outputs[2*i+1][3]>0.05 and outputs[2*i+1][0]>0.08):
            info += "1  " + "2  " + str(testdata[2*i+1][8]*20) + "\n"
            efficiency+=1
        elif (predictedClass1 == 0 and predictedClass2 == 0 and outputs[2*i][0]-outputs[2*i][1]>0.05 and outputs[2*i][0]-outputs[2*i][2]>0.07 and outputs[2*i][0]-outputs[2*i][3]>0.05 and outputs[2*i][0]>0.08):
            info += "1  " + "1  " + str(testdata[2*i][8]*20) + "\n"
            efficiency+=1
        elif (predictedClass1 != 0 and predictedClass2 != 0):
            if (testdata[2*i][8]*20 == 100 or testdata[2*i+1][8]*20 == 100):
                bad100+=1
                info += "0  " + "3  " + str(testdata[2*i][8]*20) + "\n"
            else:
                info += "0  " + "0  " + str(testdata[2*i][8]*20) + "\n"
        else:
            info += "0  " + "0  " + str(testdata[2*i][8]*20) + "\n"
    print("efficiency = " + str(100*float(efficiency)/float(allel-bad100)))
    print("100 somewhere 0 0 = " + str(100*float(bad100)/float(allel)))
    f.write(info)
    f.close()


def draw_train_results(model_name, model, aux_model=None):
    train_dataset = trainML.DESY_dataset("informationCorrected.txt")
    X = []
    Y = []
    Z = []
    if np.array(train_dataset[:][0][0, :]).shape[0] > 1:
        init_params = np.array(train_dataset[initial_static][0])
        for i in range(len(train_dataset)):
            count_as_dataset = True
            for j, val in enumerate(init_params):
                if j != x_coord and j != y_coord:
                    if val != train_dataset[i][0][j]:
                        count_as_dataset = False
            if count_as_dataset:
                X.append(train_dataset[i][0][x_coord].item())
                Y.append(train_dataset[i][0][y_coord].item())
                Z.append(train_dataset[i][1][z_coord].item())
        X = np.array(X)
        Y = np.array(Y)
        Z = np.array(Z)
        x_surf, y_surf = np.meshgrid(
            np.linspace(X.min() - (X.max() - X.min()) / 10, X.max() + (X.max() - X.min()) / 10, size_x),
            np.linspace(Y.min() - (Y.max() - Y.min()) / 10, Y.max() + (Y.max() - Y.min()) / 10, size_y))

        list_to_model = []
        for i in range(np.array(train_dataset[:][0][0, :]).shape[0]):
            if i == x_coord:
                list_to_model.append(x_surf.flatten())
            elif i == y_coord:
                list_to_model.append(y_surf.flatten())
            else:
                list_to_model.append(train_dataset[initial_static][0][i].item() * np.ones(size_x * size_y))
        model_viz = np.expand_dims(np.array(list_to_model).T, 0)
        if model_name == "NN":
            z_surf = np.array(model(torch.FloatTensor(model_viz))[0, :, z_coord].tolist())
        elif aux_model is not None:
            z_surf = (np.array(model(torch.FloatTensor(model_viz))[0, :, z_coord].tolist()) + aux_model.predict(
                model_viz)[:, z_coord]) / 2

        fig = plt.figure(figsize=(14, 8))
        ax = fig.add_subplot(111, projection='3d')
        ax.scatter(X, Y, Z, c='red', marker='o', alpha=0.5)
        ax.plot_surface(x_surf, y_surf, z_surf.reshape(x_surf.shape), alpha=0.4, cmap=cm.coolwarm, linewidth=0,
                        antialiased=False)
        _, _, names, namesX = trainML.load_dataset("informationCorrected.txt")
        ax.set_xlabel(namesX[0][x_coord])
        ax.set_ylabel(namesX[0][y_coord])
        ax.set_zlabel(names[0][z_coord])
        plt.show()


def predict(model):
    if model == "NN":
        predict_nn("Mod")
    if model == "NN draw":
        predict_nn()


print("start python predict")
#if sys.argv[2]:
#    x_coord = int(sys.argv[2])
#if sys.argv[3]:
#    y_coord = int(sys.argv[3])
#if sys.argv[4]:
#    z_coord = int(sys.argv[4])
#if sys.argv[5]:
#    initial_static = int(sys.argv[5])
predict("NN")
