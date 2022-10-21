
import torch
import torch.nn as nn
import matplotlib.pyplot as plt
import pandas
import numpy as np


def constructPredictionFrame(index,nparray,chi2Array):
    # return new DataFrame with 4 columns corresponding to probabilities of resulting classes and 5th column to chi2 
    df2 = pandas.DataFrame(index=index)
    df2['0'] = nparray[:,0].tolist()
    df2['1'] = nparray[:,1].tolist()
    df2['2'] = nparray[:,2].tolist()
    df2['3'] = nparray[:,3].tolist()
    df2['chi2'] = chi2Array.tolist()
    return df2


def loadModel(energyForm, nInputParameters):
    #load nn from file
    name = 'networks\\model_'+str("{:.3f}".format(energyForm))+'.pt'
    nn_model = nn.Sequential(
        nn.Linear(nInputParameters, 256, bias=True),
        nn.BatchNorm1d(256),
        nn.LeakyReLU(inplace=True),
        nn.Linear(256, 1024),
        nn.BatchNorm1d(1024),
        nn.LeakyReLU(inplace=True),
        nn.Linear(1024, 4),
    )
    nn_model.type(torch.FloatTensor)
    nn_model.load_state_dict(torch.load(name))
    return nn_model


def makePredicionList(mod,meanValues,stdValues,trainEnergies):
    df = pandas.read_table("testSet" + mod + ".txt",sep='	',header=None)
    dfs = dict(tuple(df.groupby(list(df.columns)[10])))
    listdf = [dfs[x] for x in dfs]
    dat_list = []

    for dft in listdf:
        columnName = list(dft.columns)[10]
        lastName = list(dft.columns)[11]
        energyForm = round(dft.reset_index(drop=True).loc[0].at[columnName]*750,3)

        #remove energy column, transform to numpy and assign type
        dfn = dft.drop([columnName,lastName], axis=1).to_numpy().astype(np.float32)
        chi2Array = dfn[:,8].copy()

        if not energyForm in trainEnergies:
            #make prediction probabilities = 0 for example idk
            dat_list.append(constructPredictionFrame(dft.index.copy(),np.zeros(dfn.shape),chi2Array))
            continue
        
        #apply mean and std to this part of dataset
        energyIndex = np.where(trainEnergies == energyForm)[0][0]
        for j in range(dfn.shape[1]):
            dfn[:,j] = (dfn[:,j]-meanValues[energyIndex][j])/stdValues[energyIndex][j]

        #load nn
        nn_model = loadModel(energyForm, dfn[0].shape[0])
            
        #predict for each row
        resultTArr = nn_model(torch.FloatTensor(dfn)).softmax(dim=1).detach().numpy()
        dat_list.append(constructPredictionFrame(dft.index.copy(),resultTArr,chi2Array))

    fullPredictionList = pandas.concat(list(dat_list)).sort_index()

    return fullPredictionList


def draw_probabilities_spread(outputs,selected):
    class_hist = outputs.to_numpy()
    class_hist_selected = selected.to_numpy()

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


def write_output(outputs, mod):

    # split output table in 3 - good, bad by nn and bad by input (no energy point in train or 100 chi2)
    badChi2Table = outputs.loc[(outputs['chi2']==100) | ((outputs['0']==0) & (outputs['1']==0))].copy()
    possGTable = outputs.drop(badChi2Table.index)
    goodTable = possGTable.loc[(possGTable['0']-possGTable['1']>0.05) & (possGTable['0']-possGTable['2']>0.07) & (possGTable['0']-possGTable['3']>0.05) & (possGTable['0']>0.1)].copy()
    badEventTable = possGTable.drop(goodTable.index)

    draw_probabilities_spread(outputs,goodTable)

    # assing final marks to each input row
    goodTable['0'] = 1
    badEventTable['0'] = 0
    badChi2Table['0'] = 3

    resultingMarks = pandas.concat([goodTable,badEventTable,badChi2Table]).sort_index().drop(['1','2','3'], axis=1)
    resultingMarks.to_csv("testMarks" + mod + ".csv",sep='	',header=False,index=False)


def predict_nn(mod):
    energies = np.loadtxt('energies.txt')
    meanValues = np.loadtxt('meanValues.txt').reshape((len(energies),10))
    stdValues = np.loadtxt('stdValues.txt').reshape((len(energies),10))

    predictionList = makePredicionList(mod,meanValues,stdValues,energies)
    print(predictionList)

    write_output(predictionList,mod)


def predict(model):
    if model == "NN":
        predict_nn("Mod")
    if model == "NN draw":
        predict_nn()


print("start python predict")
predict("NN")
