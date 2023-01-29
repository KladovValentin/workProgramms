
import sys
import torch
import torch.nn as nn
import matplotlib.pyplot as plt
import pandas
import numpy as np


class RNN(torch.nn.Module):
    def __init__(self, input_dim, embedding_dim, hidden_dim, output_dim):
        super().__init__()
        self.hidden_dim=hidden_dim
        self.input_dim = input_dim
        self.embedding_dim = embedding_dim
        self.embedding = nn.Embedding(input_dim, embedding_dim)
        self.rnn = nn.RNN(embedding_dim, hidden_dim, nonlinearity='tanh', batch_first=True)
        #self.rnn = torch.nn.LSTM(embedding_dim, hidden_dim)        
        
        self.fc = torch.nn.Linear(hidden_dim, output_dim)
        self.nn_model = nn.Sequential(
            nn.Linear(input_dim, 256, bias=True),
            nn.BatchNorm1d(256),
            #nn.Dropout(0.5),
            nn.LeakyReLU(inplace=True),
            nn.Linear(256, 1024),
            nn.BatchNorm1d(1024),
            #nn.Dropout(0.5),
            nn.LeakyReLU(inplace=True),
            nn.Linear(1024, 1024),
            nn.BatchNorm1d(1024),
            nn.LeakyReLU(inplace=True),
            nn.Linear(1024, 256),
            nn.BatchNorm1d(256),
            #nn.Dropout(0.5),
            nn.LeakyReLU(inplace=True),
            nn.Linear(256, embedding_dim),
            nn.BatchNorm1d(embedding_dim)
        )

    def forward(self, text):
        # text dim: [batch size, sentence length]
        batch_size = text.size(0)
        # Initializing hidden state for first input using method defined below
        hidden = self.init_hidden(batch_size)
        #print(text.reshape((batch_size*2, self.input_dim)).shape)
        embedded = self.nn_model(text.reshape((batch_size*2, self.input_dim))).reshape((batch_size, 2, self.embedding_dim))
        # embedded dim: [batch size, sentence length=2, embedding dim=4]
        #print(text.shape)
        output, hidden = self.rnn(embedded, hidden)
        #print(hidden.shape)
        # output dim: [batch size, sentence length=2, hidden dim]
        # hidden dim: [1, batch size, hidden dim]

        hidden.squeeze_(0)
        # hidden dim: [batch size, hidden dim]
        
        output1 = self.fc(hidden)
        #print(output1.shape)
        #print("next")
        return output1

    def init_hidden(self, batch_size):
        # This method generates the first hidden state of zeros which we'll use in the forward pass
        # We'll send the tensor holding the hidden state to the device we specified earlier as well
        hidden = torch.zeros(1, batch_size, self.hidden_dim)
        return hidden


def loadModel(energyForm, nInputParameters):
    #load nn from file
    name = 'networks\\model_'+str("{:.3f}".format(energyForm))+'.pt'
    #print(name)
    """
    nn_model = nn.Sequential(
        nn.Linear(nInputParameters, 256, bias=True),
        nn.BatchNorm1d(256),
        nn.Dropout(0.5),
        nn.LeakyReLU(inplace=True),
        nn.Linear(256, 1024),
        nn.BatchNorm1d(1024),
        nn.Dropout(0.5),
        nn.LeakyReLU(inplace=True),
        nn.Linear(1024, 4096),
        nn.BatchNorm1d(4096),
        nn.Dropout(0.5),
        nn.LeakyReLU(inplace=True),
        nn.Linear(4096, 256),
        nn.BatchNorm1d(256),
        nn.Dropout(0.5),
        nn.LeakyReLU(inplace=True),
        nn.Linear(256, 4),
    )
    """
    nn_model = RNN(input_dim=4,
               embedding_dim=256,
                  hidden_dim=512,
                  output_dim=4 
    ) # could use 1 for binary classification
    nn_model.type(torch.FloatTensor)
    nn_model.load_state_dict(torch.load(name))
    return nn_model


def constructPredictionFrame(index,nparray,chi2Array):
    # return new DataFrame with 4 columns corresponding to probabilities of resulting classes and 5th column to chi2 
    df2 = pandas.DataFrame(index=index)
    df2['0'] = nparray[:,0].tolist()
    df2['1'] = nparray[:,1].tolist()
    df2['2'] = nparray[:,2].tolist()
    df2['3'] = nparray[:,3].tolist()
    df2['chi2'] = chi2Array.tolist()
    return df2


def checkDistributions():
    dfm = pandas.read_table("testSet" + "Mod" + ".txt",sep='	',header=None)
    dfe = pandas.read_table("testSet" + "Exp" + ".txt",sep='	',header=None)
    numberOfColumns = len(dfm.columns)
    columnNames = list(dfe.columns)
    dfe[(dfe[columnNames[10]]==1000.0) & (dfe[columnNames[0]]<0.35)][columnNames[2]].plot(kind='kde')
    dfm[(dfm[columnNames[10]]==1000.0) & (dfm[columnNames[0]]<0.35)][columnNames[2]].plot(kind='kde')
    plt.show()


def makePredicionList(mod,meanValues,stdValues,trainEnergies):
    #df = pandas.read_table("testSet" + mod + ".txt",sep='	',header=None)
    dftrain = pandas.read_table("trainSet.txt", sep='	',header=None)
    trainIndexName = list(dftrain.columns)[11]
    df = dftrain.loc[((dftrain[trainIndexName]!=4))]
    dfs = dict(tuple(df.groupby(list(df.columns)[10])))
    listdf = [dfs[x] for x in dfs]
    dat_list = []
    en_list = []

    for dft in listdf:
        x2kfColName = list(dft.columns)[8]
        x24pColName = list(dft.columns)[9]
        columnName = list(dft.columns)[10]
        #energyForm = round(dft.reset_index(drop=True).loc[0].at[columnName],3)
        
        marksArr = (dft.to_numpy().astype(int))[:,11].copy()
        dftCopy = dft.copy()
        dft.drop([trainIndexName], axis=1, inplace=True)
        
        dftWithoutIndex = dft.copy().reset_index(drop=True)
        energyForm = round(dftWithoutIndex.loc[0].at[columnName],3)
        
        #remove energy column, get array of chi2, transform to numpy and assign type
        dfnChi = dftWithoutIndex.drop([columnName], axis=1).to_numpy().astype(np.float32)
        chi2Array = dfnChi[:,8].copy()
        dfn = dftWithoutIndex.drop([x2kfColName,x24pColName,columnName], axis=1).to_numpy().astype(np.float32)

        if not energyForm in trainEnergies:
            print("no such energy " + str(energyForm))
            dat_list.append(constructPredictionFrame(dft.index.copy(),np.zeros((dfn.shape[0],4)),chi2Array))
            en_list.append(energyForm)
            continue
        #print(dfn)
        #normalize dataset 
        energyIndex = np.where(trainEnergies == energyForm)[0][0]
        cond2 = (dfn[:,2]==-1)
        cond6 = (dfn[:,6]==-1)
        for j in range(dfn.shape[1]):
            dfn[:,j] = (dfn[:,j]-meanValues[energyIndex][j])/stdValues[energyIndex][j]
        dfn[cond2,2] = 0
        dfn[cond2,3] = 0
        dfn[cond6,6] = 0
        dfn[cond6,7] = 0

        #load nn and predict
        nn_model = loadModel(energyForm, dfn[0].shape[0])
        nn_model.eval()
        resultTArr = nn_model(torch.tensor(dfn.reshape(dfn.shape[0],2,4))).softmax(dim=1).detach().numpy()
        dat_list.append(constructPredictionFrame(dft.index.copy(),resultTArr,chi2Array))
        en_list.append(energyForm)

        mask = (np.argmax(resultTArr,axis=1) != marksArr)
        print(resultTArr[mask])
        print(dftCopy[mask])
            
    
    plt.show()

    #fullPredictionList = pandas.concat(list(dat_list)).sort_index()

    return list(dat_list),list(en_list)


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


def write_output(outputs, mod, enlist):
    goodTableList = []
    badEventTableList = []
    badChi2TableList = []
    i = 0
    for dft in outputs:
        # split output table in 3 - good, bad by nn and bad by input (no energy point in train or 100 chi2)
        badChi2Table = dft.loc[(dft['chi2']==100) | ((dft['0']==0) & (dft['1']==0))].copy()
        possGTable = dft.drop(badChi2Table.index)
        goodTable0 = possGTable.loc[(possGTable['0']-possGTable['1']>0.2) & (possGTable['0']-possGTable['2']>0.4) & (possGTable['0']-possGTable['3']>0.2) & (possGTable['0']>0.2)].copy()
        goodTable1 = possGTable.loc[(possGTable['1']-possGTable['0']>0) & (possGTable['1']-possGTable['2']>0) & (possGTable['1']-possGTable['3']>0) & (possGTable['1']>0.2)].copy()
        goodTable = pandas.concat([goodTable0,goodTable1])
        badEventTable = possGTable.drop(goodTable0.index).copy()
        efficiency = 0
        if (possGTable.shape[0] > 0):
            efficiency = float(goodTable.shape[0])/possGTable.shape[0]
        print(str(enlist[i]) + " efficiency equals to " + str(efficiency*100))
        goodTableList.append(goodTable0)
        badEventTableList.append(badEventTable)
        badChi2TableList.append(badChi2Table)
        i=i+1

    goodTableFull = pandas.concat(list(goodTableList)).sort_index()
    badEventTableFull = pandas.concat(list(badEventTableList)).sort_index()
    badChi2TableFull = pandas.concat(list(badChi2TableList)).sort_index()

    draw_probabilities_spread(pandas.concat(outputs).sort_index(),goodTableFull)

    # assing final marks to each input row
    goodTableFull['0'] = 1
    badEventTableFull['0'] = 0
    badChi2TableFull['0'] = 3   
    
    resultingMarks = pandas.concat([goodTableFull,badEventTableFull,badChi2TableFull]).sort_index().drop(['1','2','3'], axis=1)
    resultingMarks.to_csv("testMarks" + mod + ".csv",sep='	',header=False,index=False)


def predict_nn(mod):
    energies = np.loadtxt('energies.txt')
    meanValues = np.loadtxt('meanValues.txt').reshape((len(energies),8))
    stdValues = np.loadtxt('stdValues.txt').reshape((len(energies),8))

    predictionList, enlist = makePredicionList(mod,meanValues,stdValues,energies)
    print(pandas.concat(list(predictionList)).sort_index())

    write_output(predictionList,mod,enlist)


def predict(model):
    print(model)
    if (model == "Dis"):
        checkDistributions()
        return
    predict_nn(model)


print("start python predict")
predict(sys.argv[1])
