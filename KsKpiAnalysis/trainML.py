
import sys
import torch
import torch.nn as nn
import torch.optim as optim
from torch.utils.data import Dataset, SubsetRandomSampler
import pandas
import matplotlib.pyplot as plt
import numpy as np

device = "cuda:0" if torch.cuda.is_available() else "cpu"
print(torch.cuda.get_device_name(0))
print(f"Using {device} device")
torch.cuda.empty_cache()
print(torch.cuda.memory_allocated(0))
print(torch.cuda.memory_reserved(0))
cuda = torch.device('cuda:0')


def prepareTable(datF):
    # make datasets equal sizes for each class out of 4
    lastName = list(datF.columns)[-1]
    x = [len(datF[(datF[lastName]==i)]) for i in range(4)]
    minimumCount = np.amin(np.array(x))
    frames = [datF.loc[datF[lastName] == i].sample(minimumCount) for i in range(4)]
    return pandas.concat(frames).sort_index().reset_index(drop=True)
    #print(datF)
    #return datF.sort_index().reset_index(drop=True)


def meanAndStdTable(dataTable):
    # find mean and std values for each column of the dataset (used for train dataset)
    df = dataTable
    dfn = df.to_numpy()

    x = dfn[:,:-1].astype(np.float32)
    mean = np.array( [np.mean(x[:,j]) for j in range(x.shape[1])] )
    std  = np.array( [np.std( x[:,j]) for j in range(x.shape[1])] )
    mean[2] = np.mean(x[(x[:,2]!=-1),2])
    mean[3] = np.mean(x[(x[:,2]!=-1),3])
    mean[6] = np.mean(x[(x[:,6]!=-1),6])
    mean[7] = np.mean(x[(x[:,6]!=-1),7])
    
    std[2] = np.std(x[(x[:,2]!=-1),2])
    std[3] = np.std(x[(x[:,2]!=-1),3])
    std[6] = np.std(x[(x[:,6]!=-1),6])
    std[7] = np.std(x[(x[:,6]!=-1),7])
    return mean, std


def load_dataset(dataTable, meanValues, stdValues):
    # transform to numpy, assign types, make dataset has mean = 0 and std = 1
    df = dataTable
    dfn = df.to_numpy()

    x = dfn[:,:-1].astype(np.float32)
    y = dfn[:, -1].astype(int)
    cond2 = (x[:,2]==-1)
    cond6 = (x[:,6]==-1)
    for j in range(x.shape[1]):
        x[:,j] = (x[:,j]-meanValues[j])/stdValues[j]
    x[cond2,2] = 0
    x[cond2,3] = 0
    x[cond6,6] = 0
    x[cond6,7] = 0
    print('x shape ach equals -1 ' + str(x[x[:,2]==0].shape))
    print('x mean ach equals -1= ' + str(np.mean(x[:,2])))
    x = x.reshape(x.shape[0],2,4)
    print('x shape = ' + str(x.shape))
    print('y shape = ' + str(y.shape))
    return x,  y


class DESY_dataset(Dataset):
    def __init__(self, dataTable, meanValues, stdValues, transform=None):
        self.transform = transform
        self.datasetX, self.datasetY = load_dataset(dataTable, meanValues, stdValues)

    def __len__(self):
        return len(self.datasetY)

    def __getitem__(self, index):
        if self.transform != None:
            self.datasetX[index] = self.transform(self.datasetX[index])
        return torch.tensor(self.datasetX[index]), torch.tensor(self.datasetY[index])


def train_model(model, train_loader, loss, optimizer, num_epochs,validationX,validationY, scheduler=None):
    print("start model nn train")
    loss_history = []
    train_history = []
    validLoss_history = []
    for epoch in range(num_epochs):

        loss_accum = 0
        correct_samples = 0
        total_samples = 0
        #print(len(train_loader))
        for i_step, (x, y) in enumerate(train_loader):
            x = x.to(device)
            y = y.type(torch.LongTensor)
            y = y.to(device)

            prediction = model(x)
            loss_value = loss(prediction, y)
            optimizer.zero_grad()
            loss_value.backward()
            optimizer.step()
            _, indices = torch.max(prediction, 1)
            correct_samples += torch.sum(indices == y)
            total_samples += y.shape[0]
            
            loss_accum += loss_value
        
        ave_loss = loss_accum / i_step
        train_accuracy = float(correct_samples) / total_samples


        model.eval()
        ave_valid_loss = 0
        with torch.no_grad():
            x = torch.tensor(validationX).to(device)
            y = torch.tensor(validationY)
            y = y.type(torch.LongTensor).to(device)
            prediction = model(x)
            #print(x)
            print(prediction.softmax(dim=1))
            ave_valid_loss = loss(prediction, y)
            _, indices = torch.max(prediction, 1)
            correct_samples = torch.sum(indices == y)
            valid_accuracy = float(correct_samples) / y.shape[0]
        
            #ave_valid_loss = loss_value/len(validationY)
        model.train() 

        if scheduler is not None:
            #scheduler.step(ave_valid_loss)
            scheduler.step()

        loss_history.append(float(ave_loss))
        train_history.append(train_accuracy)
        validLoss_history.append(float(ave_valid_loss))
        ep = np.arange(1,epoch+2,1)
        lv = np.array(validLoss_history)
        lt = np.array(loss_history)
        plt.clf()
        plt.plot(ep,lt,"blue",label="train")
        plt.plot(ep,lv,"orange",label="validation")
        plt.legend(loc=[0.5,0.6])
        plt.xlabel('Epoch')
        plt.ylabel('Loss')
        if ((epoch+1)%100 == 0):
            plt.show()
        #plt.draw()

        print("Average loss: %f, valid loss: %f, Train accuracy: %f, V acc: %f, epoch: %f" % (ave_loss, ave_valid_loss, train_accuracy, valid_accuracy, epoch+1))
    
    return loss_history, train_history


def train_NN(dataTable, nnname, validTable):
    print("start nn training")
    mean, std = meanAndStdTable(dataTable)
    train_dataset = DESY_dataset(dataTable, mean, std)
    validX, validY = load_dataset(validTable, mean, std)
    print("dataset created")
    batch_size = 50

    data_size = len(train_dataset)
    indices = list(range(data_size))
    np.random.shuffle(indices)
    train_indices = indices
    train_sampler = SubsetRandomSampler(train_indices)

    dropLast = False
    if (dataTable.shape[0]%batch_size==1):
        dropLast = True
    train_loader = torch.utils.data.DataLoader(train_dataset, batch_size=batch_size,
                                               sampler=train_sampler, drop_last=dropLast)
    """
    nn_model = nn.Sequential(
        nn.Linear(train_dataset[0][0].shape[0], 256, bias=True),
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
        #nn.Softmax(dim=1),
    )
    """
    nn_model = RNN(input_dim=4,
               embedding_dim=256,
                  hidden_dim=512,
                  output_dim=4 
    ) # could use 1 for binary classification
    nn_model.type(torch.FloatTensor)
    nn_model.to(device)
    print(next(nn_model.parameters()).is_cuda)

    loss = nn.CrossEntropyLoss().type(torch.cuda.FloatTensor)

    #optimizer = optim.SGD(nn_model.parameters(), lr=0.1, momentum=0.9, weight_decay=0.05)
    optimizer = optim.Adam(nn_model.parameters(), lr=0.000005, betas=(0.5, 0.9), weight_decay=0.2)

    scheduler = optim.lr_scheduler.StepLR(optimizer, step_size=10, gamma=0.8)
    #scheduler = optim.lr_scheduler.ReduceLROnPlateau(optimizer, threshold=0.2, factor=0.2)

    print("prepared to train nn")
    train_model(nn_model, train_loader, loss, optimizer, 100, validX, validY, scheduler = scheduler)
    print("trained nn")
    torch.save(nn_model.state_dict(), nnname)


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
        return hidden.to(device)


def writeTrainData(energy,meanArr,stdArr):
    np.savetxt('energies.txt', energy, delimiter="\n", fmt='%s')
    np.savetxt('meanValues.txt', meanArr, fmt='%s')
    np.savetxt('stdValues.txt', stdArr, fmt='%s')


def train_All_NNs():
    #split by energy, split on valid and train subsets, train and save energy subsets' parameters

    fullSetTable = pandas.read_table("trainSet.txt",sep='	',header=None)
    dfs = dict(tuple(fullSetTable.groupby(list(fullSetTable.columns)[10])))
    listdf = [dfs[x] for x in dfs]
    meanArr = []
    stdArr  = []
    energy  = []
    for dft in listdf:
        dftCorr = prepareTable(dft)
        if (len(dftCorr)<100):
            continue
        x2kfColName = list(dftCorr.columns)[8]
        x24pColName = list(dftCorr.columns)[9]
        columnName = list(dftCorr.columns)[10]
        energyLoc = "{:.3f}".format(dftCorr.loc[0].at[columnName])
        #if (float(energyLoc) < 999.999):
        #    continue
        name = 'networks\\model_'+str(energyLoc)+'.pt'
        print(name)
        dftCorr.drop(columnName, axis=1, inplace=True)
        dftCorr.drop([x2kfColName,x24pColName],axis=1, inplace=True)
        dfCorrTrain = dftCorr.sample(frac=0.8).sort_index()
        dfCorrValid = dftCorr.drop(dfCorrTrain.index)
    
        mean, std = meanAndStdTable(dfCorrTrain)
        meanArr.append(mean)
        stdArr.append(std)
        energy.append(energyLoc)

        train_NN(dfCorrTrain,name,dfCorrValid)

    writeTrainData(np.array(energy),np.concatenate(list(meanArr),axis=0),np.concatenate(list(stdArr),axis=0))
    

def train(model):
    if model == "NN":
        train_All_NNs()


print("start_train_python")
train(sys.argv[1])
