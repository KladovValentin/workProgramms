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


def load_dataset(fname,lborder,hborder):
    x = []
    y = []
    local_x = []
    local_y = 0
    countLines=0
    f = io.open(fname, 'r')
    for i, line in enumerate(f):
        local_x = []
        local_y = 0
        tokens = line.strip().split(' ')
        if(float(tokens[0])>=0 and float(tokens[3])>=0 and float(tokens[10])>=lborder/750. and float(tokens[10])<hborder/750.):
            for j in range(8):
                local_x.append(float(tokens[j]))
            local_x.append(float(tokens[8])/20.0)
            local_x.append(float(tokens[8])/20.0-float(tokens[9])/20.0)
            local_x.append(float(tokens[10]))
            local_y = int(tokens[11])

        if local_x != []:
            x.append(local_x)
            y.append(local_y)
    f.close()
    print("data arrays are filled")
    
    x = np.array(x)
    y = np.array(y)
    print(x)
    print(y)
    print("transmuted into np.array")
    
    print("i am here")
    return x.astype(np.float32),  y


class DESY_dataset(Dataset):
    def __init__(self, path,lborder,hborder, transform=None):
        self.transform = transform
        self.datasetX, self.datasetY = load_dataset(path,lborder,hborder)

    def __len__(self):
        return len(self.datasetY)

    def __getitem__(self, index):
        if self.transform != None:
            self.datasetX[index] = self.transform(self.datasetX[index])
        return torch.tensor(self.datasetX[index]), torch.tensor(self.datasetY[index])


def train_model(model, train_loader, loss, optimizer, num_epochs, scheduler=None):
    print("start model nn train")
    loss_history = []
    train_history = []
    for epoch in range(num_epochs):
        model.train()  # Enter train mode

        loss_accum = 0
        correct_samples = 0
        total_samples = 0
        #print(len(train_loader))
        for i_step, (x, y) in enumerate(train_loader):
            x = x.to(device)
            y = y.type(torch.LongTensor)
            y = y.to(device)
            #print(x)
            #print(y)
            prediction = model(x)    
            loss_value = loss(prediction, y)
            optimizer.zero_grad()
            loss_value.backward()
            optimizer.step()
            
            _, indices = torch.max(prediction, 1)
            correct_samples += torch.sum(indices == y)
            total_samples += y.shape[0]
            
            loss_accum += loss_value

        if scheduler is not None:
            scheduler.step()
        ave_loss = loss_accum / i_step
        train_accuracy = float(correct_samples) / total_samples

        loss_history.append(float(ave_loss))
        train_history.append(train_accuracy)

        print("Average loss: %f, Train accuracy: %f, epoch: %f" % (ave_loss, train_accuracy, epoch))

    return loss_history, train_history


def train_NN(lborder,hborder,nnname):
    print("start nn training")
    train_dataset = DESY_dataset("trainSetCorr.txt",lborder,hborder)
    print("dataset created")
    batch_size = 50

    data_size = len(train_dataset)
    indices = list(range(data_size))
    np.random.shuffle(indices)

    train_indices = indices

    train_sampler = SubsetRandomSampler(train_indices)

    train_loader = torch.utils.data.DataLoader(train_dataset, batch_size=batch_size,
                                               sampler=train_sampler)

    nn_model = nn.Sequential(
        nn.Linear(train_dataset[0][0].shape[0], 256, bias=True),
        nn.ReLU(inplace=True),
        nn.Linear(256, 1024),
        #nn.ReLU(inplace=True),
        #nn.Linear(1024, 1024),
        nn.ReLU(inplace=True),
        nn.Linear(1024, 4),
        #nn.Softmax(dim=1),
    )
    nn_model.type(torch.FloatTensor)
    nn_model.to(device)
    print(next(nn_model.parameters()).is_cuda)

    #loss to mimimize - log(normalized probabilities)
    loss = nn.CrossEntropyLoss().type(torch.cuda.FloatTensor)

    #define optimzer like Stochastic gradient descent, which will imptove parameters based on gradient of loss function, lr - learning rate, momentum - some shit (http://www.cs.toronto.edu/~hinton/absps/momentum.pdf)
    optimizer = optim.SGD(nn_model.parameters(), lr=0.01, momentum=0.9)

    #adjust learning rate based on the number of epochs. Just aftr "step_size" epochs "lr" will be multiplied by "gamma"
    scheduler = optim.lr_scheduler.StepLR(optimizer, step_size=15, gamma=0.8)

    print("prepared to train nn")
    train_model(nn_model, train_loader, loss, optimizer, 200, scheduler)
    print("trained nn")
    torch.save(nn_model.state_dict(), nnname)



def train(model):
    if model == "NN":
        train_NN(150,750,'nnMode750.pt')
        train_NN(750,825,'nnMode825.pt')
        train_NN(825,900,'nnMode900.pt')
        train_NN(900,1100,'nnMode1100.pt')



print("start_train_python")
train(sys.argv[1])


"""
def forward(self, x):
    x = self.flatten(x)
    logits = self.linear_relu_stack(x)
    return logits

model = NeuralNetwork().to(device)
print(model)

X = torch.rand(1, 28, 28, device=device)
logits = model(X)
pred_probab = nn.Softmax(dim=1)(logits)
y_pred = pred_probab.argmax(1)
print(f"Predicted class: {y_pred}")
"""