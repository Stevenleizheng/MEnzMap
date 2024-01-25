import torch
from torch import nn
from torch.nn import functional as F

# enzyme and non_enzyme
class MLP_esm3B_1(nn.Module):
    def __init__(self):
        super().__init__()
        self.hidden0 = nn.Linear(2560,1920)
        self.hidden1 = nn.Linear(1920,720)
        self.hidden2 = nn.Linear(720,100)
        self.hidden3 = nn.Linear(100,20)
        self.out = nn.Linear(20,2)
    
    def forward(self, X):
        x = F.relu(self.hidden0(X))
        x = F.relu(self.hidden1(x))
        x = F.relu(self.hidden2(x))
        x = F.relu(self.hidden3(x))
        net = self.out(x)
        return net

# enzyme first
class MLP_esm3B_2(nn.Module):
    def __init__(self):
        super().__init__()
        self.hidden0 = nn.Linear(2560,1920)
        self.hidden1 = nn.Linear(1920,720)
        self.hidden2 = nn.Linear(720,360)
        self.hidden3 = nn.Linear(360,120)
        self.out = nn.Linear(120,7)
    
    def forward(self, X):
        x = F.relu(self.hidden0(X))
        x = F.relu(self.hidden1(x))
        x = F.relu(self.hidden2(x))
        x = F.relu(self.hidden3(x))
        net = self.out(x)
        return net

# enzyme second first(1.1)
class MLP_esm3B_3(nn.Module):
    def __init__(self):
        super().__init__()
        self.hidden0 = nn.Linear(2560,1920)
        self.hidden1 = nn.Linear(1920,720)
        self.hidden2 = nn.Linear(720,360)
        self.hidden3 = nn.Linear(360,120)
        self.out = nn.Linear(120,21)
    
    def forward(self, X):
        x = F.relu(self.hidden0(X))
        x = F.relu(self.hidden1(x))
        x = F.relu(self.hidden2(x))
        x = F.relu(self.hidden3(x))
        net = self.out(x)
        return net

# enzyme second second(2.1)
class MLP_esm3B_4(nn.Module):
    def __init__(self):
        super().__init__()
        self.hidden0 = nn.Linear(2560,1920)
        self.hidden1 = nn.Linear(1920,720)
        self.hidden2 = nn.Linear(720,360)
        self.hidden3 = nn.Linear(360,120)
        self.out = nn.Linear(120,9)
    
    def forward(self, X):
        x = F.relu(self.hidden0(X))
        x = F.relu(self.hidden1(x))
        x = F.relu(self.hidden2(x))
        x = F.relu(self.hidden3(x))
        net = self.out(x)
        return net

# enzyme second third(3.1)
class MLP_esm3B_5(nn.Module):
    def __init__(self):
        super().__init__()
        self.hidden0 = nn.Linear(2560,1920)
        self.hidden1 = nn.Linear(1920,720)
        self.hidden2 = nn.Linear(720,360)
        self.hidden3 = nn.Linear(360,90)
        self.out = nn.Linear(90,6)
    
    def forward(self, X):
        x = F.relu(self.hidden0(X))
        x = F.relu(self.hidden1(x))
        x = F.relu(self.hidden2(x))
        x = F.relu(self.hidden3(x))
        net = self.out(x)
        return net

# enzyme second fourth(4.1)
class MLP_esm3B_6(nn.Module):
    def __init__(self):
        super().__init__()
        self.hidden0 = nn.Linear(2560,1920)
        self.hidden1 = nn.Linear(1920,720)
        self.hidden2 = nn.Linear(720,360)
        self.hidden3 = nn.Linear(360,90)
        self.out = nn.Linear(90,7)
    
    def forward(self, X):
        x = F.relu(self.hidden0(X))
        x = F.relu(self.hidden1(x))
        x = F.relu(self.hidden2(x))
        x = F.relu(self.hidden3(x))
        net = self.out(x)
        return net

# enzyme second fiveth(5.1)
class MLP_esm3B_7(nn.Module):
    def __init__(self):
        super().__init__()
        self.hidden0 = nn.Linear(2560,1920)
        self.hidden1 = nn.Linear(1920,720)
        self.hidden2 = nn.Linear(720,360)
        self.hidden3 = nn.Linear(360,90)
        self.out = nn.Linear(90,6)
    
    def forward(self, X):
        x = F.relu(self.hidden0(X))
        x = F.relu(self.hidden1(x))
        x = F.relu(self.hidden2(x))
        x = F.relu(self.hidden3(x))
        net = self.out(x)
        return net

# enzyme second sixth(6.1)
class MLP_esm3B_8(nn.Module):
    def __init__(self):
        super().__init__()
        self.hidden0 = nn.Linear(2560,1920)
        self.hidden1 = nn.Linear(1920,720)
        self.hidden2 = nn.Linear(720,360)
        self.hidden3 = nn.Linear(360,60)
        self.out = nn.Linear(60,6)
    
    def forward(self, X):
        x = F.relu(self.hidden0(X))
        x = F.relu(self.hidden1(x))
        x = F.relu(self.hidden2(x))
        x = F.relu(self.hidden3(x))
        net = self.out(x)
        return net

# enzyme second seventh(7.1)
class MLP_esm3B_9(nn.Module):
    def __init__(self):
        super().__init__()
        self.hidden0 = nn.Linear(2560,1920)
        self.hidden1 = nn.Linear(1920,720)
        self.hidden2 = nn.Linear(720,360)
        self.hidden3 = nn.Linear(360,90)
        self.out = nn.Linear(90,6)
    
    def forward(self, X):
        x = F.relu(self.hidden0(X))
        x = F.relu(self.hidden1(x))
        x = F.relu(self.hidden2(x))
        x = F.relu(self.hidden3(x))
        net = self.out(x)
        return net