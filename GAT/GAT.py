import networkx as nx
import numpy as np
import torch
from sklearn.preprocessing import StandardScaler

import pandas as pd
from torch_geometric.data import InMemoryDataset, Data
from sklearn.model_selection import train_test_split
import torch_geometric.transforms as T

import torch.nn as nn
import torch.nn.functional as F
from torch_geometric.nn import GCNConv

# load graph from networkx library
# G = nx.karate_club_graph()

# retrieve the labels for each node
labels = np.asarray([G.node[i]['club'] != 'Mr. Hi' for i in G.node]).astype(np.int64)
# labels = np.asarray([x[0].astype(np.int64) for x in pd.read_csv('/home/seongwonhwang/Desktop/projects/mogrify/Statistical Consulting/BIC/data/min_test_node.txt',header=None).values])

# create edge index from 
adj = nx.to_scipy_sparse_matrix(G).tocoo()
row = torch.from_numpy(adj.row.astype(np.int64)).to(torch.long)
col = torch.from_numpy(adj.col.astype(np.int64)).to(torch.long)
edge_index = torch.stack([row, col], dim=0)
# edge_index = pd.read_csv('/home/seongwonhwang/Desktop/projects/mogrify/Statistical Consulting/BIC/data/min_test.txt', header=None, sep='\t')
# edge_index = torch.from_numpy(np.array(edge_index).T)

# using degree as embedding
embeddings = np.array(list(dict(G.degree()).values()))

# normalizing degree values
scale = StandardScaler()
embeddings = scale.fit_transform(embeddings.reshape(-1,1))
# embeddings = np.array([[int(x)] for x in labels])

my_nodes = list(range(34))

# custom dataset
class KarateDataset(InMemoryDataset):
    def __init__(self, transform=None):
        super(KarateDataset, self).__init__('.', transform, None, None)
        data = Data(edge_index=edge_index)
        # data.num_nodes = G.number_of_nodes()
        data.num_nodes = len(labels)
        # embedding 
        data.x = torch.from_numpy(embeddings).type(torch.float32)
        # labels
        y = torch.from_numpy(labels).type(torch.long)
        data.y = y.clone().detach()
        data.num_classes = 2
        # splitting the data into train, validation and test
        # X_train, X_test, y_train, y_test = train_test_split(pd.Series(list(G.nodes())), 
        X_train, X_test, y_train, y_test = train_test_split(pd.Series(my_nodes), 
                                                            pd.Series(labels),
                                                            test_size=0.30, 
                                                            random_state=42)
        # n_nodes = G.number_of_nodes()
        n_nodes = len(labels)
        # create train and test masks for data
        train_mask = torch.zeros(n_nodes, dtype=torch.bool)
        test_mask = torch.zeros(n_nodes, dtype=torch.bool)
        train_mask[X_train.index] = True
        test_mask[X_test.index] = True
        data['train_mask'] = train_mask
        data['test_mask'] = test_mask
        self.data, self.slices = self.collate([data])
    def _download(self):
        return
    def _process(self):
        return
    def __repr__(self):
        return '{}()'.format(self.__class__.__name__)
    
dataset = KarateDataset()
data = dataset[0]

##############################
# Graph Convolutional Network.

# GCN model with 2 layers 
class Net(torch.nn.Module):
    def __init__(self):
        super(Net, self).__init__()
        self.conv1 = GCNConv(data.num_features, 16)
        self.conv2 = GCNConv(16, int(data.num_classes))
    def forward(self):
        x, edge_index = data.x, data.edge_index
        x = F.relu(self.conv1(x, edge_index))
        x = F.dropout(x, training=self.training)
        x = self.conv2(x, edge_index)
        return F.log_softmax(x, dim=1)

device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
data =  data.to(device)
model = Net().to(device) 


##############################
# Train the GCN model.

torch.manual_seed(42)

optimizer_name = "Adam"
lr = 1e-1
optimizer = getattr(torch.optim, optimizer_name)(model.parameters(), lr=lr)
epochs = 200

def train():
  model.train()
  optimizer.zero_grad()
  loss = F.nll_loss(model()[data.train_mask], data.y[data.train_mask])
  loss.backward()
  optimizer.step()

@torch.no_grad()
def test():
  model.eval()
  logits = model()
  mask1 = data['train_mask']
  pred1 = logits[mask1].max(1)[1]
  acc1 = pred1.eq(data.y[mask1]).sum().item() / mask1.sum().item()
  mask = data['test_mask']
  pred = logits[mask].max(1)[1]
  acc = pred.eq(data.y[mask]).sum().item() / mask.sum().item()
  return acc1,acc

for epoch in range(1, epochs):
  train()

train_acc,test_acc = test()

print('#' * 70)
print('Train Accuracy: %s' %train_acc )
print('Test Accuracy: %s' % test_acc)
print('#' * 70)