from PyG_Link_util import *

# 1. Build graphs for training and predicting
path_files = '/home/seongwonhwang/Desktop/projects/GRN_in_general/git/GRN_inference_practice/input_data_processing/data'
graph_for_training = build_graph(path_files, 'training')

# 2. Split the graph into training and validation
split = T.RandomLinkSplit(
    num_val=0.05,
    num_test=0.,
    is_undirected=False,
    add_negative_train_samples=True,
    neg_sampling_ratio=1.0,
)
train_data, val_data, test_data = split(graph_for_training)

# 3. Build a model and train
device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
model = Net(train_data.x.shape[1], 128, 64).to(device)
optimizer = torch.optim.Adam(params=model.parameters(), lr=0.002)
criterion = torch.nn.BCEWithLogitsLoss()
model = train_link_predictor(
    model, train_data, val_data, optimizer, criterion, n_epochs=5000)

# 4. Save the final predicted model to a file
get_network(path_files, model)
