import os
os.chdir('/home/seongwonhwang/Desktop/projects/git/GRN_inference_practice/')
from GAT.GAT_util import *

ID = 'NodeAb4'
neg_ratio = 5.0

for rng_seed in (111, 123, 1234):
    # Set a fixed seed for reproducibility
    random.seed(rng_seed)
    torch.manual_seed(rng_seed)
    np.random.seed(rng_seed)
    # TEST_ID = 'NodeAb1_control_Node'
    for status in ('control', 'ablation'):
        for cell_type in ('Early_endoderm', 'Node', "Mesoderm", "Prospective_neural_plate", "Area_opaca", "Mesodermal_neural", "Non_neural_ectoderm", "Anterior_ingressing_streak", "Germina_crescent", "PGCs"):
            TEST_ID = ID+'_'+status+'_'+cell_type
            # 1. Build graphs for training and predicting
            # path_files = '/home/seongwonhwang/Desktop/projects/git/GRN_inference_practice/input_data_processing/data'
            path_files = '/home/seongwonhwang/Desktop/projects/git/Node_ablation_practice/GRN/data'
            graph_for_training = build_graph(path_files, TEST_ID, 'training')
            # 2. Split the graph into training and validation
            split = T.RandomLinkSplit(
                num_val=0.05,
                num_test=0.,
                is_undirected=False,
                add_negative_train_samples=True,
                neg_sampling_ratio=neg_ratio,
            )
            train_data, val_data, test_data = split(graph_for_training)
            # 3. Build a model and train
            device = torch.device(
                'cuda' if torch.cuda.is_available() else 'cpu')
            model = Net(train_data.x.shape[1], 128, 64).to(device)
            optimizer = torch.optim.Adam(params=model.parameters(), lr=0.002)
            criterion = torch.nn.BCEWithLogitsLoss()
            model = train_link_predictor(
                model, train_data, val_data, optimizer, criterion, TEST_ID, rng_seed, neg_ratio, n_epochs=5000)
            # 4. Save the final predicted model to a file
            get_network(path_files, TEST_ID, rng_seed, neg_ratio)
            # > png('NodeAb1_control_Node_prediction_score_distribution.png', width=500, height=500, res=150)
            # > df = read.table('NodeAb1_control_Node_prediction_score.txt')
            # > hist(df[,3],breaks=100, xlab='score', main='Score distribution')
            # > dev.off()
