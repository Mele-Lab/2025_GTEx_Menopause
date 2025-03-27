
##CNN training with early stopping and fine-tuning

import matplotlib.pyplot as plt
import matplotlib.image as mpimg
import numpy as np
import os
import pandas as pd
import seaborn as sns
import subprocess as sub
import torch
import torchvision
import random
import time
import io
import re
import csv
from collections import defaultdict
from PIL import Image
from sklearn.metrics import confusion_matrix
from torch.utils.data.sampler import SubsetRandomSampler, WeightedRandomSampler
from torchvision import datasets, models, transforms

from sklearn.metrics import roc_curve
from sklearn.metrics import auc
from sklearn.metrics import balanced_accuracy_score
from sklearn.metrics import RocCurveDisplay

# Set random seeds for reproducibility
seed = 180
random.seed(seed)
np.random.seed(seed)
torch.manual_seed(seed)
torch.backends.cudnn.deterministic = True
torch.backends.cudnn.benchmark = False

# Paths and directories
workdir = '/X/Laura/05.CNN/35yo_separation_FEMALE/'
tissue = 'Uterus'

# Set train and test paths
train_directory = workdir+tissue+'/Tiles512_05/train/'
test_directory = workdir+tissue+'/Tiles512_05/test/'

gamma_values = [0.3,0.4, 0.5, 0.6,0.7, 0.8]
drop_values = [0.3,0.4, 0.5]

lr = 0.001
patience = 15
delta = 0.001

from itertools import product
param_combinations = list(product(gamma_values, drop_values))

# Loop through combinations
for gamma, drop in param_combinations:
    print(f"Running with: lr={lr}, gamma={gamma}, drop={drop}, patience={patience}")
    log_file = f"{workdir}{tissue}/{lr}_{gamma}_{delta}_{drop}_{patience}_training_logACC.csv"
    os.makedirs(os.path.dirname(log_file), exist_ok=True) 

    with open(log_file, mode="w", newline="") as file:
        writer = csv.writer(file)
        writer.writerow(["Epoch","Phase", "Epoch loss", "Epoch acc", "Epoch bal acc"])

    # Define a set of transformations to augment the data
    image_transforms = {
        'train': transforms.Compose([
           # transforms.RandomResizedCrop(size=224, scale=(0.8, 1.0)),
          #  transforms.RandomRotation(degrees=15),
            transforms.Resize(224),
            transforms.ToTensor(),
            transforms.Normalize([0.485, 0.456, 0.406],
                                 [0.229, 0.224, 0.225])
        ]),
        'test': transforms.Compose([
            transforms.Resize(224),
            transforms.ToTensor(),
            transforms.Normalize([0.485, 0.456, 0.406],
                                 [0.229, 0.224, 0.225])
        ])
    }
    
    # Batch size
    bs = 64
    num_classes = 2
    
    # Custom dataset to include image file paths
    class ImageFolderWithPaths(torchvision.datasets.ImageFolder):
        def __getitem__(self, index):
            original_tuple = super(ImageFolderWithPaths, self).__getitem__(index)
            path = self.imgs[index][0]
            tuple_with_path = (original_tuple + (path,))
            return tuple_with_path
    
    # Load Data from folders
    data = {
        'train': ImageFolderWithPaths(root=train_directory, transform=image_transforms['train']),
        'test': ImageFolderWithPaths(root=test_directory, transform=image_transforms['test'])
    }
    
    # Size of Data, to be used for calculating Average Loss and Accuracy
    dataset_sizes = {
        'train': len(data['train']),
        'test': len(data['test'])
    }
    
    # Class distribution in train and test sets
    dataset_distributions = {
        'train': np.unique(data['train'].targets, return_counts=True)[1],
        'test': np.unique(data['test'].targets, return_counts=True)[1]
    }
    
    # Print dataset sizes and class distributions
    print('Sizes, train: {} test: {}'.format(dataset_sizes['train'], dataset_sizes['test']))
    print('train, young: {} old: {}'.format(dataset_distributions['train'][1], dataset_distributions['train'][0]))
    print('test, young: {} old: {}'.format(dataset_distributions['test'][1], dataset_distributions['test'][0]))
    
    # Load training metadata
    # train_sample_meta = pd.read_csv(workdir + tissue + '/comparison1_3_35yo' + tissue + '.csv')
    # train_sample_meta = train_sample_meta.sort_values(['phenotype']).reset_index(drop=True)
    
    # Compute class weights for the weighted sampler
    class_weights = 1. / torch.tensor(dataset_distributions['train'], dtype=torch.float) #calculate the inverse of the number of instances per class, so the class with higher number of samples has a lower weight
    sample_weights = class_weights[data['train'].targets] #assign to each sample its weight
    weighted_sampler = WeightedRandomSampler(weights=sample_weights, num_samples=len(sample_weights), replacement=True) #number of samples= size of the dataset, so each epoh will iterate over all the samples, but with
                                                                                                                        #probabilities defined by their weights, replacement=TRUE allows to take the samples from the minority class more often
    
    # Dataloaders
    dataloaders = {
        'train': torch.utils.data.DataLoader(data['train'], batch_size=bs, sampler=weighted_sampler, num_workers=10), #no need of shuffling because it is inherent of the WeightRandomSampler
        'test': torch.utils.data.DataLoader(data['test'], batch_size=bs, shuffle=True, num_workers=10)
    }
    
    # Load pre-trained model
    model = torch.load('/X/Laura/05.CNN/pretrained_vgg19_bn.pt')
    model.eval()
    
    # Freeze model parameters
    for param in model.parameters():
        param.requires_grad = False
    
    # Modify the final fully connected layer for transfer learning
    fc_inputs = model.classifier[6].in_features
    model.classifier[6] = torch.nn.Sequential(
        torch.nn.Linear(fc_inputs, 1024),
        torch.nn.ReLU(),
        torch.nn.Dropout(drop),
        torch.nn.Linear(1024, num_classes),
        torch.nn.LogSoftmax(dim=1) # For using NLLLoss()
    )
    
    # Move model to GPU if available
    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    model = model.to(device)
    
    # Define optimizer, loss function, learning rate scheduler, and epochs
    
    loss_func = torch.nn.NLLLoss()
    optimizer = torch.optim.RMSprop(model.parameters(), lr=0.001)
    lr_scheduler = torch.optim.lr_scheduler.StepLR(optimizer, step_size=5, gamma=gamma)
    epochs = 100
    
    # Function to extract the pattern
    def extract_pattern(filepath):
        # Updated regex to find the pattern GTEX- followed by letters and digits until the end or next separator
        match = re.search(r'GTEX-\w+', filepath)
        return match.group(0) if match else None
    
    def aggregate_by_pattern(true_labels, predicted_labels, filepaths):
        # Dictionary to hold lists of true labels and predicted labels for each pattern
        patterns = [extract_pattern(filepath) for filepath in filepaths]
        
        pattern_dict = defaultdict(lambda: {'true_labels': [], 'predicted_labels': []})
    
        # Populate the dictionary
        for true_label, predicted_label, pattern in zip(true_labels, predicted_labels, patterns):
            pattern_dict[pattern]['true_labels'].append(true_label)
            pattern_dict[pattern]['predicted_labels'].append(predicted_label)
    
        # Lists to hold the aggregated results
        aggregated_true_labels = []
        aggregated_predicted_labels = []
        aggregated_patterns = []
    
        # Aggregate the values
        for pattern, values in pattern_dict.items():
            true_label_list = values['true_labels']
            predicted_label_list = values['predicted_labels']
    
            # Ensure all true labels are the same for a given pattern
            assert all(label == true_label_list[0] for label in true_label_list), "Inconsistent true labels for pattern " + pattern
            
            aggregated_true_labels.append(true_label_list[0])
            aggregated_predicted_labels.append(np.mean(predicted_label_list))
            aggregated_patterns.append(pattern)
    
        return aggregated_true_labels, aggregated_predicted_labels
    
    def balanced_accuracy_score_donor(true_label, preds_npy, paths):
        true_labels_aggregated, preds_npy_aggregated = aggregate_by_pattern(true_label, preds_npy, paths)
        threshold_e = 0.5
        preds_npy_aggregated = np.array(preds_npy_aggregated)
        y_pred_binary = (preds_npy_aggregated >= threshold_e).astype(int)
        ba = balanced_accuracy_score(true_labels_aggregated, y_pred_binary)
        return ba
      
    #Early stopping implementation
    class EarlyStopping:
        def __init__(self, patience=15, min_delta=0.001):
            self.patience = patience #number of epochs to wait for an improvement in the validation loss before stopping the training.
            self.min_delta = min_delta #the minimum change in the validation loss to qualify as an improvement
            self.counter = 0 # counter to keep track of how many consecutive epochs have passed without improvement
            self.best_acc = None #the best observed validation loss. Initialized to None and updated during training
            self.early_stop = False #A flag to indicate whether training should be stopped early.
                
        def __call__(self, val_acc): 
            if self.best_acc is None: # for the first iteration
                self.best_acc = val_acc
            elif val_acc < self.best_acc + self.min_delta: 
                self.counter += 1 
                if self.counter >= self.patience: #if we reach the patience, we stop
                    self.early_stop = True
            else:
                self.best_acc = val_acc #if there is a significant improvement, we reset the counter to 0
                self.counter = 0
    
    early_stopping = EarlyStopping(patience=15, min_delta=0.001)
    
    # Training function
    def train_model(model, criterion, optimizer, scheduler, datasets, dataset_sizes, device, num_epochs=epochs):
        since = time.time()
        files_path = np.empty((0))
    
        best_model_wts = model.state_dict()
        best_acc = 0.0
        best_epoch_loss = 1.0
        metrics = {"train": {"loss": [], "acc": [], "balanced_acc": [], "balanced_acc_t": [],"tnr": [], "fpr": [], "fnr": [], "tpr": [], "auc": []}, 
                   "test": {"loss": [], "acc": [], "balanced_acc": [],"balanced_acc_t": [],  "tnr": [], "fpr": [], "fnr": [], "tpr": [], "auc": []}}
    
        for epoch in range(num_epochs):
            print("Epoch: {}/{}".format(epoch+1, num_epochs))
            print('-' * 10)
        
            
            since_epoch = time.time()
    
            for phase in ['train', 'test']:
                if phase == 'train':
                    model.train()
                else:
                    model.eval()
    
                running_loss = 0
                running_corrects = 0
                true_labels = np.empty((0), int)
                preds_npy = np.empty((0), int)
                paths_all =  np.empty((0), str)
    
                for i, (inputs, labels, paths) in enumerate(dataloaders[phase]):
                    inputs = inputs.to(device)
                    labels = labels.to(device)
                   
                    true_labels = np.append(true_labels, labels.data.cpu().numpy(), axis=0)
                    paths_all= np.append(paths_all, paths, axis=0)
    
                    optimizer.zero_grad()
    
                    with torch.set_grad_enabled(phase == 'train'):
                        outputs = model(inputs)
                        scores, preds = torch.max(outputs, 1)
                        preds_npy = np.append(preds_npy, preds.data.cpu().numpy(), axis=0)
                        loss = criterion(outputs, labels)
    
                        if phase == 'train':
                            loss.backward()
                            optimizer.step()
    
                    running_loss += loss.item() * inputs.size(0)
                    running_corrects += torch.sum(preds == labels.data)
    
                epoch_loss = running_loss / dataset_sizes[phase]
                epoch_acc = running_corrects.double() / dataset_sizes[phase]
                epoch_balanced_acc_tile = balanced_accuracy_score(true_labels, preds_npy)
                epoch_balanced_acc = balanced_accuracy_score_donor(true_labels, preds_npy, paths_all)
    
                metrics[phase]["loss"].append(epoch_loss)
                metrics[phase]["acc"].append(epoch_acc.cpu().numpy().item())
                metrics[phase]["balanced_acc"].append(epoch_balanced_acc)
                metrics[phase]["balanced_acc_t"].append(epoch_balanced_acc_tile)
    
                print('{} Loss: {:.4f} Acc: {:.4f} Bal.acc: {:.4f}'.format(
                    phase, epoch_loss, epoch_acc, epoch_balanced_acc))
                    
                with open(log_file, mode="a", newline="") as file:
                    writer = csv.writer(file)
                    writer.writerow([epoch + 1, phase, epoch_loss, epoch_acc, epoch_balanced_acc, epoch_balanced_acc_tile])
                
                #Check if new acc is higher that the previous best acc
                if phase == 'test':
                   if epoch_acc > best_acc and epoch > 8:

                #Monitor the improvement of per donor accuracy
                       print("better accuracy")
                       print(epoch_acc)
                       best_acc = epoch_acc
                       best_model_wts = model.state_dict()
                       
                   ##early stopping based on per tile accuracy not improving
                   if epoch > 8:
                       early_stopping(epoch_acc)
                      # Log to CSV
    
    
            if scheduler is not None:
                scheduler.step()
            print()
    
            if early_stopping.early_stop:
                print("Early stopping")
                break
    
        time_elapsed = time.time() - since
        print('Training complete in {:.0f}m {:.0f}s'.format(
            time_elapsed // 60, time_elapsed % 60))
        print('Best test Bal.acc: {:4f}'.format(best_acc))
    
    
        model.load_state_dict(best_model_wts)
        return model
      
    
    # Train the model
    trained_model = train_model(model, loss_func, optimizer, lr_scheduler, dataloaders, dataset_sizes, device, num_epochs=epochs)
    
    # Save the trained model
    torch.save(trained_model.state_dict(), f"{workdir}{tissue}/ACCtile_filtering_hybrid_tile_512_{lr}_{gamma}_{delta}_{drop}_{patience}_{tissue}_trained_model_vgg19_bn_1024RSM.pt")
