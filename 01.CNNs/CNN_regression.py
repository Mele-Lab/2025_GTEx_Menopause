##CNN script with early stopping and fine tunning
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
from torch.utils.data import Dataset, DataLoader

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
workdir = '/gpfs/projects/bsc83/Projects/GTEx_v8/Laura/05.CNN/RegressionModels/'
tissue = 'Ovary'
#savedir = '/gpfs/projects/bsc83/Projects/GTEx_v8/Ole/'

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
    log_file = f"{workdir}{tissue}/{lr}_{gamma}_{delta}_{drop}_{patience}_training_log_REGRESSION.csv"
    os.makedirs(os.path.dirname(log_file), exist_ok=True) 
    # filter = pd.read_csv("BSC/03.Image_processing/Second_filtering_images/Uterus_filtered_myometrium_images.csv")

    with open(log_file, mode="w", newline="") as file:
        writer = csv.writer(file)
        writer.writerow(["Epoch", "Phase", "Tile Loss", "Tile MAE", "Donor MAE"])
    
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
    

    class AgeRegressionDataset(Dataset):
        def __init__(self, csv_file, img_dir, transform=None):
            self.labels_df = pd.read_csv(csv_file)
            self.img_dir = img_dir
            self.transform = transform

        def __len__(self):
            return len(self.labels_df)

        def __getitem__(self, idx):
            img_name = self.labels_df.iloc[idx, 0]
            age = self.labels_df.iloc[idx, 1]

            img_path = os.path.join(self.img_dir, img_name)
            image = Image.open(img_path).convert("RGB")

            if self.transform:
                image = self.transform(image)

            age = torch.tensor(age, dtype=torch.float32)

            return image, age, img_name  # <-- return image name too
            
    train_dataset = AgeRegressionDataset(csv_file=workdir+tissue+'/labels_train.csv', img_dir=train_directory, transform=image_transforms["train"])
    test_dataset = AgeRegressionDataset(csv_file=workdir+tissue+'/labels_test.csv',img_dir=test_directory,transform=image_transforms["test"])
    
    

    # Dataloaders
    dataloaders = {
        'train': DataLoader(train_dataset, batch_size=bs, shuffle= True, num_workers=10), 
        'test': DataLoader(test_dataset, batch_size=bs, shuffle=False, num_workers=10)
    }
    
    
  
    # Load pre-trained model
    model = torch.load('/gpfs/projects/bsc83/Projects/GTEx_v8/Laura/05.CNN/pretrained_vgg19_bn.pt')
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
        torch.nn.Linear(1024, 1)  # single continuous output
    )
    
    # Move model to GPU if available
    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    model = model.to(device)

    loss_func = torch.nn.MSELoss()
    optimizer = torch.optim.RMSprop(model.parameters(), lr=lr)
    lr_scheduler = torch.optim.lr_scheduler.StepLR(optimizer, step_size=5, gamma=gamma)
    epochs = 100

    dataset_sizes = {'train': len(train_dataset), 'test': len(test_dataset)}

    print(f"Train size: {dataset_sizes['train']} | Test size: {dataset_sizes['test']}")


    # --------------------------
    # Early Stopping
    # --------------------------

    class EarlyStopping:
        def __init__(self, patience=15, min_delta=0.001):
            self.patience = patience
            self.min_delta = min_delta
            self.counter = 0
            self.best_loss = None
            self.early_stop = False

        def __call__(self, val_loss):
            if self.best_loss is None:
                self.best_loss = val_loss
            elif val_loss > self.best_loss - self.min_delta:
                self.counter += 1
                if self.counter >= self.patience:
                    self.early_stop = True
            else:
                self.best_loss = val_loss
                self.counter = 0

    early_stopping = EarlyStopping(patience=patience, min_delta=0.001)

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
                
            aggregated_true_labels.append(true_label_list[0])
            aggregated_predicted_labels.append(np.mean(predicted_label_list))
            aggregated_patterns.append(pattern)
    
        return aggregated_true_labels, aggregated_predicted_labels
    
    # --------------------------
    # Training function
    # --------------------------

    def train_model(model, criterion, optimizer, lr_scheduler, dataloaders, dataset_sizes, device, num_epochs=epochs):
        since = time.time()
        
        best_model_wts = model.state_dict()
        best_donor_mae = float("inf")

        metrics = {"train": {"loss": [], "mae": [], "donor_mae": []}, 
                "test": {"loss": [], "mae": [], "donor_mae": []}}

        for epoch in range(num_epochs):
            print(f"\nEpoch {epoch+1}/{num_epochs}")
            print("-" * 20)

            for phase in ['train', 'test']:
                if phase == 'train':
                    model.train()
                else:
                    model.eval()

                running_loss = 0.0
                running_mae = 0.0

                # To store for donor aggregation (only test phase)
                all_preds = []
                all_labels = []
                all_paths = []

                for inputs, labels, paths in dataloaders[phase]:
                    inputs = inputs.to(device)
                    labels = labels.to(device).unsqueeze(1)  # (batch, 1)

                    optimizer.zero_grad()

                    with torch.set_grad_enabled(phase == 'train'):
                        outputs = model(inputs)
                        loss = criterion(outputs, labels)
                        mae = torch.mean(torch.abs(outputs - labels))

                        if phase == 'train':
                            loss.backward()
                            optimizer.step()

                    running_loss += loss.item() * inputs.size(0)
                    running_mae += mae.item() * inputs.size(0)
                    if phase in ['train', 'test']:
                        all_preds.extend(outputs.detach().cpu().numpy().flatten())
                        all_labels.extend(labels.cpu().numpy().flatten())
                        all_paths.extend(paths)

                epoch_loss = running_loss / dataset_sizes[phase]
                epoch_mae = running_mae / dataset_sizes[phase]
                metrics[phase]["loss"].append(epoch_loss)
                metrics[phase]["mae"].append(epoch_mae)

                # Aggregate by donor for BOTH train and test
                true_donor, pred_donor = aggregate_by_pattern(all_labels, all_preds, all_paths)
                epoch_donor_mae = np.mean(np.abs(np.array(true_donor) - np.array(pred_donor)))
                metrics[phase]["donor_mae"].append(epoch_donor_mae)

                print(f"{phase} Tile Loss: {epoch_loss:.4f} | Tile MAE: {epoch_mae:.4f} | Donor MAE: {epoch_donor_mae:.4f}")
                
                with open(log_file, mode="a", newline="") as file:
                    writer = csv.writer(file)
                    writer.writerow([epoch + 1, phase, epoch_loss, epoch_mae, epoch_donor_mae])

                if phase == 'test':
                    # Early stopping on donor-level MAE
                    if epoch_donor_mae < best_donor_mae:
                        best_donor_mae = epoch_donor_mae
                        best_model_wts = model.state_dict()
                    early_stopping(epoch_donor_mae)


                    if early_stopping.early_stop:
                        print("Early stopping triggered!")
                        break
                else:
                    print()  # newline after train phase

            if lr_scheduler is not None:
                lr_scheduler.step()
            print()
    
            if early_stopping.early_stop:
                print("Early stopping")
                break


        time_elapsed = time.time() - since
        print('Training complete in {:.0f}m {:.0f}s'.format(time_elapsed // 60, time_elapsed % 60))
        print('Best test donor MAE: {:.4f}'.format(best_donor_mae))

        model.load_state_dict(best_model_wts)
        return model


            
    # Train the model
    trained_model = train_model(model, loss_func, optimizer, lr_scheduler, dataloaders, dataset_sizes, device, num_epochs=epochs)

    # Save the trained model
    torch.save(trained_model.state_dict(), f"{workdir}{tissue}/regression_model_512_{lr}_{gamma}_{delta}_{drop}_{patience}_{tissue}_trained_model_vgg19_bn_1024RSM.pt")



