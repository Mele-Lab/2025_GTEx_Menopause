#Code to apply temperature scaling to an already pretrained model
import numpy as np
import torch
import torch.nn as nn
import torch.nn.functional as F
from torch.optim import LBFGS
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
import numpy as np
import os
import pandas as pd
import seaborn as sns
import torch
import torchvision
import random
from PIL import Image
from torch.utils.data import DataLoader
from torchvision import datasets, models, transforms
from torch import nn
from torch.optim import LBFGS
import torch.nn.functional as F



# Set seed and device
seed = 180
random.seed(seed)
np.random.seed(seed)
torch.manual_seed(seed)
torch.backends.cudnn.deterministic = True
torch.backends.cudnn.benchmark = False
device = torch.device("cuda" if torch.cuda.is_available() else "cpu")

# --- Temperature Scaling Wrapper ---
class ModelWithTemperature(nn.Module):
    def __init__(self, model):
        super().__init__()
        self.model = model
        self.temperature = nn.Parameter(torch.ones(1) * 1.5)

    def forward(self, input):
        logits = self.model(input)
        return logits / self.temperature

    def temperature_scale(self, logits):
        return logits / self.temperature

# --- Temperature Tuning Function ---
def tune_temperature(model_with_temp, val_loader, device):
    model_with_temp.eval()
    logits_list, labels_list = [], []

    with torch.no_grad():
        for inputs, labels, _ in val_loader:
            inputs = inputs.to(device)
            labels = labels.to(device)
            logits = model_with_temp.model(inputs)# get raw logits from original model (no temp scaling)
            logits_list.append(logits)
            labels_list.append(labels)

    logits = torch.cat(logits_list)
    labels = torch.cat(labels_list)

    criterion = nn.CrossEntropyLoss()
    optimizer = LBFGS([model_with_temp.temperature], lr=0.01, max_iter=50)

    def eval():
        optimizer.zero_grad()
        loss = criterion(model_with_temp.temperature_scale(logits), labels)
        loss.backward()
        return loss

    optimizer.step(eval)
    return model_with_temp

# --- ImageFolder with file paths ---
class ImageFolderWithPaths(torchvision.datasets.ImageFolder):
    def __getitem__(self, index):
        original_tuple = super().__getitem__(index)
        path = self.imgs[index][0]
        return original_tuple + (path,)

# --- Main Execution ---
workdir = 'X/05.CNN/35yo_separation_FEMALE/'
tissues = ['Uterus', 'Ovary', 'Vagina', 'BreastMammaryTissue']
tests = ['test', 'young_post', 'middle_age', 'validationset', 'train']

for tissue in tissues:
    test = 'validationset'
    if tissue == 'Uterus':
        
        drop=0.5#for uterus
        #drop=0.3 #for myometrium
        train_dir = workdir + tissue + '/Tiles512_05/train/' 
        test_dir = workdir + tissue +'/Tiles512_05/'+ test +'/' 
        out_file = workdir+tissue+'/temp_scaling/hybridtrain'+ test +'_512_tile_temp_scaled'+ tissue +'.csv'
        model_path = workdir+tissue+'/selected_tile_ACC/hybrid_tile_512_Uterus_trained_model_vgg19_bn_1024RSM.pt'

    elif tissue == 'Ovary':
        drop=0.4
        train_dir = workdir + tissue + '/Tiles512_05/train/' 
        test_dir = workdir + tissue +'/Tiles512_05/'+ test +'/' 
        out_file = workdir+tissue+'/temp_scaling/hybridtrain'+ test +'_512_tile_temp_scaled'+ tissue +'.csv'
        model_path = workdir+tissue+'/selected_tile_ACC/hybrid_tile_512_Ovary_trained_model_vgg19_bn_1024RSM.pt'


    elif tissue == 'Vagina':
        drop=0.5
        train_dir = workdir + tissue + '/Tiles256_05/train/' 
        test_dir = workdir + tissue +'/Tiles256_05/'+ test +'/' 
        out_file = workdir+tissue+'/temp_scaling/hybridtrain'+ test +'_256_tile_temp_scaled'+ tissue +'.csv'
        model_path = workdir+tissue+'/selected_tile_ACC/hybrid_tile_256_Vagina_trained_model_vgg19_bn_1024RSM.pt'


    elif tissue == 'BreastMammaryTissue':
        drop=0.3
        train_dir = workdir + tissue + '/Tiles512_025/train/' 
        test_dir = workdir + tissue +'/Tiles512_025/'+ test +'/' 
        out_file = workdir+tissue+'/temp_scaling/hybridtrain'+ test +'_512_tile_temp_scaled'+ tissue +'.csv'
        model_path = workdir+tissue+'/selected_tile_ACC/hybrid_tile_512_BreastMammaryTissue_trained_model_vgg19_bn_1024RSM.pt'

    num_classes = 2
    bs = 64

    # Image transforms
    transforms_dict = {
        'train': transforms.Compose([
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

    # Load data
    data = {
        'train': ImageFolderWithPaths(train_dir, transform=transforms_dict['train']),
        'test': ImageFolderWithPaths(test_dir, transform=transforms_dict['test']),
    }

    dataloaders_eval = {
        'train': DataLoader(data['train'], batch_size=bs, shuffle=False),
        'test': DataLoader(data['test'], batch_size=bs, shuffle=False),
    }

    # Load pretrained model and modify
    model = torch.load('X/05.CNN/pretrained_vgg19_bn.pt')
    fc_inputs = model.classifier[6].in_features
    model.classifier[6] = nn.Sequential(
        nn.Linear(fc_inputs, 1024),
        nn.ReLU(),
        nn.Dropout(0.3),
        nn.Linear(1024, num_classes)  # logits output
    )
    model.load_state_dict(torch.load(model_path))
    model.eval()

    # Convert model to be used on GPU,
    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    model = model.to(device)

    # Wrap with temperature scaling
    model_with_temp = ModelWithTemperature(model).to(device)
    model_with_temp = tune_temperature(model_with_temp, dataloaders_eval['test'], device)
    print(f"Tuned temperature: {model_with_temp.temperature.item():.4f}")

    for test in tests:
        # Run predictions on test
        if tissue == 'Uterus':
            drop=0.5#for uterus
            #drop=0.3 #for myometrium
            train_dir = workdir + tissue + '/Tiles512_05/train/' 
            test_dir = workdir + tissue +'/Tiles512_05/'+ test +'/' 
            out_file = workdir+tissue+'/temp_scaling/hybridtrain'+ test +'_512_tile_temp_scaled'+ tissue +'.csv'
            model_path = workdir+tissue+'/selected_tile_ACC/hybrid_tile_512_Uterus_trained_model_vgg19_bn_1024RSM.pt'

        elif tissue == 'Ovary':
            drop=0.4
            train_dir = workdir + tissue + '/Tiles512_05/train/' 
            test_dir = workdir + tissue +'/Tiles512_05/'+ test +'/' 
            out_file = workdir+tissue+'/temp_scaling/hybridtrain'+ test +'_512_tile_temp_scaled'+ tissue +'.csv'
            model_path = workdir+tissue+'/selected_tile_ACC/hybrid_tile_512_Ovary_trained_model_vgg19_bn_1024RSM.pt'


        elif tissue == 'Vagina':
            drop=0.5
            train_dir = workdir + tissue + '/Tiles256_05/train/' 
            test_dir = workdir + tissue +'/Tiles256_05/'+ test +'/' 
            out_file = workdir+tissue+'/temp_scaling/hybridtrain'+ test +'_256_tile_temp_scaled'+ tissue +'.csv'
            model_path = workdir+tissue+'/selected_tile_ACC/hybrid_tile_256_Vagina_trained_model_vgg19_bn_1024RSM.pt'


        elif tissue == 'BreastMammaryTissue':
            drop=0.3
            train_dir = workdir + tissue + '/Tiles512_025/train/' 
            test_dir = workdir + tissue +'/Tiles512_025/'+ test +'/' 
            out_file = workdir+tissue+'/temp_scaling/hybridtrain'+ test +'_512_tile_temp_scaled'+ tissue +'.csv'
            model_path = workdir+tissue+'/selected_tile_ACC/hybrid_tile_512_BreastMammaryTissue_trained_model_vgg19_bn_1024RSM.pt'

        data = {
            'train': ImageFolderWithPaths(train_dir, transform=transforms_dict['train']),
            'test': ImageFolderWithPaths(test_dir, transform=transforms_dict['test']),
        }

        dataloaders_eval = {
            'train': DataLoader(data['train'], batch_size=bs, shuffle=False),
            'test': DataLoader(data['test'], batch_size=bs, shuffle=False),
        }

        model_with_temp.eval()
        true_labels = np.empty((0), int)
        preds_npy = np.empty((0), int)
        files_path = np.empty((0))
        probs_npy = np.empty((0, 2))

        for phase in ['train', 'test']:
            for inputs, labels, paths in dataloaders_eval[phase]:
                inputs = inputs.to(device)
                labels = labels.to(device)
                with torch.no_grad():
                    logits = model_with_temp(inputs)# apply temperature scaling here
                    probs = F.softmax(logits, dim=1)  # calibrated probabilities
                    scores, preds = torch.max(probs, 1)

                    true_labels = np.append(true_labels, labels.cpu().numpy(), axis=0)
                    preds_npy = np.append(preds_npy, preds.cpu().numpy(), axis=0)
                    probs_npy = np.append(probs_npy, probs.cpu().numpy(), axis=0)
                    files_path = np.append(files_path, np.array(paths), axis=0)

        predictions_df = pd.DataFrame({
            'file': files_path,
            'true_label': true_labels,
            'prediction': preds_npy,
            'prob_class0': probs_npy[:, 0],
            'prob_class1': probs_npy[:, 1],
        })
        predictions_df['tile'] = predictions_df['file'].apply(lambda x: os.path.basename(x))
        predictions_df.to_csv(out_file, index=False)
        print(f"Saved calibrated predictions to {out_file}")
