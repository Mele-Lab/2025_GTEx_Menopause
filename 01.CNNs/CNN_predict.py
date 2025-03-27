
##Prediction of tiles

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

from PIL import Image
from sklearn.metrics import confusion_matrix
from torch.utils.data.sampler import SubsetRandomSampler
from torchvision import datasets, models, transforms
from sklearn.metrics import roc_curve
from sklearn.metrics import auc
from sklearn.metrics import RocCurveDisplay

seed = 180
random.seed(seed)
np.random.seed(seed)
torch.manual_seed(seed)
torch.backends.cudnn.deterministic = True
torch.backends.cudnn.benchmark = False
workdir = '/X//Laura/05.CNN/35yo_separation_FEMALE/'
tissues = ['Uterus', 'Ovary', 'Vagina', 'BreastMammaryTissue']
tests =['test', 'young_post', 'val_and_mid']

for test in tests:
    for tissue in tissues:
        if tissue == 'Uterus':
            drop=0.5 #After the fine-tunning, each tissue model has its best parameters
            train_directory = workdir + tissue + '/Tiles512_05/train/' 
            test_directory = workdir + tissue +'/Tiles512_05/'+ test +'/' 
            out_file = workdir+tissue+'/selected_tile/hybridtrain'+ test +'_512_tile_Alinsaif_age30_vgg19_predictions_filtered'+ tissue +'.csv'
            # model_path =f"{workdir}{tissue}/tile_hybrid_tile_512_0.001_0.5_0.001_{drop}_15_{tissue}_trained_model_vgg19_bn_1024RSM.pt"
            model_path = workdir+tissue+'/selected_tile/hybrid_tile_512_Uterus_trained_model_vgg19_bn_1024RSM.pt'
    
        elif tissue == 'Ovary':
            drop=0.4
            train_directory = workdir + tissue + '/Tiles512_05/train/' 
            test_directory = workdir + tissue +'/Tiles512_05/'+ test +'/' 
            out_file = workdir+tissue+'/selected_tile/hybridtrain'+ test +'_512_tile_Alinsaif_age30_vgg19_predictions_filtered'+ tissue +'.csv'
            # model_path =f"{workdir}{tissue}/tile_hybrid_tile_512_0.001_0.5_0.001_{drop}_15_{tissue}_trained_model_vgg19_bn_1024RSM.pt"
            model_path = workdir+tissue+'/selected_tile/hybrid_tile_512_Ovary_trained_model_vgg19_bn_1024RSM.pt'
    
    
        elif tissue == 'Vagina':
            drop=0.5
            train_directory = workdir + tissue + '/Tiles256_05/train/' 
            test_directory = workdir + tissue +'/Tiles256_05/'+ test +'/' 
            # model_path =f"{workdir}{tissue}/tile_hybrid_tile_256_0.001_0.5_0.001_{drop}_15_{tissue}_trained_model_vgg19_bn_1024RSM.pt"
            out_file = workdir+tissue+'/selected_tile/hybridtrain'+ test +'_256_tile_Alinsaif_age30_vgg19_predictions_filtered'+ tissue +'.csv'
            model_path = workdir+tissue+'/selected_tile/hybrid_tile_256_Vagina_trained_model_vgg19_bn_1024RSM.pt'
    
    
        elif tissue == 'BreastMammaryTissue':
            drop=0.3
            train_directory = workdir + tissue + '/Tiles512_025_filtering_extra/train/' 
            test_directory = workdir + tissue +'/Tiles512_025_filtering_extra/'+ test +'/' #For the middle age prediction, train:middle, test:validation
            out_file = workdir+tissue+'/selected_tile/hybridtrain'+ test +'_512_tile_Alinsaif_age30_vgg19_predictions_filtered'+ tissue +'.csv'
            # model_path =f"{workdir}{tissue}/tile_filtering_hybrid_tile_512_0.001_0.5_0.001_{drop}_15_{tissue}_trained_model_vgg19_bn_1024RSM.pt"    
            model_path = workdir+tissue+'/selected_tile/hybrid_tile_512_BreastMammaryTissue_trained_model_vgg19_bn_1024RSM.pt'
    
            # Define a set of transformations to augment the data
        image_transforms = { 
            'train': transforms.Compose([
                #transforms.RandomResizedCrop(size=512, scale=(0.8, 1.0)),
                #transforms.RandomRotation(degrees=15),
                #transforms.RandomHorizontalFlip(),
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
        
        num_classes=2
        
        # ImageFolderWithPaths facilitates the retrieval of location paths for the input.
        class ImageFolderWithPaths(torchvision.datasets.ImageFolder):
            """Custom dataset that includes image file paths. 
            Extends torchvision.datasets.ImageFolder
            """
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
        
        # Print the train and test set data sizes
        #print('Sizes, train: {} test: {}'.format(
        #                dataset_sizes['train'], dataset_sizes['test']))
        
        # Class distribution in train and test sets
        dataset_distributions = {
            'train': np.unique(data['train'].targets, return_counts=True)[1],
            'test': np.unique(data['test'].targets, return_counts=True)[1]
        }
        
        model = torch.load('/X//Laura/05.CNN/pretrained_vgg19_bn.pt')
        fc_inputs = model.classifier[6].in_features
        model.classifier[6] = torch.nn.Sequential(
            torch.nn.Linear(fc_inputs, 1024),
            torch.nn.ReLU(),
            torch.nn.Dropout(0.3),
            torch.nn.Linear(1024, num_classes), 
            torch.nn.LogSoftmax(dim=1) # For using NLLLoss()
            )
        model.load_state_dict(torch.load(model_path))
        model.eval()
        
        # Convert model to be used on GPU,
        device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
        model = model.to(device)
        
        dataloaders_eval = {
            'train': torch.utils.data.DataLoader(data['train'], 
                                                    batch_size=bs, 
                                                    shuffle=False),
            'test': torch.utils.data.DataLoader(data['test'],
                                                    batch_size=bs,
                                                    shuffle=False)
            }
            
        model.eval() # Set to evaluation mode"
        
        # Setting outputs to null
        true_labels = np.empty((0), int)
        preds_npy = np.empty((0), int)
        files_path = np.empty((0))
        probs_npy = np.empty((0,2))
            
        # Each epoch has a training and test phase
        for phase in ['train', 'test']:
          # Iterate over data.
            for i, (inputs, labels, paths) in enumerate(dataloaders_eval[phase]):
            # push input to device
                inputs = inputs.to(device)
                labels = labels.to(device)
            
                true_labels = np.append(true_labels, labels.data.cpu().numpy(), axis=0)
                files_path = np.append(files_path, np.array(paths), axis=0)
            
                with torch.no_grad():
                    outputs = model(inputs)
                    probs = torch.exp(outputs.data)
                    scores, preds = torch.max(outputs.data, 1)
                    preds_npy = np.append(preds_npy, preds.data.cpu().numpy(), axis=0)
                    probs_npy = np.append(probs_npy,  probs.data.cpu().numpy(), axis=0)
        
        
        predictions_df = pd.DataFrame({'file' : files_path,
                                           'true_label':true_labels, 
                                           'prediction': preds_npy,
                                           'prob_class0': probs_npy[:,0],
                                           'prob_class1': probs_npy[:,1]})
        predictions_df['tile'] = predictions_df['file'].apply(lambda x: os.path.basename(os.path.normpath(x)))
        predictions_df.to_csv(out_file)
        predictions_df.sort_values(by='tile')





   
 

