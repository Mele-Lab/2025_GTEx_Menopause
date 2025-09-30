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
#import torchvision.models as models
from torchvision import datasets, models, transforms
from torch.utils.data import Dataset, DataLoader


from sklearn.metrics import roc_curve
from sklearn.metrics import auc
from sklearn.metrics import RocCurveDisplay

seed = 180
random.seed(seed)
np.random.seed(seed)
torch.manual_seed(seed)
torch.backends.cudnn.deterministic = True
torch.backends.cudnn.benchmark = False
workdir = 'X/05.CNN/RegressionModels/'


tissues = ['BreastMammaryTissue']
# tissues = ['Uterus', 'Ovary', 'Vagina']
tests =['test', 'young_post', 'val_and_mid']
# test = 'val_and_mid'
# tissue = 'Ovary'


for tissue in tissues:
    for test in tests:
        if tissue == 'Uterus':
            drop=0.5#for uterus
            #drop=0.3 #for myometrium
            train_directory = workdir + tissue + '/Tiles512_05/train/' 
            test_directory = workdir + tissue +'/Tiles512_05/'+ test +'/' 
            out_file = workdir+tissue+'/hybridtrain'+ test +'_512_predictions'+ tissue +'.csv'
            # model_path =f"{workdir}{tissue}/tile_hybrid_tile_512_0.001_0.5_0.001_{drop}_15_{tissue}_trained_model_vgg19_bn_1024RSM.pt"
            model_path = workdir+tissue+'/hybrid_tile_512_'+tissue+'_trained_model_vgg19_bn_1024RSM.pt'
    
        elif tissue == 'Ovary':
            drop=0.4
            train_directory = workdir + tissue + '/Tiles512_05/train/' 
            test_directory = workdir + tissue +'/Tiles512_05/'+ test +'/' 
            out_file = workdir+tissue+'/hybridtrain'+ test +'_512_predictions'+ tissue +'.csv'
            # model_path =f"{workdir}{tissue}/tile_hybrid_tile_512_0.001_0.5_0.001_{drop}_15_{tissue}_trained_model_vgg19_bn_1024RSM.pt"
            model_path = workdir+tissue+'/hybrid_tile_512_'+tissue+'_trained_model_vgg19_bn_1024RSM.pt'
    
    
        elif tissue == 'Vagina':
            drop=0.5
            train_directory = workdir + tissue + '/Tiles256_05/train/' 
            test_directory = workdir + tissue +'/Tiles256_05/'+ test +'/' 
            # model_path =f"{workdir}{tissue}/tile_hybrid_tile_256_0.001_0.5_0.001_{drop}_15_{tissue}_trained_model_vgg19_bn_1024RSM.pt"
            out_file = workdir+tissue+'/hybridtrain'+ test +'_256_predictions'+ tissue +'.csv'
            model_path = workdir+tissue+'/hybrid_tile_256_'+tissue+'_trained_model_vgg19_bn_1024RSM.pt'
    
    
        elif tissue == 'BreastMammaryTissue':
            drop=0.3
            train_directory = workdir + tissue + '/Tiles512_025/train/' 
            test_directory = workdir + tissue +'/Tiles512_025/'+ test +'/' 
            out_file = workdir+tissue+'/hybridtrain'+ test +'_512_predictions'+ tissue +'.csv'
            # model_path =f"{workdir}{tissue}/tile_filtering_hybrid_tile_512_0.001_0.5_0.001_{drop}_15_{tissue}_trained_model_vgg19_bn_1024RSM.pt"    
            model_path = workdir+tissue+'/hybrid_tile_512_'+tissue+'_trained_model_vgg19_bn_1024RSM.pt'
    
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
        test_dataset = AgeRegressionDataset(csv_file=workdir+tissue+'/labels_'+ test+'.csv',img_dir=test_directory,transform=image_transforms["test"])
        
        

        # Dataloaders
        dataloaders_eval = {
            'train': DataLoader(train_dataset, batch_size=bs, shuffle= False, num_workers=10), 
            'test': DataLoader(test_dataset, batch_size=bs, shuffle=False, num_workers=10)
        }
        
        
    
        # Load pre-trained model
        model = torch.load('X/05.CNN/pretrained_vgg19_bn.pt')
        # model = torch.load('BSC/05.CNN/pretrained_vgg19_bn.pt')

        # Modify the final fully connected layer for transfer learning
        fc_inputs = model.classifier[6].in_features
        model.classifier[6] = torch.nn.Sequential(
            torch.nn.Linear(fc_inputs, 1024),
            torch.nn.ReLU(),
            torch.nn.Dropout(drop),
            torch.nn.Linear(1024, 1)  # single continuous output
        )
        
        model.load_state_dict(torch.load(model_path))
        model.eval()
        
        # Convert model to be used on GPU,
        device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
        model = model.to(device)
    
    
            
        # Each epoch has a training and test phase

        # Evaluate both train and test sets
        for phase in ['train', 'test']:
            true_labels = []
            preds_list = []
            files_path = []

            for inputs, labels, paths in dataloaders_eval[phase]:
                inputs = inputs.to(device)
                labels = labels.to(device).unsqueeze(1)

                with torch.no_grad():
                    outputs = model(inputs)
                
                # Move tensors to CPU and convert to numpy
                preds_list.extend(outputs.cpu().numpy().flatten())
                true_labels.extend(labels.cpu().numpy().flatten())
                files_path.extend(paths)

            # Save predictions per phase
            predictions_df = pd.DataFrame({
                'file': files_path,
                'true_label': true_labels,
                'prediction': preds_list
            })
            predictions_df['tile'] = predictions_df['file'].apply(lambda x: os.path.basename(os.path.normpath(x)))

            # Save to a CSV per phase
            predictions_df.to_csv(out_file, index=False)


