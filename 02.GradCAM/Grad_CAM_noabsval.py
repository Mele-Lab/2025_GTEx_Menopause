import os
import torch
import torch.nn as nn
import torchvision.transforms as T
from torchvision.models import vgg19_bn
from PIL import Image
import numpy as np
from captum.attr import LayerGradCam, LayerAttribution
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap

# -------- PARAMETERS (edit these) --------
tile        = "GTEX-PWCY-1425_061"
#tile    = "GTEX-14DAQ-2225_38"
#tile     = "GTEX-U3ZN-1025_52"
model_path  = '/X/files/Explainability/Uterus/hybrid_tile_512_Uterus_trained_model_vgg19_bn_1024RSM.pt'
img_root    = '/X/files/Explainability/Uterus/Tiles512_05/train/young/'
img_path    = img_root + tile + ".png"
output_dir  = '/X/files/Explainability/Uterus/comparison'

os.makedirs(output_dir, exist_ok=True)

# -------- DEVICE SETUP --------
device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')

# -------- MODEL DEFINITION & LOADING --------
model = vgg19_bn(weights=None)  # no built‐in weights

# rebuild the same two‐layer head you used in training:
model.classifier[6] = nn.Sequential(
    nn.Linear(4096, 1024),
    nn.ReLU(inplace=True),
    nn.Dropout(p=0.5),
    nn.Linear(1024, 2),
)

# load your checkpoint
checkpoint = torch.load(model_path, map_location=device)
# if you saved state_dict directly:
model.load_state_dict(checkpoint)
# — or if you wrapped it in {"state_dict": ...}, do:
# model.load_state_dict(checkpoint['state_dict'])

model.to(device).eval()

# -------- IMAGE PREPROCESSING --------
mean = [0.485, 0.456, 0.406]
std  = [0.229, 0.224, 0.225]
preprocess = T.Compose([
    T.ToTensor(),
    T.Normalize(mean=mean, std=std),
])

def load_image(path):
    img = Image.open(path).convert('RGB')
    tensor = preprocess(img).unsqueeze(0).to(device)
    return img, tensor

img_path   = os.path.join(img_root, tile + ".png")
img, input_tensor = load_image(img_path)

# -------- GRAD-CAM COMPUTATION --------
# choose the last conv layer in VGG-19-BN:
# here, features[48] is the 4th conv in block5 (just before pooling)
target_layer = model.features[48]

gradcam = LayerGradCam(model, target_layer)

# target=1 for “old” class (1=young, 0=old)
mask = gradcam.attribute(input_tensor, target=1)

# upsample to original tile size
mask_upsampled = LayerAttribution.interpolate(mask, input_tensor.shape[2:])
heatmap = mask_upsampled.squeeze().cpu().detach().numpy()

# -------- SAVE & VISUALIZE --------
# save raw heatmap
np.save(os.path.join(output_dir, tile + ".npy"), heatmap)
print(f"✅ Grad-CAM heatmap (npy) saved to {output_dir}")

# -------- TRANSPARENT-TO-GREEN COLORMAP --------
transparent_green = LinearSegmentedColormap.from_list(
    'transparent_green',
    [
        (0.0, (0, 0, 0, 0)),  # fully transparent at low values
        (1.0, (0, 1, 0, 1))   # solid green at high values
    ]
)

# -------- SAVE & VISUALIZE OVERLAY --------
fig, ax = plt.subplots(figsize=(6, 6))
ax.imshow(img)
# Use custom colormap; no extra alpha needed
overlay = ax.imshow(heatmap, cmap=transparent_green)
ax.axis('off')
# colorbar to interpret intensity
cbar = fig.colorbar(overlay, ax=ax, fraction=0.046, pad=0.04)
cbar.set_label('Grad-CAM Intensity', rotation=270, labelpad=15)

out_png = os.path.join(output_dir, tile + ".png")
fig.savefig(out_png, dpi=300, bbox_inches='tight')
print(f"✅ Overlay saved to {out_png}")
plt.close(fig)
