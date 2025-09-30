import os
import torch
import torch.nn as nn
import torchvision.transforms as T
from torchvision.models import vgg19_bn
from PIL import Image
import numpy as np
from captum.attr import LayerGradCam, LayerAttribution
import matplotlib.pyplot as plt
from matplotlib import cm

# -------- PARAMETERS (edit these) --------
tile        = "GTEX-1N7R6-1825_058"
#tile    = "GTEX-14DAQ-2225_38"
#tile     = "GTEX-U3ZN-1025_52"
model_path  = '/X/files/Explainability/Ovary/hybrid_tile_512_Ovary_trained_model_vgg19_bn_1024RSM.pt'
img_root    = '/X/files/Explainability/Ovary/Tiles512_05/train/young/'
img_path    = img_root + tile + ".png"
output_dir  = '/X/files/Explainability/Ovary/comparison'


low_threshold   = 0.25
# overlay alpha for visible regions
alpha_val       = 0.6

os.makedirs(output_dir, exist_ok=True)

# -------- DEVICE SETUP --------
device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')

# -------- MODEL DEFINITION & LOADING --------
model = vgg19_bn(weights=None)
model.classifier[6] = nn.Sequential(
    nn.Linear(4096, 1024),
    nn.ReLU(inplace=True),
    nn.Dropout(p=0.5),
    nn.Linear(1024, 2),
)
checkpoint = torch.load(model_path, map_location=device)
model.load_state_dict(checkpoint)
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

img_path, _ = os.path.join(img_root, tile + ".png"), None
img, input_tensor = load_image(img_path)

# -------- GRAD-CAM COMPUTATION --------
target_layer = model.features[48]
gradcam = LayerGradCam(model, target_layer)
mask = gradcam.attribute(input_tensor, target=1)
mask_up = LayerAttribution.interpolate(mask, input_tensor.shape[2:])
heatmap = np.abs(mask_up.squeeze().cpu().detach().numpy())

# -------- NORMALIZE & APPLY THRESHOLD --------
# normalize to [0,1]
heat_min, heat_max = heatmap.min(), heatmap.max()
heat_norm = (heatmap - heat_min) / (heat_max - heat_min + 1e-8)
# create alpha mask based on threshold
alpha_mask = np.where(heat_norm > low_threshold, alpha_val, 0)

# -------- CREATE RGBA HEATMAP --------
greens_cmap = cm.get_cmap('Greens')
rgba_heatmap = greens_cmap(heat_norm)
rgba_heatmap[..., 3] = alpha_mask

# -------- SAVE RAW HEATMAP --------
np.save(os.path.join(output_dir, tile + ".npy"), heatmap)
print(f"✅ Grad-CAM heatmap (npy) saved to {output_dir}")

# -------- OVERLAY & SAVE FIGURE --------
fig, ax = plt.subplots(figsize=(6, 6))
ax.imshow(img)
ax.imshow(rgba_heatmap, interpolation='nearest')
ax.axis('off')

# colorbar for visual reference
norm = plt.Normalize(vmin=0, vmax=heat_max)
sm = plt.cm.ScalarMappable(cmap='Greens', norm=norm)
sm.set_array([])
cbar = fig.colorbar(sm, ax=ax, fraction=0.046, pad=0.04)
cbar.set_label('Grad-CAM Intensity', rotation=270, labelpad=15)

out_png = os.path.join(output_dir, tile + ".png")
fig.savefig(out_png, dpi=300, bbox_inches='tight')
print(f"✅ Overlay saved to {out_png}")
plt.close(fig)
