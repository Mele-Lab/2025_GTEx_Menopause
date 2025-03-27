
#LIME code to interpret CNNs classification

import os
import torch
import numpy as np
from torchvision import models, transforms
from PIL import Image
from PIL import ImageEnhance
from lime import lime_image
from skimage.segmentation import mark_boundaries
from matplotlib import cm
from matplotlib.colors import ListedColormap

workdir='~/X/Laura/05.CNN/'
tissue = 'Vagina'
num_classes = 2
class_labels = ["0", "1"]

#Load pretrained VGG19 model
model = torch.load(workdir+'_pretrained_vgg19_bn.pt')

fc_inputs = model.classifier[6].in_features
model.classifier[6] = torch.nn.Sequential(
    torch.nn.Linear(fc_inputs, 1024),
    torch.nn.ReLU(),
    torch.nn.Dropout(0.4),
    torch.nn.Linear(1024, num_classes),
    torch.nn.LogSoftmax(dim=1) # For using NLLLoss()
    )
#Load out tissue model
model.load_state_dict(torch.load('~/X/Laura/05.CNN/35yo_separation_FEMALE/Vagina/ACCtile_hybrid_tile_256_0.001_0.4_0.001_0.4_15_Vagina_trained_model_vgg19_bn_1024RSM.pt'))
model.eval()


# Preprocess input image
transform = transforms.Compose([
    transforms.Resize((224, 224)),
    transforms.ToTensor(),
    transforms.Normalize(mean=[0.485, 0.456, 0.406], std=[0.229, 0.224, 0.225])
])


# Load and preprocess an image for interpretation
image_path = "~/X/Laura/05.CNN/35yo_separation_FEMALE/Vagina/Tiles512_05/validationset/val/GTEX-14PHW-2325_066.png"
image = Image.open(image_path).convert('RGB')
preprocessed_image = transform(image)
print(preprocessed_image.size())

# Convert the image to a tensor
preprocessed_image = preprocessed_image.double()
print(preprocessed_image.size())

def predict(image_array):
    image_tensor = torch.from_numpy(image_array).permute(0, 3, 1, 2)
    print("image_tensor_size", image_tensor.size())
    with torch.no_grad():
        outputs = model(image_tensor.float())
        probabilities = torch.nn.functional.softmax(outputs, dim=1)
        scores, predictions = torch.max(outputs.data, 1)
    print("type_probabilities", probabilities.dtype)
    print(probabilities.cpu().numpy())

    return probabilities.cpu().numpy()

# Create an explainer using LIME
explainer = lime_image.LimeImageExplainer()

print(len(class_labels))
print(class_labels)
print("image_input_to_explainer", preprocessed_image[0].numpy().astype('double'))

# Explain the image classification using LIME
explanation = explainer.explain_instance(
    image=preprocessed_image[0].numpy().astype('double'),
    classifier_fn=predict,
    top_labels=2,
    hide_color=0,
    num_samples=200)
print(explanation)


# Get the LIME explanation image
temp, mask = explanation.get_image_and_mask(
    explanation.top_labels[0],
    positive_only=False,
    num_features=15,
    hide_rest=False)

temp = np.squeeze(temp)

# Adjust the normalization and data type
temp = (temp / np.max(np.abs(temp))) * 255
temp = temp.astype(np.uint8)

# Convert the explanation image to PIL format
explanation_image = Image.fromarray(temp).convert('RGBA')

# Save the explanation image to a file
output_path = "~/X/Laura/06.LIME/Vagina/old/explained_back.png"
explanation_image.save(output_path)

output_path2 = "~/X/Laura/06.LIME/Vagina/old/explained2_only_mask.png"
img_boundry = mark_boundaries(preprocessed_image[0], mask, color=(1, 0, 0))
img_boundry = Image.fromarray((img_boundry * 255).astype(np.uint8))
img_boundry.save(output_path2)

print(f"Mask image saved at {output_path2}")

# Select the same class explained on the figures above.
ind = explanation.top_labels[0]
# Map each explanation weight to the corresponding superpixel
dict_heatmap = dict(explanation.local_exp[ind])
heatmap = np.vectorize(dict_heatmap.get)(explanation.segments)

# Normalize the heatmap
heatmap_normalized = (heatmap - heatmap.min()) / (heatmap.max() - heatmap.min())


# Apply a colormap to the normalized heatmap
cmap = cm.get_cmap('coolwarm')  # You can choose any colormap you prefer
heatmap_colored = cmap(heatmap_normalized)

# Convert the colored heatmap to PIL format
heatmap_image = Image.fromarray((heatmap_colored * 255).astype(np.uint8), mode='RGBA')

# Save the heatmap image to a file
output_path3 = "~/X/Laura/LIME/Vagina/old/heatmap.png"
heatmap_image.save(output_path3)

print(f"Heatmap image saved at {output_path3}")


