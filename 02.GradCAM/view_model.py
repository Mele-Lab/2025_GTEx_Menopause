import torch

# This is the same .pt you use for inference
CKPT_FILE = "/X/Explainability/Ovary/hybrid_tile_512_Ovary_trained_model_vgg19_bn_1024RSM.pt"

# Load the saved state_dict (no shards, no persistent_load issues)
state_dict = torch.load(CKPT_FILE, map_location="cpu")

# Print out every parameter and its shape
for name, tensor in state_dict.items():
    print(f"{name:40s} → {tuple(tensor.shape)}")
