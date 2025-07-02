import nibabel as nib
import numpy as np
import pandas as pd
import os
import glob
import torch
import torch.nn as nn
import torch.nn.functional as F
from torch.optim.lr_scheduler import StepLR
import torch.optim as optim
import random
import time
from torch.utils.data import DataLoader, TensorDataset, Dataset
from tqdm import tqdm
import matplotlib.pyplot as plt
from math import floor
import scipy
from scipy.ndimage import center_of_mass
from sklearn.model_selection import train_test_split
from sklearn.metrics import accuracy_score
from functools import partial
import torch.nn.functional as F
import torchvision
import os
import umap
# from umap.parametric_umap import load_ParametricUMAP
import pickle
import joblib
lesionls = []
hcls = []
tr_split = 0.8
directory_path = 'XXXXXXXX/NewMRF'
all_items = os.listdir(directory_path)
files = [os.path.join(directory_path, item) for item in all_items]
id_split = int(floor(len(files)*tr_split))
trainDS = files[0:id_split]
testDS = files[id_split:]

for idx in range(0,len(testDS)):
  n1_img = nib.load(testDS[idx])
  n1 = n1_img.get_fdata()
  n1[:,:,:,-5] = (n1[:,:,:,-5]>=3)&(n1[:,:,:,-5]<=8)*1
  n1[:,:,:,-4] = (n1[:,:,:,-4]>=3)&(n1[:,:,:,-4]<=8)*1
  n1[:,:,:,-3] = (n1[:,:,:,-3]>=3)&(n1[:,:,:,-3]<=8)*1
  n1[:,:,:,-2] = (n1[:,:,:,-2]>=3)&(n1[:,:,:,-2]<=8)*1
  n1[:,:,:,-1] = n1[:,:,:,-1]*(np.sum(n1[:,:,:,-5:-1], axis = -1) > 0)

  image = np.array(n1[12:172,14:206,12:172, :])
  mask = np.array(n1[12:172,14:206,12:172,-1:]).squeeze()
  lesionls.append(image[(mask>=0),:])
hc_stack = np.vstack(lesionls)

hc_stack = hc_stack[(hc_stack[:, 0] > 0)]
hc_stack = hc_stack[hc_stack[:, 5] < 0.2]
Xy = hc_stack
X = np.hstack([Xy[:, 1:6], Xy[:, -5:-1]])
y = Xy[:,-1]
X = np.nan_to_num(X)
# joblib.dump(y, 'C:/Users/zofdi/Downloads/test_label.pkl')
mapper = umap.UMAP(n_neighbors=20, n_components=2).fit(X, y)
test_embedding = mapper.transform(X)
joblib.dump(test_embedding, 'XXXXXXXtest_embedding.pkl')
