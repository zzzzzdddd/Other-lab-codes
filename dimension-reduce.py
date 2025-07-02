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
from sklearn.manifold import TSNE
from sklearn.metrics import accuracy_score
from sklearn.svm import SVC
import joblib

lesionls = []
hcls = []
tr_split = 0.8
directory_path = 'C:/Users/zofdi/Downloads/NewMRF/NewMRF'
all_items = os.listdir(directory_path)
files = [os.path.join(directory_path, item) for item in all_items]
id_split = int(floor(len(files)*tr_split))
trainDS = files[0:id_split]
testDS = files[id_split:]

# for idx in range(0,len(trainDS)):
for idx in range(0,len(trainDS)):
  n1_img = nib.load(trainDS[idx])
  n1 = n1_img.get_fdata()
  n1[:,:,:,-5] = (n1[:,:,:,-5]>=3)&(n1[:,:,:,-5]<=8)*1
  n1[:,:,:,-4] = (n1[:,:,:,-4]>=3)&(n1[:,:,:,-4]<=8)*1
  n1[:,:,:,-3] = (n1[:,:,:,-3]>=3)&(n1[:,:,:,-3]<=8)*1
  n1[:,:,:,-2] = (n1[:,:,:,-2]>=3)&(n1[:,:,:,-2]<=8)*1
  n1[:,:,:,-1] = n1[:,:,:,-1]*(np.sum(n1[:,:,:,-5:-1], axis = -1) > 0)

  image = np.array(n1[12:172,14:206,12:172, :])
  mask = np.array(n1[12:172,14:206,12:172,-1:]).squeeze()
  lesionls.append(image[(mask==1),:])
  if center_of_mass(mask)[1]>image.shape[1]/2:
    hcls.append(image[:,image.shape[1]//2:,:,:][(mask[:,image.shape[1]//2:,:]==0),:])
  else:
    hcls.append(image[:,:image.shape[1]//2,:,:][(mask[:,:image.shape[1]//2,:]==0),:])
lesion_stack = np.vstack(lesionls)
hc_stack = np.vstack(hcls)

n_r_keep = hc_stack.shape[0] //10
selected_indices = np.random.choice(hc_stack.shape[0], n_r_keep, replace=False)
hc_stack = hc_stack[selected_indices, :]

n_r_keep = lesion_stack.shape[0]
selected_indices = np.random.choice(lesion_stack.shape[0], n_r_keep, replace=False)
lesion_stack = lesion_stack[selected_indices, :]

hc_stack = hc_stack[(hc_stack[:, 0] > 0)]
hc_stack = hc_stack[hc_stack[:, 5] < 0.2]
Xy = np.vstack((lesion_stack, hc_stack))
np.random.shuffle(Xy)
# X = Xy[:,-5:-1]
X = np.hstack([Xy[:, 1:6], Xy[:, -5:-1]])
y = Xy[:,-1]
X[:, 0] = X[:, 0] / 3000
X[:, 1] = X[:, 1] / 2000
X = np.nan_to_num(X)
mapper = umap.UMAP(n_neighbors=20, n_components=3).fit(X, y)
embedding_3d = mapper.transform(X)

# joblib.dump(mapper, 'C:/mapper')
# joblib.dump(embedding_3d, 'C:/train_embedding3d.pkl')
# joblib.dump(y, 'C:/train_label3d.pkl')


fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
scatter = ax.scatter(embedding_3d[:, 0], embedding_3d[:, 1], embedding_3d[:, 2], c=y, cmap='Spectral', alpha = 0.2)
legend = ax.legend(*scatter.legend_elements(), loc="best", title="Classes")
ax.add_artist(legend)
plt.title('3D UMAP Embedding')
plt.show()

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
X[:, 0] = X[:, 0] / 3000
X[:, 1] = X[:, 1] / 2000
X = np.nan_to_num(X)
test_embedding = mapper.transform(X)
# joblib.dump(test_embedding, 'C:/test_embedding3d.pkl')
# joblib.dump(y, 'C:/test_label3d.pkl')

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
scatter = ax.scatter(test_embedding[:, 0], test_embedding[:, 1], test_embedding[:, 2], c=y, cmap='Spectral', alpha = 0.2)
# legend = ax.legend(*scatter.legend_elements(), loc="best", title="Classes")
# ax.add_artist(legend)
plt.title('3D UMAP Embedding test set')
plt.show()

