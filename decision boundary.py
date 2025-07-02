import torch
import numpy as np
import matplotlib.pyplot as plt
from sklearn.datasets import make_moons
from torch import nn
import joblib
import torchvision


# Generate synthetic data
X = joblib.load('C:/Users/zofdi/Downloads/train_embedding.pkl')
y = joblib.load('C:/Users/zofdi/Downloads/train_label.pkl')
X = torch.from_numpy(X).float()
y = torch.from_numpy(y).long()
alpha = 0.9
gamma = 4
# # Visualizing the dataset
# plt.scatter(X[:, 0], X[:, 1], c=y, cmap='viridis')
# plt.title("Synthetic Dataset")
# plt.show()

# Neural network model
class SimpleNN(nn.Module):
    def __init__(self):
        super(SimpleNN, self).__init__()
        self.layers = nn.Sequential(
            nn.Linear(2, 4),
            nn.ReLU(),
            nn.Linear(4, 4),
            nn.ReLU(),
            nn.Linear(4, 4),
            nn.ReLU(),
            nn.Linear(4, 2),
        )

    def forward(self, x):
        return self.layers(x)

def plot_decision_boundary():
    x_min, x_max = X[:, 0].min() - 0.1, X[:, 0].max() + 0.1
    y_min, y_max = X[:, 1].min() - 0.1, X[:, 1].max() + 0.1
    xx, yy = np.meshgrid(np.linspace(x_min, x_max, 100),
                         np.linspace(y_min, y_max, 100))

    Z = model(torch.from_numpy(np.c_[xx.ravel(), yy.ravel()]).float())
    Z = Z.detach().numpy()
    Z = np.argmax(Z, axis=1).reshape(xx.shape)

    plt.contourf(xx, yy, Z, alpha=0.8, cmap='viridis')
    plt.scatter(X[:, 0], X[:, 1], c=y, alpha=0.2, cmap='viridis')
    plt.title("Decision Boundary")
    plt.show()

model = SimpleNN()

# criterion = nn.CrossEntropyLoss()
criterion = nn.CrossEntropyLoss()
optimizer = torch.optim.Adam(model.parameters(), lr=0.01)

# Training loop
for epoch in range(10):
    # Forward pass
    outputs = model(X)
    # loss = torch.nan_to_num(torchvision.ops.sigmoid_focal_loss(outputs.squeeze(), y, alpha = alpha, gamma = gamma).sum())
    loss = criterion(outputs, y)
    # Backward and optimize
    optimizer.zero_grad()
    loss.backward()
    optimizer.step()

    if (epoch+1) % 100 == 0:
        print(f'Epoch [{epoch+1}/1000], Loss: {loss.item():.4f}')
    plot_decision_boundary()