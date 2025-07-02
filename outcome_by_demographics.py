import pandas as pd
import numpy as np
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import StandardScaler, LabelEncoder
from sklearn.metrics import classification_report, confusion_matrix
from tensorflow.keras.models import Sequential
from tensorflow.keras.layers import Dense, Dropout
from tensorflow.keras.utils import to_categorical
import shap

df = pd.read_csv('XXXXXX.csv')

df = df.drop(columns=[''])

df['sex'] = df['sex'].map({'M': 0, 'F': 1})
df['handedness'] = df['handedness'].map({'L': 0, 'R': 1, 'A': 2})
df['seizure freedom'] = df['seizure freedom'].map({'ILAE-1': 1, 'ILAE-2': 2', ILAE-3': 3,'ILAE-4': 4})
df = df.dropna()

X = df[['age', 'sex', 'handedness', 'epilepsy age of onset', 'epilepsy duration', 'percentage resection volume']]
y = df['seizure freedom']

X = df[['age', 'sex', 'handedness', 'epilepsy age of onset',
        'epilepsy duration', 'percentage resection volume']]
y = df['seizure freedom']

le = LabelEncoder()
y = le.fit_transform(y)
num_classes = len(np.unique(y))
y_cat = to_categorical(y, num_classes)

X_train, X_test, y_train, y_test = train_test_split(X, y_cat, test_size=0.2, random_state=42, stratify=y)

scaler = StandardScaler()
X_train_scaled = scaler.fit_transform(X_train)
X_test_scaled = scaler.transform(X_test)

# don't use these
# clf_rf = RandomForestClassifier(n_estimators=100, random_state=42)
# clf_svc = SVC(kernel='rbf', probability=True, random_state=42)
# clf_lr = LogisticRegression(max_iter=1000, random_state=42)


# training
model = Sequential()
model.add(Dense(64, activation='relu', input_dim=X_train_scaled.shape[1]))
model.add(Dropout(0.3))
model.add(Dense(32, activation='relu'))
model.add(Dense(num_classes, activation='softmax'))

model.compile(optimizer='adam', loss='categorical_crossentropy', metrics=['accuracy'])
model.fit(X_train_scaled, y_train, epochs=50, batch_size=16, validation_split=0.1, verbose=1)

y_pred_ann = model.predict(X_test_scaled)
y_pred_classes = np.argmax(y_pred_ann, axis=1)
y_test_classes = np.argmax(y_test, axis=1)

print(confusion_matrix(y_test_classes, y_pred_classes))
print(classification_report(y_test_classes, y_pred_classes, target_names=le.classes_))

explainer = shap.DeepExplainer(model, X_train_scaled[:100])  # use a subset to speed up
shap_values = explainer.shap_values(X_test_scaled[:50])

# summary for each class
for i in range(num_classes):
    shap.summary_plot(shap_values[i], X_test_scaled[:50], feature_names=X.columns, show=False)
