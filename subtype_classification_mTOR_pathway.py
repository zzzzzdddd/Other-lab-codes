import pandas as pd
import numpy as np
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import StandardScaler, LabelEncoder
from sklearn.metrics import classification_report, confusion_matrix
from tensorflow.keras.models import Sequential
from tensorflow.keras.layers import Dense, Dropout
from tensorflow.keras.utils import to_categorical
import shap

df = pd.read_csv('XXXXXXXXX.csv')

df['gender'] = df['gender'].map({'M': 0, 'F': 1})
df = df.dropna()

# feature testing
X = df[['AKT3', 'PIK3CA', 'MTOR', 'RHEB', 'TSC1', 'TSC2',
        'DEPDC5', 'NPRL2', 'NPRL3', 'age', 'gender']]
y = df['subtype']

if 'MRF_T1' in df.columns and 'MRF_T2' in df.columns:
    X_tested = pd.concat([X_base, df[['MRF_T1', 'MRF_T2']]], axis=1)
else:
    X_tested = X_base

le = LabelEncoder()
y_encoded = le.fit_transform(y)
num_classes = len(np.unique(y_encoded))
y_cat = to_categorical(y_encoded, num_classes)
X_train, X_test, y_train, y_test = train_test_split(X, y_cat, test_size=0.2, random_state=42, stratify=y_encoded)


scaler = StandardScaler()
X_train_scaled = scaler.fit_transform(X_train)
X_test_scaled = scaler.transform(X_test)


# clf_rf = RandomForestClassifier(n_estimators=100, random_state=42)
# clf_rf.fit(X_train_scaled, np.argmax(y_train, axis=1))
# y_pred_rf = clf_rf.predict(X_test_scaled)

# clf_svc = SVC(kernel='rbf', probability=True, random_state=42)
# clf_svc.fit(X_train_scaled, np.argmax(y_train, axis=1))
# y_pred_svc = clf_svc.predict(X_test_scaled)

# clf_lr = LogisticRegression(max_iter=1000, random_state=42)
# clf_lr.fit(X_train_scaled, np.argmax(y_train, axis=1))
# y_pred_lr = clf_lr.predict(X_test_scaled)


model = Sequential()
model.add(Dense(64, activation='relu', input_dim=X_train_scaled.shape[1]))
model.add(Dropout(0.3))
model.add(Dense(32, activation='relu'))
model.add(Dense(num_classes, activation='softmax'))

model.compile(optimizer='adam', loss='categorical_crossentropy', metrics=['accuracy'])
model.fit(X_train_scaled, y_train, epochs=300, batch_size=16, validation_split=0.1, verbose=1)

y_pred_ann = model.predict(X_test_scaled)
y_pred_classes = np.argmax(y_pred_ann, axis=1)
y_test_classes = np.argmax(y_test, axis=1)

print(confusion_matrix(y_test_classes, y_pred_classes))
print(classification_report(y_test_classes, y_pred_classes, target_names=le.classes_))

# Feature Importance
explainer = shap.DeepExplainer(model, X_train_scaled[:100])  # subset to speed up
shap_values = explainer.shap_values(X_test_scaled[:50])

for i in range(num_classes):
    shap.summary_plot(shap_values[i], X_test_scaled[:50], feature_names=X.columns, show=False)
