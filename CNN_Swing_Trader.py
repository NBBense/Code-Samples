import yfinance as yf
import pandas as pd
import numpy as np
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import StandardScaler
import tensorflow as tf
from tensorflow.keras.models import Sequential
from tensorflow.keras.layers import Conv1D, MaxPooling1D, Flatten, Dense, Dropout
import matplotlib.pyplot as plt

# Fetch historical data for a list of stocks
def fetch_stock_data(tickers, start_date, end_date):
    data = {}
    for ticker in tickers:
        stock_data = yf.download(ticker, start=start_date, end=end_date)
        data[ticker] = stock_data
    return data

# Load a list of stock symbols from a CSV file
stock_csv = 'stock_symbols.csv'
df = pd.read_csv(stock_csv)
stocks = df['symbol'].tolist()

# Define time period
start_date = '2022-01-01'
end_date = '2023-01-01'

data = fetch_stock_data(stocks, start_date, end_date)

def create_features(data, window_size):
    X, y = [], []
    for ticker, df in data.items():
        df = df[['Open', 'High', 'Low', 'Close', 'Volume']]
        df = df.pct_change().fillna(0)  # Normalize features with percentage change
        
        for i in range(len(df) - window_size):
            X.append(df.iloc[i:i + window_size].values)
            # Define target: predicting if the stock price will go up in the next day
            y.append(1 if df['Close'].iloc[i + window_size] > df['Close'].iloc[i + window_size - 1] else 0)
    
    return np.array(X), np.array(y)

window_size = 20
X, y = create_features(data, window_size)

# Split the data
X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.2, random_state=42)

# Normalize the data
scaler = StandardScaler()
X_train = scaler.fit_transform(X_train.reshape(-1, X_train.shape[-1])).reshape(X_train.shape)
X_test = scaler.transform(X_test.reshape(-1, X_test.shape[-1])).reshape(X_test.shape)

# Define the CNN model
model = Sequential([
    Conv1D(filters=64, kernel_size=3, activation='relu', input_shape=(window_size, 5)),
    MaxPooling1D(pool_size=2),
    Conv1D(filters=128, kernel_size=3, activation='relu'),
    MaxPooling1D(pool_size=2),
    Flatten(),
    Dense(64, activation='relu'),
    Dropout(0.5),
    Dense(1, activation='sigmoid')
])

model.compile(optimizer='adam', loss='binary_crossentropy', metrics=['accuracy'])

# Train the model
history = model.fit(X_train, y_train, epochs=10, batch_size=32, validation_split=0.1)

# Evaluate the model
loss, accuracy = model.evaluate(X_test, y_test)
print(f'Test Accuracy: {accuracy:.2f}')

# Plot training history
plt.plot(history.history['accuracy'], label='Accuracy')
plt.plot(history.history['val_accuracy'], label = 'Validation Accuracy')
plt.xlabel('Epoch')
plt.ylabel('Accuracy')
plt.legend()
plt.show()

# Make predictions
predictions = model.predict(X_test)

# Convert predictions to binary class labels
predicted_classes = (predictions > 0.5).astype(int)
