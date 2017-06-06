import numpy as np
import math
from keras.models import Sequential
from keras.layers import Dense
from keras.layers import LSTM
from sklearn.preprocessing import normalize
from sklearn.metrics import mean_squared_error

def generate_data(data):
    '''
    Splits the data into two parts. 
    Two thirds of them are training data.
    The other one third are testing data.
    '''
    num = data.shape[0]
    split = int(num / 3 * 2)
    train_data = data[0:split, :]
    test_data = data[split:num, :]
    return train_data, test_data
    
def shift_data(X, Y, look_back=1):
    dataX = []
    dataY = []
    num = X.shape[0]
    for i in range(num-look_back):
        x = X[i:(i+look_back), :].flatten()
        y = Y[i+look_back, :].flatten()
        dataX.append(x)
        dataY.append(y)
    return np.array(dataX), np.array(dataY)

# load dataset
X = np.load('X.npy')
Y = np.load('Y.npy')

X = X.astype('float32')
Y = Y.astype('float32')

X = normalize(X)

trainX, testX = generate_data(X)
trainY, testY = generate_data(Y)

#print('trainX.shape:', trainX.shape)
#print('trainY.shape:', trainX.shape)
#print('testX.shape:', testX.shape)
#print('testY.shape:', testX.shape)

look_back = 2
trainX, trainY = shift_data(trainX, trainY, look_back=look_back)
testX, testY = shift_data(testX, testY, look_back=look_back)
#print('trainX.shape:', trainX.shape)
#print('trainY.shape:', trainX.shape)
#print('testX.shape:', testX.shape)
#print('testY.shape:', testX.shape)

trainX = np.reshape(trainX, (trainX.shape[0], 1, trainX.shape[1]))
testX = np.reshape(testX, (testX.shape[0], 1, testX.shape[1]))
print('trainX.shape:', trainX.shape)
print('trainY.shape:', trainY.shape)
print('testX.shape:', testX.shape)
print('testY.shape:', testY.shape)

# create and fit the LSTM network
model = Sequential()
model.add(LSTM(4, input_shape=(1, look_back * 4453)))
model.add(Dense(66))
model.compile(loss='mean_squared_error', optimizer='adam')
model.fit(trainX, trainY, epochs=250, batch_size=1, verbose=2)

# make predictions
trainPredict = model.predict(trainX)
testPredict = model.predict(testX)

print('trainPredcit.shape:', trainPredict.shape)
print('testPredict.shape:', testPredict.shape)
np.save('trainPredict.npy', trainPredict)
np.save('testPredict.npy', testPredict)

