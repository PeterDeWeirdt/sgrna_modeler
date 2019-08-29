from tensorflow import keras as k
import xgboost as xgb
import sklearn as sk

#xgb.XGBRegressor().fit()

#k.wrappers.scikit_learn.KerasRegressor().fit()
def build_kim2018(input_shape=(34, 4)):
    Input_SEQ = k.layers.Input(shape=input_shape)
    C1 = k.layers.Convolution1D(80, 5, activation='relu')(Input_SEQ)
    P1 = k.layers.AveragePooling1D(2)(C1)
    F = k.layers.Flatten()(P1)
    DO1 = k.layers.Dropout(0.3)(F)
    D1 = k.layers.Dense(80, activation='relu')(DO1)
    DO2 = k.layers.Dropout(0.3)(D1)
    D2 = k.layers.Dense(40, activation='relu')(DO2)
    DO3 = k.layers.Dropout(0.3)(D2)
    D3 = k.layers.Dense(40, activation='relu')(DO3)
    DO4 = k.layers.Dropout(0.3)(D3)
    Output = k.layers.Dense(1, activation='linear')(DO4)
    model = k.models.Model(inputs = Input_SEQ, outputs = Output)
    return model

def build_kim2019(input_shape=(30,4)):
    pass
