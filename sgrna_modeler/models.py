from sgrna_modeler import features as fe
from sklearn.model_selection import train_test_split
from sklearn import ensemble
from tensorflow import keras as k
import pandas as pd
import os
from joblib import load

def curr_path():
    return os.path.dirname(__file__)

def get_deepcpf1_weights():
    path = os.path.join(curr_path(), 'data/saved_models/Seq_deepCpf1_weights_tf.h5')
    return path

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

class Keras_sgrna_Model(object):
    def __init__(self, random_state = 7, val_frac = 0.1, base_arc = None):
        self.base_name = 'M_Kim_2018'
        self.val_frac = val_frac
        self.random_state = random_state
        if base_arc is None:
            self.base_arc = build_kim2018
        else:
            self.base_arc = base_arc
        self.train_dataset = None
        self.enzyme = None
        self.model = None
        self.model_history = None
        self.train_name = None

    def load_weights(self, weights = None, name = None):
        model = self.base_arc()
        if weights is None:
            deepcpf1_weights = get_deepcpf1_weights()
            model.load_weights(deepcpf1_weights)
            self.train_name = 'Seq-DeepCpf1'
        else:
            model.load_weights(weights)
            self.train_name = name
        self.model = model
        return self

    def fit(self, train_dataset):
        self.train_dataset = train_dataset
        self.train_name = train_dataset.name
        self.enzyme = train_dataset.enzyme
        train_val_x, y = train_dataset.get_xy()
        encoded_train_val_x = fe.encode_seqs(train_val_x)
        train_x, val_x, train_y, val_y = train_test_split(encoded_train_val_x, y, test_size=self.val_frac,
                                                          random_state=self.random_state)
        model = self.base_arc(input_shape = (self.enzyme['context_length'],4))
        model.compile(optimizer='RMSprop',loss='mse',metrics=['mae'])
        self.model_history = model.fit(train_x, train_y, epochs = 200,
                                       validation_data = (val_x, val_y),
                                       callbacks = [k.callbacks.EarlyStopping(patience=20,restore_best_weights=True),
                                                    k.callbacks.History()],
                                       verbose = 0)
        self.model = model
        return self

    def predict(self, test_dataset):
        x, y = test_dataset.get_xy()
        encoded_x = fe.encode_seqs(x)
        predictions = self.model.predict(encoded_x)
        out_data = pd.DataFrame({'kmer': x, 'y': y})
        if test_dataset.group_column:
            out_data['group'] = test_dataset.data[test_dataset.group_column]
        else:
            out_data['group'] = ''
        out_data['prediction'] = predictions
        out_data['model'] = self.base_name
        out_data['training_data'] = self.train_name
        out_data['test_data'] = test_dataset.name
        return out_data

    def predict_seqs(self, seqs):
        featurized_x = fe.encode_seqs(seqs)
        predictions = self.model.predict(featurized_x).flatten()
        return predictions

def get_rs2():
    #path = os.path.join(curr_path(), 'data/saved_models/')
    #return path
    pass

def get_enPAM_GB():
    path = os.path.join(curr_path(), 'data/saved_models/enPAM_GB.joblib')
    return path

class SKLearn_GB_sgrna_Model(object):
    """SKLearn gradient boosting for modeling sgRNA activity

    Model the activity of guides using gradient boosting, as in:
        Doench, John G., et al. "Optimized sgRNA design to maximize activity and minimize off-target \
        effects of CRISPR-Cas9." Nature biotechnology 34.2 (2016): 184.

    Parameters
    ----------
    val_frac : float, optional (default = 0.1)
        fraction of data to use as a validation set for stopping. If set to None, then no stopping will occur
    """

    def __init__(self, val_frac = 0.1, model = None, features = None):
        self.base_name = 'M_Doench_2016'
        self.val_frac = val_frac
        if model is None:
            # Gradient boosted model
            self.model = ensemble.GradientBoostingRegressor(n_iter_no_change=20,
                                                            validation_fraction = self.val_frac)
        else:
            self.model = model
        if features is None:
            # Default features for RuleSet2
            self.features = ['Pos. Ind. 1mer', 'Pos. Ind. 2mer', 'Pos. Dep. 1mer', 'Pos. Dep. 2mer', 'GC content', 'Tm']
        else:
            self.features = features
        self.enzyme = None
        self.train_dataset = None
        self.train_name = None


    def load_model(self, model, enzyme, name):
        self.enzyme = enzyme
        self.model = load(model)
        self.train_name = name
        return self

    def fit(self, train_dataset):
        self.train_name = train_dataset.name
        self.enzyme = train_dataset.enzyme
        train_val_x, y = train_dataset.get_xy()
        featurized_train_val_x = fe.featurize_guides(train_val_x, features=self.features,
                                                     guide_start = self.enzyme['guide_start'],
                                                     guide_length = self.enzyme['guide_length'])
        self.model.fit(featurized_train_val_x, y)
        return self

    def predict(self, test_dataset):
        x, y = test_dataset.get_xy()
        featurized_x = fe.featurize_guides(x, features=self.features,
                                           guide_start=test_dataset.enzyme['guide_start'],
                                           guide_length=test_dataset.enzyme['guide_length'])
        predictions = self.model.predict(featurized_x)
        out_data = pd.DataFrame({'kmer': x, 'y': y})
        if test_dataset.group_column:
            out_data['group'] = test_dataset.data[test_dataset.group_column]
        else:
            out_data['group'] = ''
        out_data['prediction'] = predictions
        out_data['model'] = self.base_name
        out_data['training_data'] = self.train_name
        out_data['test_data'] = test_dataset.name
        return out_data

    def predict_seqs(self, seqs):
        featurized_x = fe.featurize_guides(seqs, features=self.features,
                                           guide_start=self.enzyme['guide_start'],
                                           guide_length=self.enzyme['guide_length'])
        predictions = self.model.predict(featurized_x)
        return predictions

