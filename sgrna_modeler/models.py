from sgrna_modeler import architectures as ar
from sgrna_modeler import features as fe
from sklearn.model_selection import train_test_split
from sklearn import ensemble
from tensorflow import keras as k
import pandas as pd

class Model_Kim2018(object):
    def __init__(self, random_state = 7, val_frac = 0.1):
        self.base_name = 'M_Kim_2018'
        self.val_frac = val_frac
        self.random_state = random_state
        self.base_arc = ar.build_kim2018

    def fit(self, train_dataset):
        self.train_datset = train_dataset
        train_val_x, y = train_dataset.get_xy()
        encoded_train_val_x = fe.encode_seqs(train_val_x)
        train_x, val_x, train_y, val_y = train_test_split(encoded_train_val_x, y, test_size=self.val_frac,
                                                          random_state=self.random_state)
        model = self.base_arc(input_shape = (train_dataset.enzyme.context_length,4))
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
        out_data['training_data'] = self.train_datset.name
        out_data['test_data'] = test_dataset.name
        return out_data

class Model_Doench2016(object):
    def __init__(self, random_state = 7, val_frac = 0.1):
        self.base_name = 'M_Doench_2016'
        self.random_state = random_state
        self.val_frac = val_frac
        self.model = ensemble.GradientBoostingRegressor(n_iter_no_change=20,
                                                        validation_fraction = self.val_frac)
        self.features = ['Pos. Ind. 1mer', 'Pos. Ind. 2mer',
                         'Pos. Dep. 1mer', 'Pos. Dep. 2mer',
                         'GC content', 'Tm']

    def fit(self, train_dataset):
        self.train_datset = train_dataset
        train_val_x, y = train_dataset.get_xy()
        featurized_train_val_x = fe.featurize_guides(train_val_x, features=self.features,
                                                     guide_start = train_dataset.enzyme.guide_start,
                                                     guide_length = train_dataset.enzyme.guide_length)
        self.model.fit(featurized_train_val_x, y)
        return self

    def predict(self, test_dataset):
        x, y = test_dataset.get_xy()
        featurized_x = fe.featurize_guides(x, features=self.features,
                                           guide_start=test_dataset.enzyme.guide_start,
                                           guide_length=test_dataset.enzyme.guide_length)
        predictions = self.model.predict(featurized_x)
        out_data = pd.DataFrame({'kmer': x, 'y': y})
        if test_dataset.group_column:
            out_data['group'] = test_dataset.data[test_dataset.group_column]
        else:
            out_data['group'] = ''
        out_data['prediction'] = predictions
        out_data['model'] = self.base_name
        out_data['training_data'] = self.train_datset.name
        out_data['test_data'] = test_dataset.name
        return out_data

