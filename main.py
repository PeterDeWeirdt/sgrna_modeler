import sgrna_modeler.models as sg
import sgrna_modeler.datasets as da
from scipy import stats
import pandas as pd
import os

if __name__ == '__main__':
    models = [sg.Model_Doench2016(), sg.Model_Kim2018()]
    train_datum = [da.load_doench_2016_train(), da.load_meyers_2017_train(), da.load_kim_2019_train()]
    test_datum = [da.load_doench_2016_test(), da.load_meyers_2017_test(), da.load_kim_2019_test()]

    predictions = []
    for model in models:
        print(model.base_name)
        for train_data in train_datum:
            print('\t' + train_data.name)
            model.fit(train_data)
            for test_data in test_datum:
                print('\t\t' + test_data.name)
                predicted_test_data = model.predict(test_data)
                predictions.append(predicted_test_data)
    all_predictions = pd.concat(predictions)
    all_predictions.to_csv(os.path.join('test_predictions', '2019-08-29_predictions.csv'))
