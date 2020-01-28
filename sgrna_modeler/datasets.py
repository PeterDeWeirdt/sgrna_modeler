"""Training and testing datasets for modeling guide activity.

Includes class for generating new dataset objects, so new data can be easily modeled and tested

    :Example:

    >>> import pandas as pd
    >>> import sgrna_modeler.enzyme as en
    >>> new_data = pd.read_csv('new dataset')
    >>> new_dataset = Activity_Data(new_data, en.cas9, '30mer', 'activity', 'new data')
"""
import pandas as pd
from sgrna_modeler import enzymes as en
import os


def curr_path():
    return os.path.dirname(__file__)


class ActivityData(object):
    """Store information about activity data


    :param data: data to model
    :type data: pandas dataframe
    :param enzyme: cas9 or cas12a
    :type enzyme: dict
    :param kmer_column: sequences to model
    :type kmer_column: str
    :param name: name of the dataset
    :type name: str
    """

    def __init__(self, data, enzyme, kmer_column, activity_column, name):
        """Inits Activity data"""
        self.data = data
        self.enzyme = enzyme
        self.kmer_column = kmer_column
        self.activity_column = activity_column
        self.name = 'D_' + name

    def get_xy(self):
        """Gets modeling matrix (x) and output matrix (y)

        :return two series, x and y
        :rtype pandas series
        """
        x = self.data[self.kmer_column]
        y = self.data[self.activity_column]
        return x, y


# SpCas9 Datasets

def load_doench_2016():
    """Data from """
    data = pd.read_csv(os.path.join(curr_path(), 'data/datasets/Doench_2016.csv.zip'))
    data_class = ActivityData(data=data, enzyme=en.cas9, kmer_column='30mer',
                              activity_column='score_drug_gene_rank',
                              name='Doench_2016',
                              group_column='Target gene')
    return data_class


def load_meyers_2017_train():
    data = pd.read_csv(os.path.join(curr_path(), 'data/datasets/Meyers_2017_Train.csv.zip'))
    data_class = ActivityData(data=data, enzyme=en.cas9, kmer_column='sgRNA context sequence',
                              activity_column='mean_activity',
                              name='Meyers_2017_Train',
                              group_column='Gene Symbol')
    return data_class


def load_meyers_2017_test():
    data = pd.read_csv(os.path.join(curr_path(), 'data/datasets/Meyers_2017_Test.csv.zip'))
    data_class = ActivityData(data=data, enzyme=en.cas9, kmer_column='sgRNA context sequence',
                              activity_column='mean_activity',
                              name='Meyers_2017_Test',
                              group_column='Gene Symbol')
    return data_class


def load_kim_2019_train():
    data = pd.read_csv(os.path.join(curr_path(), 'data/datasets/Kim_2019_Train.csv.zip'))
    data_class = ActivityData(data=data, enzyme=en.cas9,
                              kmer_column='Target context sequence (4+20+3+3)',
                              activity_column='Background subtracted indel (%)',
                              name='Kim_2019_Train')
    return data_class


def load_kim_2019_test():
    data = pd.read_csv(os.path.join(curr_path(), 'data/datasets/Kim_2019_Test.csv.zip'))
    data_class = ActivityData(data=data, enzyme=en.cas9,
                              kmer_column='Target context sequence (4+20+3+3)',
                              activity_column='Background subtracted indel frequencies\r(average, %)',
                              name='Kim_2019_Test')
    return data_class


# AsCas12a datasets
def load_kim_2018_train():
    data = pd.read_csv(os.path.join(curr_path(), 'data/datasets/Kim_2018_Train.csv.zip'))
    data_class = ActivityData(data=data, enzyme=en.cas12a,
                              kmer_column='Context Sequence',
                              activity_column='Indel frequency',
                              name='Kim_2018_Train')
    return data_class


def load_kim_2018_test():
    data = pd.read_csv(os.path.join(curr_path(), 'data/datasets/Kim_2018_Test.csv.zip'))
    data_class = ActivityData(data=data, enzyme=en.cas12a,
                              kmer_column='Context Sequence',
                              activity_column='Indel frequency',
                              name='Kim_2018_Test')
    return data_class

# enAsCas12a datasets
