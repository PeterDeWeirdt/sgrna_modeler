import pandas as pd
from sgrna_modeler import enzymes as en
import os

def curr_path():
    return os.path.dirname(__file__)

class Activity_Data(object):
    def __init__(self, data, enzyme, kmer_column, activity_column, name,
                 group_column = ''):
        self.data = data
        self.enzyme = enzyme
        self.kmer_column = kmer_column
        self.activity_column = activity_column
        self.name = 'D_' + name
        self.group_column = group_column
    def get_xy(self):
        x = self.data[self.kmer_column]
        y = self.data[self.activity_column]
        return x,y
    def get_groups(self):
        if not self.group_column:
            assert ValueError('No Group Column Supplied')
        return self.data[self.group_column]

# SpCas9 Datasets

def load_doench_2016():
    data = pd.read_csv(os.path.join(curr_path(), 'data/datasets/Doench_2016.csv.zip'))
    data_class = Activity_Data(data = data, enzyme = en.cas9, kmer_column='30mer',
                               activity_column='score_drug_gene_rank',
                               name = 'Doench_2016',
                               group_column='Target gene')
    return data_class

def load_meyers_2017_train():
    data = pd.read_csv(os.path.join(curr_path(), 'data/datasets/Meyers_2017_Train.csv.zip'))
    data_class = Activity_Data(data = data, enzyme = en.cas9, kmer_column='sgRNA context sequence',
                               activity_column='mean_activity',
                               name = 'Meyers_2017_Train',
                               group_column = 'Gene Symbol')
    return data_class

def load_meyers_2017_test():
    data = pd.read_csv(os.path.join(curr_path(),'data/datasets/Meyers_2017_Test.csv.zip'))
    data_class = Activity_Data(data = data, enzyme=en.cas9, kmer_column='sgRNA context sequence',
                               activity_column='mean_activity',
                               name='Meyers_2017_Test',
                               group_column='Gene Symbol')
    return data_class

def load_kim_2019_train():
    data = pd.read_csv(os.path.join(curr_path(),'data/datasets/Kim_2019_Train.csv.zip'))
    data_class = Activity_Data(data = data, enzyme = en.cas9,
                               kmer_column='Target context sequence (4+20+3+3)',
                               activity_column='Background subtracted indel (%)',
                               name = 'Kim_2019_Train')
    return data_class

def load_kim_2019_test():
    data = pd.read_csv(os.path.join(curr_path(),'data/datasets/Kim_2019_Test.csv.zip'))
    data_class = Activity_Data(data = data, enzyme = en.cas9,
                               kmer_column='Target context sequence (4+20+3+3)',
                               activity_column='Background subtracted indel frequencies\r(average, %)',
                               name = 'Kim_2019_Test')
    return data_class

# AsCas12a datasets
def load_kim_2018_train():
    data = pd.read_csv(os.path.join(curr_path(),'data/datasets/Kim_2018_Train.csv.zip'))
    data_class = Activity_Data(data = data, enzyme = en.cas12a,
                               kmer_column='Context Sequence',
                               activity_column='Indel frequency',
                               name = 'Kim_2018_Train')
    return data_class

def load_kim_2018_test():
    data = pd.read_csv(os.path.join(curr_path(),'data/datasets/Kim_2018_Test.csv.zip'))
    data_class = Activity_Data(data = data, enzyme = en.cas12a,
                               kmer_column='Context Sequence',
                               activity_column='Indel frequency',
                               name = 'Kim_2019_Test')
    return data_class

# enAsCas12a datasets


