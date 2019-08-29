import pandas as pd
from sgrna_modeler import enzymes as en
import os

def curr_path():
    return os.path.dirname(__file__)

class Activity_Data(object):
    def __init__(self, data, enzyme, kmer_column, activity_column, name, type,
                 group_column = ''):
        self.data = data
        self.enzyme = enzyme
        self.kmer_column = kmer_column
        self.activity_column = activity_column
        self.type = type
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

# Sp Datasets
def load_doench_2016_train():
    data = pd.read_csv(os.path.join(curr_path(), 'data/Doench_2016_Train.csv'))
    data_class = Activity_Data(data = data, enzyme = en.cas9, kmer_column='30mer',
                               activity_column='score_drug_gene_rank',
                               name = 'Doench_2016_Train',
                               type = 'train',
                               group_column='Target gene')
    return data_class

def load_doench_2016_test():
    data = pd.read_csv(os.path.join(curr_path(), 'data/Doench_2016_Test.csv'))
    data_class = Activity_Data(data = data, enzyme = en.cas9, kmer_column='30mer',
                               activity_column='score_drug_gene_rank',
                               name = 'Doench_2016_Test',
                               type = 'test',
                               group_column='Target gene')
    return data_class

def load_meyers_2017_train():
    data = pd.read_csv(os.path.join(curr_path(), 'data/Meyers_2017_Train.csv'))
    data_class = Activity_Data(data = data, enzyme = en.cas9, kmer_column='sgRNA context sequence',
                               activity_column='mean_activity',
                               name = 'Meyers_2017_Train',
                               type = 'train',
                               group_column='Gene Symbol')
    return data_class

def load_meyers_2017_test():
    data = pd.read_csv(os.path.join(curr_path(),'data/Meyers_2017_Test.csv'))
    data_class = Activity_Data(data = data, enzyme = en.cas9, kmer_column='sgRNA context sequence',
                               activity_column='mean_activity',
                               name = 'Meyers_2017_Test',
                               type = 'test',
                               group_column='Gene Symbol')
    return data_class

def load_kim_2019_train():
    data = pd.read_csv(os.path.join(curr_path(),'data/Kim_2019_Train.csv'))
    data_class = Activity_Data(data = data, enzyme = en.cas9,
                               kmer_column='Target context sequence (4+20+3+3)',
                               activity_column='Background subtracted indel (%)',
                               name = 'Kim_2019_Train',
                               type = 'train')
    return data_class

def load_kim_2019_test():
    data = pd.read_csv(os.path.join(curr_path(),'data/Kim_2019_Test.csv'))
    data_class = Activity_Data(data = data, enzyme = en.cas9,
                               kmer_column='Target context sequence (4+20+3+3)',
                               activity_column='Background subtracted indel frequencies\r(average, %)',
                               name = 'Kim_2019_Test',
                               type = 'test')
    return data_class

