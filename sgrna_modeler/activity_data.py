# library

class Activity_Data(object):
    def __init__(self, data, enzyme, kmer_column, activity_column,
                 group_column = ''):
        self.data = data
        self.enzyme = enzyme
        self.kmer_column = kmer_column
        self.activity_column = activity_column
        self.group_column = group_column
    def get_xy(self):
        x = self.data[self.kmer_column]
        y = self.data[self.activity_column]
        return x,y
    def get_groups(self):
        if not self.group_column:
            assert ValueError('No Group Column Supplied')
        return self.data[self.group_column]

