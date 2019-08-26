import matplotlib.pyplot as plt
import pandas as pd

def plot_model_history(model_history):
    plt.plot(pd.DataFrame(model_history.history).loc[:, ['loss', 'val_loss']])
    plt.xlabel('Epochs')
    plt.ylabel('Loss')
