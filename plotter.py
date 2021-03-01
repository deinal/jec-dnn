
import matplotlib.pyplot as plt

def plot_loss(history):
    plt.plot(history['loss'])
    plt.xlabel('epoch')
    plt.ylabel('training loss')
    plt.savefig('plots/loss.png')
