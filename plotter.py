import matplotlib.pyplot as plt


def plot_loss(history):
    plt.ylim([0, history['loss'][1]])
    plt.plot(history['loss'], label='training loss')
    plt.plot(history['val_loss'], label='validation_loss')
    plt.xlabel('epoch')
    plt.ylabel('loss')
    plt.legend()
    plt.savefig('plots/loss.png')
