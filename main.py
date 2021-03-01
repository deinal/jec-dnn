import os
os.environ['TF_CPP_MIN_LOG_LEVEL'] = '2'

import tensorflow as tf
import yaml
from model import get_model
from plotter import plot_loss
from data import create_dataset


if __name__ == '__main__':

    for device in tf.config.list_physical_devices('GPU'):
        print(device)

    with open('config.yaml') as f:
        config = yaml.safe_load(f)

    dataset = create_dataset(features=config['features'], batch_size=config['batch_size'])

    dnn = get_model(num_features=len(config['features']['constituents']))

    dnn.compile(optimizer='adam', loss='mean_absolute_error')

    fit = dnn.fit(dataset, epochs=config['epochs'])

    plot_loss(fit.history)
