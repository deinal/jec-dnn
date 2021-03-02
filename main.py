import os
os.environ['TF_CPP_MIN_LOG_LEVEL'] = '2'

import tensorflow as tf
import yaml
from model import get_model
from plotter import plot_loss
from data import create_datasets


if __name__ == '__main__':

    for device in tf.config.list_physical_devices('GPU'):
        print(device)

    with open('config.yaml') as f:
        config = yaml.safe_load(f)

    train, test, validation = create_datasets(
        features=config['features'], batch_size=config['batch_size'],
        train_size=config['train_size'], test_size=config['test_size']
    )

    train = train.shuffle(100)

    dnn = get_model(num_features=len(config['features']['constituents']))

    dnn.compile(optimizer='adam', loss='mean_absolute_error')

    fit = dnn.fit(train, validation_data=validation, epochs=config['epochs'])

    plot_loss(fit.history)

    predictions = dnn.predict(test)
