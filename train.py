import os
os.environ['TF_CPP_MIN_LOG_LEVEL'] = '2'

import tensorflow as tf
import pickle
import yaml
import multiprocessing
from model import get_model
from data import create_datasets


def get_callbacks(config):
    # Reduce learning rate when nearing convergence
    reduce_lr_on_plateau = tf.keras.callbacks.ReduceLROnPlateau(
        monitor='val_loss', factor=config['factor'], patience=config['patience'], min_lr=config['min_lr'],
        mode='auto', min_delta=config['min_delta'], cooldown=0, verbose=1
    )
    # Stop early if the network stops improving
    early_stopping = tf.keras.callbacks.EarlyStopping(
        monitor='val_loss', min_delta=config['min_delta'], patience=config['patience'], 
        mode='auto', baseline=None, restore_best_weights=True, verbose=1
    )

    return [reduce_lr_on_plateau, early_stopping]


if __name__ == '__main__':
    for device in tf.config.list_physical_devices('GPU'):
        print(device)

    with open('config.yaml') as f:
        config = yaml.safe_load(f)

    outdir = config['outdir']
    try:
        os.mkdir(f'results/{outdir}')
    except FileExistsError:
        pass

    train, validation, test, test_files = create_datasets(config['indir'], config['data'])

    train = train.shuffle(100)

    num_constituents = len(config['data']['features']['jet_pf_cands']) + len(config['data']['features']['pf_cands'])
    num_globals = len(config['data']['features']['jets'])

    strategy = tf.distribute.MirroredStrategy()
    with strategy.scope():
        dnn = get_model(outdir, num_constituents, num_globals, config['model'])
        dnn.compile(optimizer=config['optimizer'], loss=config['loss'])

    callbacks = get_callbacks(config['callbacks'])

    fit = dnn.fit(train, validation_data=validation, epochs=config['epochs'], callbacks=callbacks)

    predictions = dnn.predict(test, use_multiprocessing=True, workers=multiprocessing.cpu_count())

    # Save predictions and corresponding test files
    with open(f'./results/{outdir}/predictions.pkl', 'wb') as f:
        pickle.dump((predictions, test_files), f)

    # Save training history
    with open(f'./results/{outdir}/history.pkl', 'wb') as f:
        pickle.dump(fit.history, f)

    # Save model
    dnn.save(f'./results/{outdir}/dnn')
