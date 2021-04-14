import os
import tensorflow as tf
import argparse
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
    arg_parser = argparse.ArgumentParser(description=__doc__)
    arg_parser.add_argument('-i', '--indir', required=True, help='Directory containing jet data')
    arg_parser.add_argument('-o', '--outdir', required=True, help='Where to store outputs')
    args = arg_parser.parse_args()

    try:
        os.mkdir(f'{args.outdir}')
    except FileExistsError:
        pass

    with open('config.yaml') as f:
        config = yaml.safe_load(f)

    train_ds, val_ds, test_ds, test_files = create_datasets(args.indir, config['data'])

    # train_ds = train_ds.shuffle(100)

    num_constituents = len(config['data']['features']['pf_cands'])
    num_globals = len(config['data']['features']['jets'])

    dnn = get_model(num_constituents, num_globals, config['model'])
    dnn.compile(optimizer=config['optimizer'], loss=config['loss'])

    callbacks = get_callbacks(config['callbacks'])

    fit = dnn.fit(train_ds, validation_data=val_ds, epochs=config['epochs'], callbacks=callbacks)

    predictions = dnn.predict(test_ds, use_multiprocessing=True, workers=multiprocessing.cpu_count())

    # Save predictions and corresponding test files
    with open(f'{args.outdir}/predictions.pkl', 'wb') as f:
        pickle.dump((predictions, test_files), f)

    # Save training history
    with open(f'{args.outdir}/history.pkl', 'wb') as f:
        pickle.dump(fit.history, f)

    # Save model
    dnn.save(f'{args.outdir}/dnn')
