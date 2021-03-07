import os
os.environ['TF_CPP_MIN_LOG_LEVEL'] = '2'

import tensorflow as tf
import pickle
import yaml
from model import get_model
from data import create_datasets


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

    train, validation, test, test_files = create_datasets(
        indir=config['indir'], features=config['features'], batch_size=config['batch_size'],
        train_size=config['train_size'], test_size=config['test_size']
    )

    train = train.shuffle(100)

    num_features = len(config['features']['jet_pf_cands']) + len(config['features']['pf_cands'])

    dnn = get_model(outdir=outdir, num_features=num_features)

    dnn.compile(optimizer='adam', loss='mean_absolute_error')

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

    fit = dnn.fit(train, validation_data=validation, epochs=config['epochs'], callbacks=[reduce_lr_on_plateau, early_stopping])

    predictions = dnn.predict(test)

    # Save predictions and corresponding test files
    with open(f'./results/{outdir}/predictions.pkl', 'wb') as f:
        pickle.dump((predictions, test_files), f)

    # Save training history
    with open(f'./results/{outdir}/history.pkl', 'wb') as f:
        pickle.dump(fit.history, f)
