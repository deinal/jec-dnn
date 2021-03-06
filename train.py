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

    fit = dnn.fit(train, validation_data=validation, epochs=config['epochs'])

    predictions = dnn.predict(test)

    # Save predictions and corresponding test files
    with open(f'./results/{outdir}/predictions.pkl', 'wb') as f:
        pickle.dump((predictions, test_files), f)

    # Save training history
    with open(f'./results/{outdir}/history.pkl', 'wb') as f:
        pickle.dump(fit.history, f)
