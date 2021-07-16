import os
os.environ['TF_CPP_MIN_LOG_LEVEL'] = '2'
import glob
import shutil
import tensorflow as tf
import argparse
import pickle
import yaml
from src.deepset import get_deepset
from src.particle_net import get_particle_net
from src.data import create_datasets


def get_callbacks(config):
    # Reduce learning rate when nearing convergence
    reduce_lr_on_plateau = tf.keras.callbacks.ReduceLROnPlateau(
        monitor='val_loss', factor=config['reduce_lr_on_plateau']['factor'], 
        patience=config['reduce_lr_on_plateau']['patience'], min_lr=config['reduce_lr_on_plateau']['min_lr'],
        mode='auto', min_delta=config['reduce_lr_on_plateau']['min_delta'], cooldown=0, verbose=1
    )
    # Stop early if the network stops improving
    early_stopping = tf.keras.callbacks.EarlyStopping(
        monitor='val_loss', min_delta=config['early_stopping']['min_delta'], 
        patience=config['early_stopping']['patience'], mode='auto', baseline=None, 
        restore_best_weights=True, verbose=1
    )

    return [reduce_lr_on_plateau, early_stopping]


if __name__ == '__main__':
    arg_parser = argparse.ArgumentParser(description=__doc__)
    arg_parser.add_argument('-i', '--indir', required=True, help='Directory containing jet data')
    arg_parser.add_argument('-o', '--outdir', required=True, help='Where to store outputs')
    arg_parser.add_argument('-c', '--config', required=True, help='Config file')
    arg_parser.add_argument('--gpus', nargs='+', required=True, help='GPUs to run on in the form 0 1 etc.')
    arg_parser.add_argument('--save-model', action='store_true', help='If model should be saved')
    arg_parser.add_argument('--unzip', action='store_true', help='Whether to unzip the input dataset or not')
    args = arg_parser.parse_args()

    os.environ['CUDA_VISIBLE_DEVICES'] = ','.join(args.gpus)
    print('GPU devices:', tf.config.list_physical_devices('GPU'))

    try:
        os.makedirs(args.outdir)
    except FileExistsError:
        pass

    with open(args.config) as f:
        config = yaml.safe_load(f)
        net = config['net']

    shutil.copyfile(args.config, f'{args.outdir}/config.yaml')

    if len(glob.glob(os.path.join(args.indir, '*.root'))):
        train_ds, val_ds, test_ds, metadata = create_datasets(net, args.indir, config['data'])
    else:
        with open(f'{args.indir}/metadata.pkl', 'rb') as f:
            metadata = pickle.load(f)

        if args.unzip:
            compression = 'GZIP'
        else:
            compression = None

        train_ds = tf.data.experimental.load(
            os.path.join(args.indir, 'train'),
            element_spec=metadata['element_spec'], compression=compression
        )
        val_ds = tf.data.experimental.load(
            os.path.join(args.indir, 'val'),
            element_spec=metadata['element_spec'], compression=compression
        )
        test_ds = tf.data.experimental.load(
            os.path.join(args.indir, 'test'),
            element_spec=metadata['element_spec'], compression=compression
        )

    train_ds = train_ds.shuffle(config['shuffle_buffer'])

    strategy = tf.distribute.MirroredStrategy()
    with strategy.scope():
        if net == 'deepset':
            dnn = get_deepset(
                metadata['num_constituents'], metadata['num_globals'], 
                config['model']['deepset']
            )
        if net == 'particlenet':
            dnn = get_particle_net(
                metadata['num_constituents'], metadata['num_points'], metadata['num_globals'], 
                config['model']['particlenet']
            )

        dnn.compile(optimizer=config['optimizer'], loss=config['loss'])
        dnn.optimizer.lr.assign(config['lr'])

    tf.keras.utils.plot_model(dnn, os.path.join(args.outdir, 'model.png'), dpi=100, show_shapes=True, expand_nested=True)

    callbacks = get_callbacks(config['callbacks'])

    fit = dnn.fit(train_ds, validation_data=val_ds, epochs=config['epochs'], callbacks=callbacks)

    predictions = dnn.predict(test_ds, use_multiprocessing=True, workers=os.cpu_count()-1)

    # Save predictions and corresponding test files
    with open(os.path.join(args.outdir, 'predictions.pkl'), 'wb') as f:
        pickle.dump((predictions, metadata['test_files']), f)

    # Save training history
    with open(os.path.join(args.outdir, 'history.pkl'), 'wb') as f:
        pickle.dump(fit.history, f)

    # Save model
    if args.save_model:
        dnn.save(os.path.join(args.outdir, 'dnn'))
