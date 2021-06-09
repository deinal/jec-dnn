import os
os.environ['TF_CPP_MIN_LOG_LEVEL'] = '2'
import shutil
import tensorflow as tf
import argparse
import pickle
import yaml
from deepset import get_deepset
from particle_net import get_particle_net
from data import create_datasets


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
    arg_parser.add_argument('--gpus', nargs='+', required=True, help='GPUs to run on in the form 0 1 etc.')
    args = arg_parser.parse_args()

    os.environ['CUDA_VISIBLE_DEVICES'] = ','.join(args.gpus)
    print('GPU devices:', tf.config.list_physical_devices('GPU'))

    try:
        os.mkdir(args.outdir)
    except FileExistsError:
        pass

    with open('config.yaml') as f:
        config = yaml.safe_load(f)
        net = config['net']

    shutil.copyfile('config.yaml', f'{args.outdir}/config.yaml')

    train_ds, val_ds, test_ds, test_files, metadata = create_datasets(net, args.indir, config['data'])
    num_constituents, num_globals, num_points = metadata

    train_ds = train_ds.shuffle(config['shuffle_buffer'])

    strategy = tf.distribute.MirroredStrategy()
    with strategy.scope():
        if net == 'deepset':
            dnn = get_deepset(num_constituents, num_globals, config['model']['deepset'])
        if net == 'particlenet':
            dnn = get_particle_net(num_constituents, num_points, num_globals, config['model']['particlenet'])

        dnn.compile(optimizer=config['optimizer'], loss=config['loss'])
        dnn.optimizer.lr.assign(config['lr'])

    tf.keras.utils.plot_model(dnn, f'{args.outdir}/model.png', dpi=100, show_shapes=True, expand_nested=True)

    callbacks = get_callbacks(config['callbacks'])

    fit = dnn.fit(train_ds, validation_data=val_ds, epochs=config['epochs'], callbacks=callbacks)

    predictions = dnn.predict(test_ds, use_multiprocessing=True, workers=os.cpu_count()-1)

    # Save predictions and corresponding test files
    with open(f'{args.outdir}/predictions.pkl', 'wb') as f:
        pickle.dump((predictions, test_files), f)

    # Save training history
    with open(f'{args.outdir}/history.pkl', 'wb') as f:
        pickle.dump(fit.history, f)

    # Save model
    dnn.save(f'{args.outdir}/dnn')
