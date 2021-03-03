import glob
import pickle
import tensorflow as tf
import numpy as np


def create_datasets(indir, features, batch_size, train_size, test_size):

    pickle_paths = glob.glob(f'./data/{indir}/*.pkl')
    num_files = len(pickle_paths)
    train_split = int(train_size * num_files)
    test_split = int(test_size * num_files) + train_split
    
    train_files = pickle_paths[:train_split]
    test_files = pickle_paths[train_split:test_split]
    validation_files = pickle_paths[test_split:]

    train = _create_dataset(train_files, features, batch_size)
    test = _create_dataset(test_files, features, batch_size)
    validation = _create_dataset(validation_files, features, batch_size)

    return train, validation, test, test_files


def _create_dataset(files, features, batch_size):
    
    dataset = tf.data.Dataset.from_tensor_slices(files)

    dataset = dataset.map(
        lambda path: _retrieve_data(path, features),
        num_parallel_calls=tf.data.experimental.AUTOTUNE,
    )

    dataset = dataset.unbatch().batch(batch_size)

    dataset = dataset.prefetch(tf.data.experimental.AUTOTUNE)

    return dataset


def _retrieve_data(path, features):
    inp = [path, features['globals'], features['constituents']]
    Tout = [tf.int32] + [tf.float32] + [tf.float32] * len(features['globals']) + [tf.float32] * len(features['constituents'])

    data = tf.numpy_function(_read_pickle, inp=inp, Tout=Tout)

    columns = ['row_lengths'] + ['target'] + features['globals'] + features['constituents']
    data = {key: value for key, value in zip(columns, data)}

    # change each column's shape from <unknown> to (None,)
    for values in data.values():
       values.set_shape((None,))

    row_lengths = data.pop('row_lengths')
    target = data.pop('target')

    for constituent in features['constituents']:
        data[constituent] = tf.RaggedTensor.from_row_lengths(data[constituent], row_lengths=row_lengths)
        data[constituent] = tf.expand_dims(data[constituent], axis=2)

    inputs = tf.concat([data[name] for name in features['constituents']], axis=2)

    return (inputs, target)


def _read_pickle(path, globals_cols, constituents_cols):

    path = path.decode()
    globals_cols = globals_cols.astype(str)
    constituents_cols = constituents_cols.astype(str)

    with open(path, 'rb') as f:
        data = pickle.load(f)

    globals = data['globals'][globals_cols]
    constituents = data['constituents'][constituents_cols]

    target = globals['GenJet_pt'] / globals['Jet_pt']

    _, row_lengths = np.unique(constituents.index.get_level_values(0), return_counts=True)

    data = [row_lengths.astype(np.int32)] + [target.astype(np.float32)]

    for col in globals_cols:
        data.append(globals[col].astype(np.float32))

    for col in constituents_cols:
        data.append(constituents[col].astype(np.float32))

    return data
