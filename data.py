import glob
import tensorflow as tf
import awkward as ak
from coffea.nanoevents import NanoEventsFactory, NanoAODSchema


def create_datasets(indir, config):
    root_paths = glob.glob(f'./data/{indir}/*.root')
    num_files = len(root_paths)
    train_split = int(config['train_size'] * num_files)
    test_split = int(config['test_size'] * num_files) + train_split
    
    train_files = root_paths[:train_split]
    test_files = root_paths[train_split:test_split]
    validation_files = root_paths[test_split:]

    train = _create_dataset(train_files, config['features'], config['batch_size'])
    test = _create_dataset(test_files, config['features'], config['batch_size'])
    validation = _create_dataset(validation_files, config['features'], config['batch_size'])

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
    inp = [
        path, features['jets'], features['pf_cands']
    ]
    Tout = (
        [tf.int64] + [tf.float32] +
        [tf.float32] * len(features['jets']) +
        [tf.float32] * len(features['pf_cands'])
    )

    data = tf.numpy_function(_read_nanoaod, inp=inp, Tout=Tout)

    globals_cols = [f'jet_{field}' for field in features['jets']]
    constituent_cols = [f'pf_cand_{field}' for field in features['pf_cands']]

    cols = ['row_lengths'] + ['target'] + globals_cols + constituent_cols

    data = {key: value for key, value in zip(cols, data)}

    for values in data.values():
        # shape from <unknown> to (None,)
        values.set_shape((None,))

    row_lengths = data.pop('row_lengths')
    target = data.pop('target')

    inputs = {}

    for col in globals_cols:
        # shape from (None,) to (None, 1)
        data[col] = tf.expand_dims(data[col], axis=1)

    inputs['globals'] = tf.concat([data[col] for col in globals_cols], axis=1)

    for col in constituent_cols:
        # shape from (None,) to (None, None)
        data[col] = tf.RaggedTensor.from_row_lengths(data[col], row_lengths=row_lengths)
        # shape from (None, None) to (None, None, 1)
        data[col] = tf.expand_dims(data[col], axis=2)

    inputs['constituents'] = tf.concat([data[col] for col in constituent_cols], axis=2)

    return (inputs, target)


def _read_nanoaod(path, jets_cols, pf_cands_cols):
    # Decode bytestrings
    path = path.decode()
    jets_cols = [col.decode() for col in jets_cols]
    pf_cands_cols = [col.decode() for col in pf_cands_cols]

    events = NanoEventsFactory.from_root(path, schemaclass=NanoAODSchema).events()

    sorted_jets = events.Jet[ak.argsort(events.Jet.pt, ascending=False, axis=1)]

    leading_jets = ak.concatenate((sorted_jets[:,0:1], sorted_jets[:,1:2]), axis=0)
    leading_jets = ak.flatten(leading_jets, axis=1)

    globals = leading_jets[jets_cols]
    target = leading_jets.matched_gen.pt / leading_jets.pt

    constituents = leading_jets.constituents.pf[pf_cands_cols]
    row_lengths = ak.num(constituents, axis=1)
    constituents = ak.flatten(constituents, axis=1)

    # Construct a list of numpy arrays with the desired data
    data = [ak.to_numpy(row_lengths)] + [ak.to_numpy(target)]

    for col in jets_cols:
        data.append(ak.to_numpy(globals[col]))

    for col in pf_cands_cols:
        data.append(ak.to_numpy(constituents[col]))

    return data
