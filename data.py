import glob
import warnings
import tensorflow as tf
import awkward as ak
from coffea.nanoevents import NanoEventsFactory, PFNanoAODSchema
import tensorflow_transform as tft


def create_datasets(indir, config):
    root_paths = glob.glob(f'{indir}/*.root')
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
    global_fields = [f'jet_{field}' for field in features['jets']['numerical']]
    constituent_fields = [f'pf_cand_{field}' for field in features['pf_cands']['numerical']]

    dataset = tf.data.Dataset.from_tensor_slices(files)

    dataset = dataset.map(
        lambda path: _retrieve_data(path, features, global_fields, constituent_fields),
        num_parallel_calls=tf.data.experimental.AUTOTUNE,
    )

    dataset = dataset.map(
        lambda data, target: (_preprocess(data, global_fields, constituent_fields), target),
        num_parallel_calls=tf.data.experimental.AUTOTUNE,
    )

    dataset = dataset.unbatch().batch(batch_size)

    dataset = dataset.prefetch(tf.data.experimental.AUTOTUNE)

    return dataset


def _preprocess(data, global_fields, constituent_fields):
    for field in constituent_fields:
        flat_constituent = data[field].flat_values
        normalized_constituent = tft.scale_to_gaussian(flat_constituent)
        data[field] = tf.RaggedTensor.from_row_splits(normalized_constituent, data[field].row_splits)

    for field in global_fields:
        data[field] = tft.scale_to_gaussian(data[field])
    
    constituents = tf.concat([data[field] for field in constituent_fields], axis=2)

    globals = tf.concat([data[field] for field in global_fields], axis=1)

    # mind the order of the inputs when constructing the neural net 
    inputs = (constituents, globals)

    return inputs


def _retrieve_data(path, features, global_fields, constituent_fields):
    inp = [
        path, features['jets']['numerical'], features['pf_cands']['numerical']
    ]
    Tout = (
        [tf.int64] + [tf.float32] +
        [tf.float32] * len(features['jets']['numerical']) +
        [tf.float32] * len(features['pf_cands']['numerical'])
    )

    data = tf.numpy_function(_read_nanoaod, inp=inp, Tout=Tout)

    fields = ['row_lengths'] + ['target'] + global_fields + constituent_fields

    data = {key: value for key, value in zip(fields, data)}

    for values in data.values():
        # shape from <unknown> to (None,)
        values.set_shape((None,))

    row_lengths = data.pop('row_lengths')
    target = data.pop('target')

    for field in global_fields:
        # shape from (None,) to (None, 1)
        data[field] = tf.expand_dims(data[field], axis=1)

    for field in constituent_fields:
        # shape from (None,) to (None, None)
        data[field] = tf.RaggedTensor.from_row_lengths(data[field], row_lengths=row_lengths)
        # shape from (None, None) to (None, None, 1)
        data[field] = tf.expand_dims(data[field], axis=2)

    return (data, target)


def _read_nanoaod(path, jet_fields, pf_cand_fields):
    # Decode bytestrings
    path = path.decode()
    jet_fields = [field.decode() for field in jet_fields]
    pf_cand_fields = [field.decode() for field in pf_cand_fields]

    with warnings.catch_warnings():
        warnings.filterwarnings('ignore', message='Found duplicate branch')
        warnings.filterwarnings('ignore', message='Missing cross-reference index')
        events = NanoEventsFactory.from_root(path, schemaclass=PFNanoAODSchema).events()

    jets = events.Jet[(ak.count(events.Jet.matched_gen.pt, axis=1) >= 2)]

    sorted_jets = jets[ak.argsort(jets.matched_gen.pt, ascending=False, axis=1)]

    leading_jets = ak.concatenate((sorted_jets[:,0], sorted_jets[:,1]), axis=0)

    selected_jets = leading_jets[(leading_jets.matched_gen.pt > 30) & (abs(leading_jets.matched_gen.eta) < 2.5)]

    valid_jets = selected_jets[~ak.is_none(selected_jets.matched_gen.pt)]

    target = valid_jets.matched_gen.pt / valid_jets.pt

    globals = valid_jets[jet_fields]

    constituents = valid_jets.constituents.pf[pf_cand_fields]
    row_lengths = ak.num(constituents, axis=1)
    flat_constituents = ak.flatten(constituents, axis=1)

    # Construct a list of numpy arrays with the desired data
    data = [ak.to_numpy(row_lengths)] + [ak.to_numpy(target)]

    for field in jet_fields:
        data.append(ak.to_numpy(globals[field]))

    for field in pf_cand_fields:
        data.append(ak.to_numpy(flat_constituents[field]))

    return data
