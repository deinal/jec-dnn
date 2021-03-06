import glob
import tensorflow as tf
import numpy as np
import pandas as pd
import uproot
import awkward as ak

def create_datasets(indir, features, batch_size, train_size, test_size):

    pickle_paths = glob.glob(f'./data/{indir}/*.root')
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

    inp = [
        path, features['jets'], features['gen_jets'], features['jet_pf_cands'], features['pf_cands']
    ]
    Tout = (
        [tf.int32] + [tf.float32] +
        [tf.float32] * len(features['jets']) + [tf.float32] * len(features['gen_jets']) +
        [tf.float32] * len(features['jet_pf_cands']) + [tf.float32] * len(features['pf_cands'])
    )

    data = tf.numpy_function(_read_nanoaod, inp=inp, Tout=Tout)

    cols = (
        ['row_lengths'] + ['target'] + 
        features['jets'] + features['gen_jets'] + 
        features['jet_pf_cands'] + features['pf_cands']
    )
    data = {key: value for key, value in zip(cols, data)}

    # change each column's shape from <unknown> to (None,)
    for values in data.values():
       values.set_shape((None,))

    row_lengths = data.pop('row_lengths')
    target = data.pop('target')

    constituent_cols = features['jet_pf_cands'] + features['pf_cands']

    for constituent in constituent_cols:
        data[constituent] = tf.RaggedTensor.from_row_lengths(data[constituent], row_lengths=row_lengths)
        data[constituent] = tf.expand_dims(data[constituent], axis=2)

    inputs = tf.concat([data[name] for name in constituent_cols], axis=2)

    return (inputs, target)


def _read_nanoaod(path, jets_cols, gen_jets_cols, jet_pf_cands_cols, pf_cands_cols):

    # Decode bytestrings
    path = path.decode()
    jets_cols = [col.decode() for col in jets_cols]
    gen_jets_cols = [col.decode() for col in gen_jets_cols]
    jet_pf_cands_cols = [col.decode() for col in jet_pf_cands_cols]
    pf_cands_cols = [col.decode() for col in pf_cands_cols]

    all_features = jets_cols + gen_jets_cols + jet_pf_cands_cols + pf_cands_cols

    # open root file
    with uproot.open(path) as file:
        events = file['Events'].arrays(all_features)

    # Select global observables
    globals = events[jets_cols + gen_jets_cols]
    
    # Reorder gen jets according to the jet collection
    for field in gen_jets_cols:
        globals[field] = globals[field][events.Jet_genJetIdx]

    # Select leading jets and concatenate vertically to have one jet's globals per row
    globals = ak.concatenate((globals[:,0:1], globals[:,1:2]), axis=0)
    globals = ak.to_pandas(globals)

    # Select pf candidates
    constituents = events[jet_pf_cands_cols + pf_cands_cols]

    # Reorder pf cands according to the jet pf cands collection
    for field in pf_cands_cols:
        constituents[field] = constituents[field][events.JetPFCands_pFCandsIdx]
    
    # Select leading jets and concatenate vertically to have one jet's constituents per row
    constituents = ak.concatenate(
        (constituents[constituents.JetPFCands_jetIdx == 0], constituents[constituents.JetPFCands_jetIdx == 1]), axis=0
    )
    constituents = ak.to_pandas(constituents)

    # Calculate regression target
    target = globals['GenJet_pt'] / globals['Jet_pt']

    # Calculate the lengths of every set of constituents
    _, row_lengths = np.unique(constituents.index.get_level_values(0), return_counts=True)

    # Construct a list of arrays with the desired data
    data = [row_lengths.astype(np.int32)] + [target.astype(np.float32)]

    for col in jets_cols + gen_jets_cols:
        data.append(globals[col].astype(np.float32))

    for col in jet_pf_cands_cols + pf_cands_cols:
        data.append(constituents[col].astype(np.float32))

    return data
