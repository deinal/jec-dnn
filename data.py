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

    # Concatenate leading jets vertically
    jets = events[jets_cols]
    jets = ak.concatenate((jets[:,0:1], jets[:,1:2]), axis=0)
    jets_df = ak.to_pandas(jets)

    # Sort gen jets according to the jet collection and concatenate leading gen jets vertically
    gen_jets = events[gen_jets_cols][events.Jet_genJetIdx]
    gen_jets = ak.concatenate((gen_jets[:,0:1], gen_jets[:,1:2]), axis=0)
    gen_jets_df = ak.to_pandas(gen_jets)

    # Concatentate jets and gen jets horizontally
    globals = pd.concat((jets_df, gen_jets_df), axis=1)

    # Select jet pf cands
    jet_pf_cands = events[jet_pf_cands_cols]
    jet_pf_cands_df = ak.to_pandas(jet_pf_cands)

    # Sort pf cands according to jet pf cands collection
    pf_cands = events[pf_cands_cols][events.JetPFCands_pFCandsIdx]
    pf_cands_df = ak.to_pandas(pf_cands)

    # Concatenate pf cands and jet pf cands horizontally
    constituents = pd.concat((pf_cands_df, jet_pf_cands_df), axis=1)

    # Make a tuple of leading jet constituents and next-to-leading jet constituents 
    constituents = (constituents[constituents.JetPFCands_jetIdx == 0], constituents[constituents.JetPFCands_jetIdx == 1])

    # Create a new outer multi index for next-to-leading jet constituents in preparation for concatenation
    shift = constituents[0].index[-1][0] + 1
    outer_index = constituents[1].index.get_level_values(0) + shift
    inner_index = constituents[1].index.get_level_values(1)
    constituents[1].index = pd.MultiIndex.from_tuples(zip(outer_index, inner_index), names=('entry', 'subentry'))

    # Concatenate leading jet constituents and next-to-leading jet constituents vertically
    constituents = pd.concat(constituents, axis=0)

    # Calculate regression target
    target = globals['GenJet_pt'] / globals['Jet_pt']

    # Calculate the lengths of every constituent subset
    _, row_lengths = np.unique(constituents.index.get_level_values(0), return_counts=True)

    # Construct a list of arrays with the desired data
    data = [row_lengths.astype(np.int32)] + [target.astype(np.float32)]

    for col in jets_cols + gen_jets_cols:
        data.append(globals[col].astype(np.float32))

    for col in jet_pf_cands_cols + pf_cands_cols:
        data.append(constituents[col].astype(np.float32))

    return data
