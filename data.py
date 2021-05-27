import glob
import warnings
import tensorflow as tf
import awkward as ak
import numpy as np
import math
from coffea.nanoevents import NanoEventsFactory, PFNanoAODSchema


def create_datasets(net, indir, config):
    root_paths = glob.glob(f'{indir}/*.root')
    num_files = len(root_paths)
    train_split = int(config['train_size'] * num_files)
    test_split = int(config['test_size'] * num_files) + train_split
    
    train_files = root_paths[:train_split]
    test_files = root_paths[train_split:test_split]
    val_files = root_paths[test_split:]

    train = _create_dataset(net, train_files, config['features'], config['batch_size'], config['num_points'])
    test = _create_dataset(net, test_files, config['features'], config['batch_size'], config['num_points'])
    val = _create_dataset(net, val_files, config['features'], config['batch_size'], config['num_points'])
    
    metadata = _get_metadata(config)

    return train, val, test, test_files, metadata


def _get_metadata(config):
    features = config['features']
    num_constituents = (
        len(features['pf_cands']['numerical']) + 
        len(sum(features['pf_cands']['categorical'].values(), [])) + 
        len(features['pf_cands']['synthetic'])
    )
    num_globals = (
        len(features['jets']['numerical']) + 
        len(sum(features['jets']['categorical'].values(), []))
    )
    num_points = config['num_points']

    return (num_constituents, num_globals, num_points)


def _create_dataset(net, files, features, batch_size, num_points):
    global_numerical = features['jets']['numerical']
    global_categorical = features['jets']['categorical']
    constituent_numerical = features['pf_cands']['numerical']
    constituent_categorical = features['pf_cands']['categorical']
    constituent_synthetic = features['pf_cands']['synthetic']

    dataset = tf.data.Dataset.from_tensor_slices(files)

    dataset = dataset.map(
        lambda path: _retrieve_data(
            net, path, num_points, global_numerical, list(global_categorical.keys()), 
            constituent_numerical, list(constituent_categorical.keys())
        ),
        num_parallel_calls=4 # a fixed number instead of autotune limits the RAM usage
    )

    global_tables = _create_category_tables(global_categorical)
    constituent_tables = _create_category_tables(constituent_categorical)

    dataset = dataset.map(
        lambda data, target: (_preprocess(
            net, data, global_tables, constituent_tables, global_numerical, global_categorical, 
            constituent_numerical, constituent_categorical, constituent_synthetic), target
        ),
        num_parallel_calls=tf.data.AUTOTUNE
    )

    dataset = dataset.unbatch().batch(batch_size)

    dataset = dataset.prefetch(tf.data.AUTOTUNE)

    return dataset


def _preprocess(
        net, data, global_tables, constituent_tables, global_numerical, global_categorical, 
        constituent_numerical, constituent_categorical, constituent_synthetic
    ):

    if 'rel_pt' in constituent_synthetic:
        data['pf_cand_rel_pt'] = data['pf_cand_pt'] / tf.expand_dims(data['jet_pt'], axis=1)

    if 'rel_mass' in constituent_synthetic:
        data['pf_cand_rel_mass'] = data['pf_cand_mass'] / tf.expand_dims(data['jet_mass'], axis=1)

    if 'rel_eta' in constituent_synthetic:
        jet_eta = tf.expand_dims(data['jet_eta'], axis=1)
        data['pf_cand_rel_eta'] = (data['pf_cand_eta'] - jet_eta) * tf.math.sign(jet_eta)

    if 'rel_phi' in constituent_synthetic:
        phi_diff = data['pf_cand_phi'] - tf.expand_dims(data['jet_phi'], axis=1)
        data['pf_cand_rel_phi'] = (phi_diff + math.pi) % (2 * math.pi) - math.pi

    # create ParticleNet inputs
    if net == 'particlenet':
        mask = tf.cast(tf.math.not_equal(data['pf_cand_pt'], 0), dtype=tf.float32) # 1 if valid
        coord_shift = tf.multiply(1e6, tf.cast(tf.math.equal(mask, 0), dtype=tf.float32))
        points = tf.concat([data['pf_cand_rel_eta'], data['pf_cand_rel_phi']], axis=2)

    for field in global_numerical:
        name = f'jet_{field}'
        data[name] = (data[name] - tf.math.reduce_mean(data[name])) / tf.math.reduce_variance(data[name])

    for field in constituent_numerical + constituent_synthetic:
        name = f'pf_cand_{field}'
        if isinstance(data[name], tf.RaggedTensor):
            flat_constituent = data[name].flat_values
            normalized_constituent = (
                (flat_constituent - tf.math.reduce_mean(flat_constituent)) / tf.math.reduce_variance(flat_constituent)
            )
            data[name] = tf.RaggedTensor.from_row_splits(normalized_constituent, data[name].row_splits)
        else:
            data[name] = (data[name] - tf.math.reduce_mean(data[name])) / tf.math.reduce_variance(data[name])
    
    for field in global_categorical:
        name = f'jet_{field}'
        categories = global_categorical[field]
        encoded_feature = _one_hot_encode(data[name], global_tables[field], categories)
        data[name] = tf.squeeze(encoded_feature, axis=1) # remove excess dimension created by tf.one_hot

    for field in constituent_categorical:
        name = f'pf_cand_{field}'
        categories = constituent_categorical[field]
        encoded_feature = _one_hot_encode(data[name], constituent_tables[field], categories)
        data[name] = tf.squeeze(encoded_feature, axis=2) # remove excess dimension created by tf.one_hot

    globals = tf.concat(
        [data[f'jet_{field}'] for field in global_numerical + list(global_categorical.keys())], axis=1
    )
    constituents = tf.concat(
        [data[f'pf_cand_{field}'] for field in 
        constituent_numerical + list(constituent_categorical.keys()) + constituent_synthetic], axis=2
    )
    
    # mind the order of the inputs when constructing the neural net
    if net == 'deepset':
        inputs = (constituents, globals)
    if net == 'particlenet':
        inputs = (points, constituents, mask, coord_shift)

    return inputs


def _one_hot_encode(feature, table, categories):
    cardinality = len(categories)

    # Map integer categories to an ordered list e.g. charge from [-1, 0, 1] to [0, 1, 2]
    if isinstance(feature, tf.RaggedTensor):
        feature = tf.ragged.map_flat_values(lambda x: table.lookup(x), feature)
    else:
        feature = table.lookup(feature)

    # One-hot encode categories to orthogonal vectors e.g. [0, 1, 2] to [[1, 0, 0], [0, 1, 0], [0, 0, 1]]
    return tf.one_hot(tf.cast(feature, tf.int32), depth=cardinality, dtype=tf.float32)


def _create_category_tables(feature):
    tables = {}

    for field in feature:
        keys_tensor = tf.constant(feature[field])
        vals_tensor = tf.range(len(feature[field]))

        table = tf.lookup.StaticHashTable(
            tf.lookup.KeyValueTensorInitializer(keys_tensor, vals_tensor),
            default_value=-1
        )

        tables[field] = table

    return tables


def _retrieve_data(
        net, path, num_points, global_numerical, global_categorical, 
        constituent_numerical, constituent_categorical
    ):

    global_fields = [
        f'jet_{field}' for field in global_numerical + global_categorical
    ]
    constituent_fields = [
        f'pf_cand_{field}' for field in constituent_numerical + constituent_categorical
    ]
    fields =  ['target'] + global_fields + constituent_fields

    inp = [
        net, path, num_points, global_numerical, global_categorical, 
        constituent_numerical, constituent_categorical
    ]
    Tout = (
        [tf.float32] +
        [tf.float32] * len(global_numerical) +
        [tf.int32] * len(global_categorical) +
        [tf.float32] * len(constituent_numerical) +
        [tf.int32] * len(constituent_categorical)
    )

    if net == 'deepset':
        Tout.append(tf.int32)
        fields.append('row_lengths')

    data = tf.numpy_function(_read_nanoaod, inp=inp, Tout=Tout)

    data = {key: value for key, value in zip(fields, data)}

    target = data.pop('target')
    target.set_shape((None,))

    for field in global_fields:
        # shape from <unknown> to (None,)
        data[field].set_shape((None,))
        # shape from (None,) to (None, 1)
        data[field] = tf.expand_dims(data[field], axis=1)

    if net == 'deepset':
        row_lengths = data.pop('row_lengths')
        row_lengths.set_shape((None,))

        for field in constituent_fields:
            # shape from <unknown> to (None,)
            data[field].set_shape((None,))
            # shape from (None,) to (None, None)
            data[field] = tf.RaggedTensor.from_row_lengths(data[field], row_lengths=row_lengths)
            # shape from (None, None) to (None, None, 1)
            data[field] = tf.expand_dims(data[field], axis=2)
    
    if net == 'particlenet':
        for field in constituent_fields:
            # shape from <unknown> to (None, P)
            data[field].set_shape((None, num_points))
            # shape from (None, P) to (None, P, 1)
            data[field] = tf.expand_dims(data[field], axis=2)

    return (data, target)


def _read_nanoaod(
        net, path, num_points, global_numerical, global_categorical,
        constituent_numerical, constituent_categorical
    ):

    # Decode bytestrings
    net = net.decode()
    path = path.decode()
    global_numerical = [field.decode() for field in global_numerical]
    global_categorical = [field.decode() for field in global_categorical]
    constituent_numerical = [field.decode() for field in constituent_numerical]
    constituent_categorical = [field.decode() for field in constituent_categorical]

    with warnings.catch_warnings():
        warnings.filterwarnings('ignore', message='found duplicate branch')
        warnings.filterwarnings('ignore', message='missing cross-reference index')
        events = NanoEventsFactory.from_root(path, schemaclass=PFNanoAODSchema).events()

    jets = events.Jet[(ak.count(events.Jet.matched_gen.pt, axis=1) >= 2)]

    sorted_jets = jets[ak.argsort(jets.matched_gen.pt, ascending=False, axis=1)]

    leading_jets = ak.concatenate((sorted_jets[:,0], sorted_jets[:,1]), axis=0)

    selected_jets = leading_jets[(leading_jets.matched_gen.pt > 30) & (abs(leading_jets.matched_gen.eta) < 2.5)]

    valid_jets = selected_jets[~ak.is_none(selected_jets.matched_gen.pt)]

    target = valid_jets.matched_gen.pt / valid_jets.pt

    globals = valid_jets[global_numerical + global_categorical]

    constituents = valid_jets.constituents.pf[constituent_numerical + constituent_categorical]
    
    data = [ak.to_numpy(target).astype(np.float32)]

    for field in global_numerical:
        data.append(ak.to_numpy(globals[field]).astype(np.float32))

    for field in global_categorical:
        data.append(ak.to_numpy(globals[field]).astype(np.int32))

    if net == 'deepset':
        row_lengths = ak.num(constituents, axis=1)
        flat_constituents = ak.flatten(constituents, axis=1)

        for field in constituent_numerical:
            data.append(ak.to_numpy(flat_constituents[field]).astype(np.float32))

        for field in constituent_categorical:
            data.append(ak.to_numpy(flat_constituents[field]).astype(np.int32))

        data.append(ak.to_numpy(row_lengths).astype(np.int32))

    if net == 'particlenet':
        for field in constituent_numerical:
            none_padded_constituent = ak.pad_none(constituents[field], target=num_points, clip=True, axis=1)
            zero_padded_constituent = ak.to_numpy(none_padded_constituent).filled(0)
            data.append(zero_padded_constituent.astype(np.float32))

        for field in constituent_categorical:
            none_padded_constituent = ak.pad_none(constituents[field], target=num_points, clip=True, axis=1)
            zero_padded_constituent = ak.to_numpy(none_padded_constituent).filled(0)
            data.append(zero_padded_constituent.astype(np.int32))

    return data
