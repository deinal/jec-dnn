import argparse
import pickle
import yaml
import os
import tensorflow as tf
from src.data import create_datasets


if __name__ == '__main__':
    arg_parser = argparse.ArgumentParser(description=__doc__)
    arg_parser.add_argument('-i', '--indir', required=True, help='Directory containing training data in nanoaod root format')
    arg_parser.add_argument('-o', '--outdir', required=True, help='Directory to save tensorflow datasets')
    arg_parser.add_argument('-c', '--config', required=True, help='Config file')
    arg_parser.add_argument('--zip', action='store_true', help='Whether or not to zip the output')
    args = arg_parser.parse_args()

    with open(args.config) as f:
        config = yaml.safe_load(f)
        net = config['net']

    train_ds, val_ds, test_ds, metadata = create_datasets(net, args.indir, config['data'])

    if args.zip:
        compression = 'GZIP'
    else:
        compression = None

    tf.data.experimental.save(train_ds, path=os.path.join(args.outdir, 'train'), compression=compression)
    tf.data.experimental.save(val_ds, path=os.path.join(args.outdir, 'val'), compression=compression)
    tf.data.experimental.save(test_ds, path=os.path.join(args.outdir, 'test'), compression=compression)

    with open(os.path.join(args.outdir, 'metadata.pkl'), 'wb') as f:
        pickle.dump(metadata, f)
