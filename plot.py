import os
import argparse
import warnings
import pickle
import awkward as ak
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from coffea.nanoevents import NanoEventsFactory, PFNanoAODSchema


def read_data(paths, predictions):
    dfs = []
    for path in paths:
        with warnings.catch_warnings():
            warnings.filterwarnings('ignore', message='Found duplicate branch')
            warnings.filterwarnings('ignore', message='Missing cross-reference index')
            events = NanoEventsFactory.from_root(path, schemaclass=PFNanoAODSchema).events()

        jets = events.Jet[(ak.count(events.Jet.matched_gen.pt, axis=1) >= 2)]

        sorted_jets = jets[ak.argsort(jets.matched_gen.pt, ascending=False, axis=1)]

        leading_jets = ak.concatenate((sorted_jets[:,0], sorted_jets[:,1]), axis=0)

        selected_jets = leading_jets[(leading_jets.matched_gen.pt > 30) & (abs(leading_jets.matched_gen.eta) < 2.5)]

        valid_jets = selected_jets[~ak.is_none(selected_jets.matched_gen.pt)]

        jet_pt = ak.to_pandas(valid_jets.pt)
        gen_jet_pt = ak.to_pandas(valid_jets.matched_gen.pt)
        parton_flavour = ak.to_pandas(valid_jets.matched_gen.partonFlavour)
        hadron_flavour = ak.to_pandas(valid_jets.matched_gen.hadronFlavour)

        df = pd.concat((jet_pt, gen_jet_pt, parton_flavour, hadron_flavour), axis=1)
        df.columns = ['Jet_pt', 'GenJet_pt', 'GenJet_partonFlavour', 'GenJet_hadronFlavour']

        flavour = df.GenJet_hadronFlavour.where(df.GenJet_hadronFlavour != 0, other=np.abs(df.GenJet_partonFlavour))
        df = df.drop(columns=['GenJet_partonFlavour', 'GenJet_hadronFlavour'])
        df['flavour'] = flavour
       
        dfs.append(df)

    df = pd.concat(dfs, axis=0)

    df['response'] = df.Jet_pt / df.GenJet_pt

    corrected_pt = predictions.flatten() * df.Jet_pt
    df['dnn_response'] = corrected_pt / df.GenJet_pt

    return df


def plot_loss(outdir, history):
    plt.plot(history['loss'][1:], label='Training loss')
    plt.plot(history['val_loss'][1:], label='Validation loss')
    plt.xlabel('Epoch')
    plt.ylabel('Loss')
    plt.legend()
    plt.savefig(f'{outdir}/loss.png')


def plot_mean_response(outdir, df, flavour):
    binning = np.geomspace(30, 3000, 20)
    bin_centers = np.sqrt(binning[:-1] * binning[1:])

    bins = df.groupby(pd.cut(df.GenJet_pt, binning))

    response_mean = bins.response.mean()
    response_std = bins.response.std()

    dnn_response_mean = bins.dnn_response.mean() 
    dnn_response_std = bins.dnn_response.std()

    fig, ax = plt.subplots(
        nrows=2, ncols=1, sharex=True, figsize=(10, 6),
        gridspec_kw={'height_ratios': [2, 1], 'hspace': 0.1}
    )
    fig.suptitle('Mean ' + flavour + '-jet response w.r.t. gen p$_{T}$')

    ax[0].errorbar(bin_centers, dnn_response_mean, yerr=dnn_response_std, ms=4, fmt='o', alpha=.7, capsize=3, capthick=1, label='DNN')
    ax[0].fill_between(bin_centers, dnn_response_mean - dnn_response_std, dnn_response_mean + dnn_response_std, alpha=.2)
    ax[0].errorbar(bin_centers, response_mean, yerr=response_std, ms=4, fmt='o', alpha=.7, capsize=3, capthick=1, label='Standard')
    ax[0].fill_between(bin_centers, response_mean - response_std, response_mean + response_std, alpha=.2)
    ax[0].axhline(1, ls='dashed', c='gray', alpha=.7)
    ax[0].set_ylabel('Jets/bin')
    ax[0].set_ylabel('Mean response')
    ax[0].legend()

    ax[1].hist(df.GenJet_pt, bins=binning, alpha=.7)
    ax[1].set_xscale('log')
    ax[1].set_ylabel('Jets/bin')
    ax[1].set_xlabel('gen p$_{T}$')

    fig.savefig(f'{outdir}/{flavour}_mean_response.png')


def plot_mean_residual(outdir, df, flavour_1, flavour_2):
    binning = np.geomspace(30, 3000, 20)
    bin_centers = np.sqrt(binning[:-1] * binning[1:])

    df1 = df[df.flavour.isin(flavour_1[1])]
    bins_1 = df1.groupby(pd.cut(df1.GenJet_pt, binning))
    response_mean_1 = bins_1.response.mean()
    response_std_1 = bins_1.response.std()
    dnn_response_mean_1 = bins_1.dnn_response.mean() 
    dnn_response_std_1 = bins_1.dnn_response.std()

    df2 = df[df.flavour.isin(flavour_2[1])]
    bins_2 = df2.groupby(pd.cut(df2.GenJet_pt, binning))
    response_mean_2 = bins_2.response.mean()
    response_std_2 = bins_2.response.std()
    dnn_response_mean_2 = bins_2.dnn_response.mean() 
    dnn_response_std_2 = bins_2.dnn_response.std()

    
    dnn_response_mean_diff = dnn_response_mean_1 - dnn_response_mean_2
    dnn_err = np.sqrt(dnn_response_std_1 ** 2 + dnn_response_std_2 ** 2)
    response_mean_diff = response_mean_1 - response_mean_2
    err = np.sqrt(response_std_1 ** 2 + response_std_2 ** 2)

    fig, ax = plt.subplots(
        nrows=2, ncols=1, sharex=True, figsize=(10, 6),
        gridspec_kw={'height_ratios': [2, 1], 'hspace': 0.1}
    )
    fig.suptitle('Mean response residuals w.r.t. gen p$_{T}$')

    ax[0].errorbar(bin_centers, dnn_response_mean_diff, yerr=dnn_err, ms=4, fmt='o', alpha=.7, capsize=3, capthick=1, label='DNN')
    ax[0].fill_between(bin_centers, dnn_response_mean_diff - dnn_err, dnn_response_mean_diff + dnn_err, alpha=.2)
    ax[0].errorbar(bin_centers, response_mean_diff, yerr=err, ms=4, fmt='o', alpha=.7, capsize=3, capthick=1, label='Standard')
    ax[0].fill_between(bin_centers, response_mean_diff - err, response_mean_diff + err, alpha=.2)
    ax[0].axhline(0, ls='dashed', c='gray', alpha=.7)
    ax[0].set_ylabel('Jets/bin')
    ax[0].set_ylabel('R$_{' + flavour_1[0] + '}$-R$_{' + flavour_2[0] + '}$')
    ax[0].legend()

    ax[1].hist(df.GenJet_pt, bins=binning, alpha=.7)
    ax[1].set_xscale('log')
    ax[1].set_ylabel('Jets/bin')
    ax[1].set_xlabel('gen p$_{T}$')

    fig.savefig(f'{outdir}/{flavour_1[0]}_{flavour_2[0]}_mean_residual.png')


if __name__ == '__main__':
    arg_parser = argparse.ArgumentParser(description=__doc__)
    arg_parser.add_argument('-i', '--indir', required=True, help='Directory with dnn training results')
    arg_parser.add_argument('-o', '--outdir', required=True, help='Where to store plots')
    args = arg_parser.parse_args()

    try:
        os.mkdir(f'{args.outdir}')
    except FileExistsError:
        pass

    with open(f'{args.indir}/history.pkl', 'rb') as f:
        history = pickle.load(f)

    plot_loss(args.outdir, history)

    with open(f'{args.indir}/predictions.pkl', 'rb') as f:
        predictions, test_files = pickle.load(f)
    
    df = read_data(test_files, predictions)

    for flavour, ids in [
        ('uds', {1, 2, 3}), ('b', {5}), ('g', {21}), ('all', {0, 1, 2, 3, 4, 5, 21})
    ]:
        plot_mean_response(args.outdir, df[df.flavour.isin(ids)], flavour)

    plot_mean_residual(args.outdir, df, ('g', {21}), ('uds', {1, 2, 3}))