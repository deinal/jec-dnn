import os
import argparse
import pickle
import itertools
import awkward as ak
import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
from data import read_nanoaod


def read_data(paths, predictions):
    dfs = []
    for path in paths:
        valid_jets = read_nanoaod(path)

        jet_pt = ak.to_pandas(valid_jets.pt)
        gen_jet_pt = ak.to_pandas(valid_jets.matched_gen.pt)
        parton_flavour = ak.to_pandas(valid_jets.matched_gen.partonFlavour)
        hadron_flavour = ak.to_pandas(valid_jets.matched_gen.hadronFlavour)
        qgl_axis2 = ak.to_pandas(valid_jets.qgl_axis2)

        df = pd.concat((jet_pt, gen_jet_pt, parton_flavour, hadron_flavour, qgl_axis2), axis=1)
        df.columns = ['Jet_pt', 'GenJet_pt', 'GenJet_partonFlavour', 'GenJet_hadronFlavour', 'Jet_qgl_axis2']

        flavour = df.GenJet_hadronFlavour.where(df.GenJet_hadronFlavour != 0, other=np.abs(df.GenJet_partonFlavour))
        df = df.drop(columns=['GenJet_partonFlavour', 'GenJet_hadronFlavour'])
        df['flavour'] = flavour
       
        dfs.append(df)

    df = pd.concat(dfs, axis=0)

    df['response'] = df.Jet_pt / df.GenJet_pt

    corrected_pt = predictions.flatten() * df.Jet_pt
    df['dnn_response'] = corrected_pt / df.GenJet_pt

    return df


def plot_loss(history, outdir):
    plt.plot(history['loss'][1:], label='Training loss')
    plt.plot(history['val_loss'][1:], label='Validation loss')

    changes = np.where(np.roll(history['lr'], 1) != history['lr'])[0]
    ymin, ymax = plt.gca().get_ylim()
    plt.vlines(changes[1:], ymin=ymin, ymax=ymax, ls='dashed', lw=0.8, colors='gray')

    plt.xlabel('Epoch')
    plt.ylabel('Loss')
    plt.legend()
    plt.savefig(f'{outdir}/loss.png')
    plt.close()


def plot_mean_response(df, flavour, outdir):
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


def plot_mean_residual(df, flavour_1, flavour_2, outdir):
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


def plot_distrs(dataframe, fig_dir):
    """Plot distributions of response in a few representative bins."""

    binning = np.linspace(0.5, 1.5, num=101)
    pt_bins = [(30., 40.), (100., 110.), (1000., 1100.)]

    ref_histograms, nn_histograms = {}, {}
    for ipt, pt_bin in enumerate(pt_bins):
        df_bin = dataframe[
            (dataframe.GenJet_pt >= pt_bin[0]) & (dataframe.GenJet_pt < pt_bin[1])
        ]
        for label, selection in [
            ('uds', (df_bin.flavour <= 3) & (df_bin.flavour != 0)),
            ('b', df_bin.flavour == 5),
            ('g', df_bin.flavour == 21)
        ]:
            h, _ = np.histogram(df_bin.response[selection], bins=binning)
            ref_histograms[ipt, label] = h
            h, _ = np.histogram(df_bin.dnn_response[selection], bins=binning)
            nn_histograms[ipt, label] = h

    for ipt, flavour in itertools.product(
        range(len(pt_bins)), ['uds', 'b', 'g']
    ):
        fig = plt.figure()
        axes = fig.add_subplot()
        axes.hist(
            binning[:-1], weights=ref_histograms[ipt, flavour],
            bins=binning, histtype='step', label='Standard')
        axes.hist(
            binning[:-1], weights=nn_histograms[ipt, flavour],
            bins=binning, histtype='step', label='DNN')
        axes.axvline(1., ls='dashed', lw=0.8, c='gray')
        axes.margins(x=0)
        axes.set_xlabel(
            r'$p_\mathrm{T}^\mathrm{corr}\//\/p_\mathrm{T}^\mathrm{ptcl}$')
        axes.set_ylabel('Jets')
        axes.legend()
        axes.text(
            1., 1.002,
            r'${}$, ${:g} < p_\mathrm{{T}}^\mathrm{{ptcl}} < {:g}$ GeV, '
            r'$|\eta^\mathrm{{ptcl}}| < 2.5$'.format(
                flavour, pt_bins[ipt][0], pt_bins[ipt][1]
            ),
            ha='right', va='bottom', transform=axes.transAxes
        )
        fig.savefig(os.path.join(
            fig_dir, f'{flavour}_pt{ipt + 1}.png'
        ))
        plt.close(fig)


def bootstrap(x, num=30):
    """Compute errors on median and IQR with bootstrapping."""

    if len(x) == 0:
        return np.nan, np.nan

    medians, iqrs = [], []
    for _ in range(num):
        x_resampled = np.random.choice(x, len(x))
        medians.append(np.median(x_resampled))
        quantiles = np.percentile(x_resampled, [25, 75])
        iqrs.append(quantiles[1] - quantiles[0])
    return np.std(medians), np.std(iqrs)


def compute_iqr(groups):
    """Compute IQR from series GroupBy."""
    
    q = groups.quantile([0.25, 0.75])
    iqr = q[1::2].values - q[0::2].values

    return iqr


def plot_summary(dataframe, fig_dir):
    """Plot median response and its IQR as a function of pt."""

    pt_binning = np.geomspace(20., 3e3, num=20)
    pt_centres = np.sqrt(pt_binning[:-1] * pt_binning[1:])
    for flavour_label, flavours in [
            ('uds', {1, 2, 3}), ('b', {5}), ('g', {21}),
            ('all', {0, 1, 2, 3, 4, 5, 21})
        ]:

        df_bin = dataframe[dataframe.flavour.isin(flavours)]
        bins = df_bin.groupby(pd.cut(df_bin.GenJet_pt, pt_binning))

        ref_median = bins.response.median().to_numpy()
        ref_iqr = compute_iqr(bins.response)
        ref_median_error = np.empty_like(ref_median)
        ref_iqr_error = np.empty_like(ref_median)
        for i, (_, df) in enumerate(bins):
            ref_median_error[i], ref_iqr_error[i] = bootstrap(
                df.response.to_numpy())

        nn_median = bins.dnn_response.median().to_numpy()
        nn_iqr = compute_iqr(bins.dnn_response)
        nn_median_error = np.empty_like(ref_median)
        nn_iqr_error = np.empty_like(ref_median)
        for i, (_, df) in enumerate(bins):
            nn_median_error[i], nn_iqr_error[i] = bootstrap(
                df.dnn_response.to_numpy())

        fig = plt.figure()
        axes = fig.add_subplot()
        axes.errorbar(
            pt_centres, ref_median, yerr=ref_median_error,
            marker='o', lw=0, elinewidth=0.8, label='Standard')
        axes.errorbar(
            pt_centres, nn_median, yerr=nn_median_error,
            marker='o', lw=0, elinewidth=0.8, label='DNN')
        axes.axhline(1., ls='dashed', lw=0.8, c='gray')
        axes.set_xlim(pt_binning[0], pt_binning[-1])
        axes.set_xscale('log')
        axes.set_xlabel(r'$p_\mathrm{T}^\mathrm{ptcl}$')
        axes.set_ylabel('Median response')
        axes.legend()
        axes.text(
            1., 1.002,
            '{}$|\\eta^\\mathrm{{ptcl}}| < 2.5$'.format(
                f'${flavour_label}$, ' if flavour_label != 'all' else ''
            ),
            ha='right', va='bottom', transform=axes.transAxes
        )
        fig.savefig(os.path.join(
            fig_dir, f'{flavour_label}_median.png'))
        plt.close(fig)

        fig = plt.figure()
        gs = mpl.gridspec.GridSpec(2, 1, hspace=0.02, height_ratios=[4, 1])
        axes_upper = fig.add_subplot(gs[0, 0])
        axes_lower = fig.add_subplot(gs[1, 0])

        axes_upper.errorbar(
            pt_centres, ref_iqr / ref_median, yerr=ref_iqr_error / ref_median,
            marker='o', lw=0, elinewidth=0.8, label='Standard')
        axes_upper.errorbar(
            pt_centres, nn_iqr / nn_median, yerr=nn_iqr_error / nn_median,
            marker='o', lw=0, elinewidth=0.8, label='DNN')
        axes_lower.plot(
            pt_centres, (nn_iqr / nn_median) / (ref_iqr / ref_median),
            marker='o', lw=0, color='black')

        axes_upper.set_ylim(0., None)
        axes_lower.set_ylim(0.85, 1.02)
        for axes in [axes_upper, axes_lower]:
            axes.set_xscale('log')
            axes.set_xlim(pt_binning[0], pt_binning[-1])
        axes_upper.xaxis.set_major_formatter(mpl.ticker.NullFormatter())
        axes_upper.xaxis.set_minor_formatter(mpl.ticker.NullFormatter())
        axes_upper.legend()
        axes_upper.text(
            1., 1.002,
            '{}$|\\eta^\\mathrm{{ptcl}}| < 2.5$'.format(
                f'${flavour_label}$, ' if flavour_label != 'all' else ''
            ),
            ha='right', va='bottom', transform=axes_upper.transAxes
        )
        axes_upper.set_ylabel('IQR / median for response')
        axes_lower.set_ylabel('Ratio')
        axes_lower.set_xlabel(r'$p_\mathrm{T}^\mathrm{ptcl}$')
        fig.align_ylabels()

        fig.savefig(os.path.join(
            fig_dir, f'{flavour_label}_iqr.png'))
        plt.close(fig)


def compare_flavours(dataframe, fig_dir):
    """Plot median response as a function of jet flavour."""

    pt_cut = 30
    df_pteta = dataframe[(dataframe.GenJet_pt > pt_cut)]

    ref_median, ref_median_error = [], []
    nn_median, nn_median_error = [], []
    flavours = [('g', {21}), ('uds', {1, 2, 3}), ('c', {4}), ('b', {5})]
    for _, pdg_ids in flavours:
        df = df_pteta[df_pteta.flavour.isin(pdg_ids)]
        ref_median.append(df.response.median())
        ref_median_error.append(bootstrap(df.response)[0])
        nn_median.append(df.dnn_response.median())
        nn_median_error.append(bootstrap(df.dnn_response)[0])

    fig = plt.figure()
    axes = fig.add_subplot()
    axes.errorbar(
        np.arange(len(flavours)) - 0.02, ref_median, yerr=ref_median_error,
        marker='o', ms=2, lw=0, elinewidth=0.8, label='Standard'
    )
    axes.errorbar(
        np.arange(len(flavours)) + 0.02, nn_median, yerr=nn_median_error,
        marker='o', ms=2, lw=0, elinewidth=0.8, label='DNN'
    )
    axes.set_xlim(-0.5, len(flavours) - 0.5)
    axes.axhline(1, ls='dashed', lw=0.8, c='gray')
    axes.set_xticks(np.arange(len(flavours)))
    axes.set_xticklabels([f[0] for f in flavours])
    axes.legend()
    axes.set_ylabel('Median response')
    axes.text(
        1., 1.002,
        r'$p_\mathrm{{T}}^\mathrm{{ptcl}} > {:g}$ GeV, '
        r'$|\eta^\mathrm{{ptcl}}| < 2.5$'.format(pt_cut),
        ha='right', va='bottom', transform=axes.transAxes
    )
    fig.savefig(os.path.join(fig_dir, 'flavours.png'))
    plt.close(fig)


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

    plot_loss(history, args.outdir)

    with open(f'{args.indir}/predictions.pkl', 'rb') as f:
        predictions, test_files = pickle.load(f)
    
    df = read_data(test_files, predictions)

    for subdir in ['distributions', 'summary', 'flavours', 'response', 'residual']:
        try:
            os.makedirs(os.path.join(args.outdir, subdir))
        except FileExistsError:
            pass

    plot_distrs(df, os.path.join(args.outdir, 'distributions'))
    plot_summary(df, os.path.join(args.outdir, 'summary'))
    compare_flavours(df, os.path.join(args.outdir, 'flavours'))

    for flavour, ids in [
        ('uds', {1, 2, 3}), ('b', {5}), ('g', {21}), ('all', {0, 1, 2, 3, 4, 5, 21})
    ]:
        plot_mean_response(df[df.flavour.isin(ids)], flavour, os.path.join(args.outdir, 'response'))

    plot_mean_residual(df, ('g', {21}), ('uds', {1, 2, 3}), os.path.join(args.outdir, 'residual'))