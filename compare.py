import os
import argparse
import pickle
import itertools
import awkward as ak
import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
from src.data import read_nanoaod


def read_data(paths, ds_predictions, pn_predictions):
    dfs = []
    for path in paths:
        valid_jets = read_nanoaod(path)

        jet_pt = ak.to_pandas(valid_jets.pt)
        gen_jet_pt = ak.to_pandas(valid_jets.matched_gen.pt)
        gen_jet_eta = ak.to_pandas(valid_jets.matched_gen.eta)
        parton_flavour = ak.to_pandas(valid_jets.matched_gen.partonFlavour)
        hadron_flavour = ak.to_pandas(valid_jets.matched_gen.hadronFlavour)

        df = pd.concat((jet_pt, gen_jet_pt, gen_jet_eta, parton_flavour, hadron_flavour), axis=1)
        df.columns = ['Jet_pt', 'GenJet_pt', 'GenJet_eta', 'GenJet_partonFlavour', 'GenJet_hadronFlavour']

        flavour = df.GenJet_hadronFlavour.where(df.GenJet_hadronFlavour != 0, other=np.abs(df.GenJet_partonFlavour))
        df = df.drop(columns=['GenJet_partonFlavour', 'GenJet_hadronFlavour'])
        df['flavour'] = flavour
       
        dfs.append(df)

    df = pd.concat(dfs, axis=0)

    df['response'] = df.Jet_pt / df.GenJet_pt
    df['ds_response'] = ds_predictions.flatten() * df.Jet_pt / df.GenJet_pt
    df['pn_response'] = pn_predictions.flatten() * df.Jet_pt / df.GenJet_pt

    return df


def plot_distrs(dataframe, fig_dir):
    """Plot distributions of response in a few representative bins."""

    binning = np.linspace(0.5, 1.5, num=101)
    pt_bins = [(30, 40), (100, 110), (1000, 1100)]
    eta_bins = [(0., 2.5), (2.5, 5)]

    ref_histograms, ds_histograms, pn_histograms = {}, {}, {}
    for (ipt, pt_bin), (ieta, eta_bin) in itertools.product(
        enumerate(pt_bins), enumerate(eta_bins)
    ):
        df_bin = dataframe[
            (dataframe.GenJet_pt >= pt_bin[0]) & (dataframe.GenJet_pt < pt_bin[1])
            & (np.abs(dataframe.GenJet_eta) >= eta_bin[0])
            & (np.abs(dataframe.GenJet_eta) < eta_bin[1])
        ]
        for label, selection in [
            ('uds', (df_bin.flavour <= 3) & (df_bin.flavour != 0)),
            ('b', df_bin.flavour == 5),
            ('g', df_bin.flavour == 21)
        ]:
            h, _ = np.histogram(df_bin.response[selection], bins=binning)
            ref_histograms[ipt, ieta, label] = h
            h, _ = np.histogram(df_bin.ds_response[selection], bins=binning)
            ds_histograms[ipt, ieta, label] = h
            h, _ = np.histogram(df_bin.pn_response[selection], bins=binning)
            pn_histograms[ipt, ieta, label] = h

    for ipt, ieta, flavour in itertools.product(
        range(len(pt_bins)), range(len(eta_bins)), ['uds', 'b', 'g']
    ):
        fig = plt.figure()
        ax = fig.add_subplot()
        ax.hist(
            binning[:-1], weights=ref_histograms[ipt, ieta, flavour],
            bins=binning, histtype='step', label='Standard')
        ax.hist(
            binning[:-1], weights=ds_histograms[ipt, ieta, flavour],
            bins=binning, histtype='step', label='Deep Sets')
        ax.hist(
            binning[:-1], weights=pn_histograms[ipt, ieta, flavour],
            bins=binning, histtype='step', label='ParticleNet')
        ax.axvline(1., ls='dashed', lw=0.8, c='gray')
        ax.margins(x=0)
        ax.set_xlabel(
            r'$p_\mathrm{T}^\mathrm{corr}\//\/p_\mathrm{T}^\mathrm{gen}$')
        ax.set_ylabel('Jets')
        ax.legend(loc='upper right')
        ax.text(
            1., 1.002,
            r'${}$, ${:g} < p_\mathrm{{T}}^\mathrm{{gen}} < {:g}$ GeV, '
            r'${:g} < |\eta^\mathrm{{gen}}| < {:g}$'.format(
                flavour, pt_bins[ipt][0], pt_bins[ipt][1],
                eta_bins[ieta][0], eta_bins[ieta][1]
            ),
            ha='right', va='bottom', transform=ax.transAxes
        )
        ax.tick_params(
            axis='both', which='both', direction='in', 
            bottom=True, top=True, left=True, right=True
        )

        fig.savefig(os.path.join(
            fig_dir, f'{flavour}_pt{ipt + 1}_eta{ieta + 1}.png'
        ))
        plt.close(fig)


def bootstrap_median(x, num=30):
    """Compute errors on median with bootstrapping."""

    if len(x) == 0:
        return np.nan

    medians = []
    for _ in range(num):
        x_resampled = np.random.choice(x, len(x))
        medians.append(np.median(x_resampled))
    return np.std(medians)


def compare_flavours(dataframe, fig_dir):
    """Plot median response as a function of jet flavour."""
    
    pt_cut = 30
    for ieta, eta_bin in enumerate([(0, 2.5), (2.5, 5)], start=1):
        df_pteta = dataframe[
            (np.abs(dataframe.GenJet_eta) >= eta_bin[0])
            & (np.abs(dataframe.GenJet_eta) < eta_bin[1])
        ]
        ref_median, ref_median_error = [], []
        ds_median, ds_median_error = [], []
        pn_median, pn_median_error = [], []
        flavours = [('g', {21}), ('uds', {1, 2, 3}), ('c', {4}), ('b', {5})]
        for _, pdg_ids in flavours:
            df = df_pteta[df_pteta.flavour.isin(pdg_ids)]
            ref_median.append(df.response.median())
            ref_median_error.append(bootstrap_median(df.response))
            ds_median.append(df.ds_response.median())
            ds_median_error.append(bootstrap_median(df.ds_response))
            pn_median.append(df.pn_response.median())
            pn_median_error.append(bootstrap_median(df.pn_response))

        fig = plt.figure()
        ax = fig.add_subplot()
        ax.errorbar(
            np.arange(len(flavours)) - 0.04, ref_median, yerr=ref_median_error,
            marker='o', ms=3, lw=0, elinewidth=0.8, label='Standard'
        )
        ax.errorbar(
            np.arange(len(flavours)), ds_median, yerr=ds_median_error,
            marker='^', ms=3, lw=0, elinewidth=0.8, label='Deep Sets'
        )
        ax.errorbar(
            np.arange(len(flavours)) + 0.04, pn_median, yerr=pn_median_error,
            marker='s', ms=3, lw=0, elinewidth=0.8, label='ParticleNet'
        )
        ax.set_xlim(-0.5, len(flavours) - 0.5)
        ax.axhline(1, ls='dashed', lw=0.8, c='gray')
        ax.set_xticks(np.arange(len(flavours)))
        ax.set_xticklabels([f[0] for f in flavours])
        ax.legend()
        ax.set_ylabel('Median response')
        ax.text(
            1., 1.002,
            r'$p_\mathrm{{T}}^\mathrm{{gen}} > {:g}$ GeV, '
            r'${:g} < |\eta^\mathrm{{gen}}| < {:g}$'.format(
                pt_cut, eta_bin[0], eta_bin[1]
            ),
            ha='right', va='bottom', transform=ax.transAxes
        )
        ax.tick_params(
            axis='both', which='both', direction='in', 
            bottom=True, top=True, left=True, right=True
        )

        fig.savefig(os.path.join(fig_dir, f'eta{ieta}.png'))
        plt.close(fig)


def plot_median_response(outdir, flavour_label, bins, bin_centers, eta_bin, ieta):
    """Plot median response as a function of pt."""

    ref_median = bins.response.median().to_numpy()
    ref_median_error = np.empty_like(ref_median)
    for i, (_, df) in enumerate(bins):
        ref_median_error[i] = bootstrap_median(df.response.to_numpy())

    ds_median = bins.ds_response.median().to_numpy()
    ds_median_error = np.empty_like(ref_median)
    for i, (_, df) in enumerate(bins):
        ds_median_error[i] = bootstrap_median(df.ds_response.to_numpy())

    pn_median = bins.pn_response.median().to_numpy()
    pn_median_error = np.empty_like(ref_median)
    for i, (_, df) in enumerate(bins):
        pn_median_error[i] = bootstrap_median(df.pn_response.to_numpy())

    fig = plt.figure()
    fig.suptitle('Median ' + flavour_label + '-jet response w.r.t. gen p$_{T}$')
    
    ax = fig.add_subplot()

    vals = np.geomspace(0.5, 50, 20)
    shift = np.sqrt(vals[:-1] * vals[1:])
    
    ax.errorbar(
        bin_centers - shift, ref_median, yerr=ref_median_error,
        ms=3, fmt='o', elinewidth=0.8, label='Standard'
    )
    ax.errorbar(
        bin_centers, ds_median, yerr=ds_median_error, 
        ms=3, fmt='^', elinewidth=0.8, label='Deep Sets'
    )
    ax.errorbar(
        bin_centers + shift, pn_median, yerr=pn_median_error, 
        ms=3, fmt='s', elinewidth=0.8, label='ParticleNet'
    )
    ax.axhline(1, ls='dashed', c='gray', alpha=.7)
    ax.set_xlabel('$p^\\mathrm{{gen}}_{T}$')
    ax.set_ylabel('Median response')
    ax.text(
        1., 1.002,
        '{}${:g} < |\\eta^\\mathrm{{gen}}| < {:g}$'.format(
            f'${flavour_label}$, ' if flavour_label != 'all' else '',
            eta_bin[0], eta_bin[1]
        ),
        ha='right', va='bottom', transform=ax.transAxes
    )
    ax.legend(loc='upper right')
    ax.tick_params(
        axis='both', which='both', direction='in', 
        bottom=True, top=True, left=True, right=True
    )
    ax.set_xscale('log')

    fig.savefig(os.path.join(outdir, f'{flavour_label}_eta{ieta}.png'))
    plt.close(fig)


def bootstrap_iqr(x, num=30):
    """Compute errors on IQR with bootstrapping."""

    if len(x) == 0:
        return np.nan

    iqrs = []
    for _ in range(num):
        x_resampled = np.random.choice(x, len(x))
        quantiles = np.percentile(x_resampled, [25, 75])
        iqrs.append(quantiles[1] - quantiles[0])
    return np.std(iqrs)


def compute_iqr(groups):
    """Compute IQR from series GroupBy."""
    
    q = groups.quantile([0.25, 0.75])
    iqr = q[1::2].values - q[0::2].values

    return iqr


def plot_resolution(outdir, flavour_label, bins, bin_centers, eta_bin, ieta):
    ref_median = bins.response.median().to_numpy()
    ref_iqr = compute_iqr(bins.response)
    ref_iqr_error = np.empty_like(ref_iqr)
    for i, (_, df) in enumerate(bins):
        ref_iqr_error[i] = bootstrap_iqr(df.response.to_numpy())

    ds_median = bins.ds_response.median().to_numpy()
    ds_iqr = compute_iqr(bins.ds_response)
    ds_iqr_error = np.empty_like(ref_iqr)
    for i, (_, df) in enumerate(bins):
        ds_iqr_error[i] = bootstrap_iqr(df.ds_response.to_numpy())

    pn_median = bins.pn_response.median().to_numpy()
    pn_iqr = compute_iqr(bins.pn_response)
    pn_iqr_error = np.empty_like(ref_iqr)
    for i, (_, df) in enumerate(bins):
        pn_iqr_error[i] = bootstrap_iqr(df.pn_response.to_numpy())

    fig = plt.figure()
    gs = mpl.gridspec.GridSpec(2, 1, hspace=0.02, height_ratios=[4, 1])
    axes_upper = fig.add_subplot(gs[0, 0])
    axes_lower = fig.add_subplot(gs[1, 0])

    axes_upper.errorbar(
        bin_centers, ref_iqr / ref_median, yerr=ref_iqr_error / ref_median,
        ms=3, marker='o', lw=0, elinewidth=0.8, label='Standard'
    )
    axes_upper.errorbar(
        bin_centers, ds_iqr / ds_median, yerr=ds_iqr_error / ds_median,
        ms=3, marker='^', lw=0, elinewidth=0.8, label='Deep Sets'
    )
    axes_upper.errorbar(
        bin_centers, pn_iqr / pn_median, yerr=pn_iqr_error / pn_median,
        ms=3, marker='s', lw=0, elinewidth=0.8, label='ParticleNet'
    )
    axes_lower.plot(
        bin_centers, (ds_iqr / ds_median) / (ref_iqr / ref_median),
        ms=3, marker='^', lw=0, color='tab:orange'
    )
    axes_lower.plot(
        bin_centers, (pn_iqr / pn_median) / (ref_iqr / ref_median),
        ms=3, marker='s', lw=0, color='tab:green'
    )

    axes_upper.set_ylim(0, None)
    axes_lower.set_ylim(0.85, 1.02)
    for axes in [axes_upper, axes_lower]:
        axes.set_xscale('log')
        axes.set_xlim(binning[0], binning[-1])
    axes_upper.xaxis.set_major_formatter(mpl.ticker.NullFormatter())
    axes_upper.xaxis.set_minor_formatter(mpl.ticker.NullFormatter())
    axes_upper.legend()
    axes_upper.text(
        1, 1.002,
        '{}${:g} < |\\eta^\\mathrm{{gen}}| < {:g}$'.format(
            f'${flavour_label}$, ' if flavour_label != 'all' else '',
            eta_bin[0], eta_bin[1]
        ),
        ha='right', va='bottom', transform=axes_upper.transAxes
    )
    axes_upper.set_ylabel('IQR / median for response')
    axes_lower.set_ylabel('Ratio')
    axes_lower.set_xlabel(r'$p_\mathrm{T}^\mathrm{gen}$')
    axes_upper.tick_params(
        axis='both', which='both', direction='in', 
        bottom=True, top=True, left=True, right=True
    )
    axes_lower.tick_params(
        axis='both', which='both', direction='in', 
        bottom=True, top=True, left=True, right=True
    )
    fig.align_ylabels()

    fig.savefig(os.path.join(outdir, f'{flavour_label}_eta{ieta}_iqr.png'))
    plt.close(fig)


def plot_median_residual(outdir, bin_centers, flavour_labels, bins, eta_bin, ieta):
    """Plot difference in median response between flavours as a function of pt."""

    ref_median_1 = bins[0].response.median().to_numpy()
    ref_median_error_1 = np.empty_like(ref_median_1)
    for i, (_, df) in enumerate(bins[0]):
        ref_median_error_1[i] = bootstrap_median(df.response.to_numpy())

    ds_median_1 = bins[0].ds_response.median().to_numpy()
    ds_median_error_1 = np.empty_like(ref_median_1)
    for i, (_, df) in enumerate(bins[0]):
        ds_median_error_1[i] = bootstrap_median(df.ds_response.to_numpy())

    pn_median_1 = bins[0].pn_response.median().to_numpy()
    pn_median_error_1 = np.empty_like(ref_median_1)
    for i, (_, df) in enumerate(bins[0]):
        pn_median_error_1[i] = bootstrap_median(df.pn_response.to_numpy())
    
    ref_median_2 = bins[1].response.median().to_numpy()
    ref_median_error_2 = np.empty_like(ref_median_2)
    for i, (_, df) in enumerate(bins[1]):
        ref_median_error_2[i] = bootstrap_median(df.response.to_numpy())

    ds_median_2 = bins[1].ds_response.median().to_numpy()
    ds_median_error_2 = np.empty_like(ref_median_2)
    for i, (_, df) in enumerate(bins[1]):
        ds_median_error_2[i] = bootstrap_median(df.ds_response.to_numpy())

    pn_median_2 = bins[1].pn_response.median().to_numpy()
    pn_median_error_2 = np.empty_like(ref_median_2)
    for i, (_, df) in enumerate(bins[1]):
        pn_median_error_2[i] = bootstrap_median(df.pn_response.to_numpy())

    diff = ref_median_1 - ref_median_2
    err = np.sqrt(ref_median_error_1 ** 2 + ref_median_error_2 ** 2)
    ds_diff = ds_median_1 - ds_median_2
    ds_err = np.sqrt(ds_median_error_1 ** 2 + ds_median_error_2 ** 2)
    pn_diff = pn_median_1 - pn_median_2
    pn_err = np.sqrt(pn_median_error_1 ** 2 + pn_median_error_2 ** 2)

    fig = plt.figure()
    ax = fig.add_subplot()

    fig.suptitle('Median response residuals w.r.t. gen p$_{T}$')

    vals = np.geomspace(0.5, 50, 20)
    shift = np.sqrt(vals[:-1] * vals[1:])

    ax.errorbar(
        bin_centers - shift, diff, yerr=err, 
        ms=3, fmt='o', elinewidth=0.8, label='Standard'
    )
    ax.errorbar(
        bin_centers, ds_diff, yerr=ds_err, 
        ms=3, fmt='^', elinewidth=0.8, label='Deep Sets'
    )
    ax.errorbar(
        bin_centers + shift, pn_diff, yerr=pn_err, 
        ms=3, fmt='s', lw=0, elinewidth=0.8, label='ParticleNet'
    )
    ax.axhline(0, ls='dashed', c='gray', alpha=.7)
    ax.set_xlabel('$p^\\mathrm{{gen}}_{T}$')
    ax.set_ylabel('$R_{' + flavour_labels[0] + '}-R_{' + flavour_labels[1] + '}$')
    ax.text(
        1., 1.002,
        '${:g} < |\\eta^\\mathrm{{gen}}| < {:g}$'.format(eta_bin[0], eta_bin[1]),
        ha='right', va='bottom', transform=ax.transAxes
    )
    ax.set_xscale('log')
    ax.legend(loc='upper right')
    ax.tick_params(
        axis='both', which='both', direction='in', 
        bottom=True, top=True, left=True, right=True
    )

    fig.savefig(os.path.join(outdir, f'{flavour_labels[0]}-{flavour_labels[1]}_eta{ieta}.png'))
    plt.close(fig)


if __name__ == '__main__':
    arg_parser = argparse.ArgumentParser(description=__doc__)
    arg_parser.add_argument('-d', '--deepsets', required=True, help='Deep Sets results directory')
    arg_parser.add_argument('-p', '--particlenet', required=True, help='ParticleNet results directory')
    arg_parser.add_argument('-o', '--outdir', required=True, help='Where to store plots')
    args = arg_parser.parse_args()

    try:
        os.mkdir(args.outdir)
    except FileExistsError:
        pass

    with open(os.path.join(args.deepsets, 'predictions.pkl'), 'rb') as f:
        ds_predictions, ds_test_files = pickle.load(f)

    with open(os.path.join(args.particlenet, 'predictions.pkl'), 'rb') as f:
        pn_predictions, pn_test_files = pickle.load(f)

    if pn_test_files != ds_test_files:
        raise RuntimeError('Test files are different.')

    df = read_data(ds_test_files, ds_predictions, pn_predictions)

    for subdir in ['distributions', 'flavours', 'response', 'resolution', 'residual']:
        try:
            os.makedirs(os.path.join(args.outdir, subdir))
        except FileExistsError:
            pass

    plot_distrs(df, os.path.join(args.outdir, 'distributions'))
    compare_flavours(df, os.path.join(args.outdir, 'flavours'))

    binning = np.geomspace(30, 3000, 20)
    bin_centers = np.sqrt(binning[:-1] * binning[1:])

    for (ieta, eta_bin), (flavour_label, flavour_ids) in itertools.product(
        enumerate([(0, 2.5), (2.5, 5)], start=1),
        [
            ('uds', {1, 2, 3}), ('b', {5}), ('g', {21}),
            ('all', {0, 1, 2, 3, 4, 5, 21})
        ]
    ):
        df_bin = df[
            (np.abs(df.GenJet_eta) >= eta_bin[0])
            & (np.abs(df.GenJet_eta) < eta_bin[1])
            & df.flavour.isin(flavour_ids)
        ]
        bins = df_bin.groupby(pd.cut(df_bin.GenJet_pt, binning))

        plot_median_response(
            os.path.join(args.outdir, 'response'),
            flavour_label, bins, bin_centers, eta_bin, ieta
        )

        plot_resolution(
            os.path.join(args.outdir, 'resolution'),
            flavour_label, bins, bin_centers, eta_bin, ieta
        )
    
    for (ieta, eta_bin), flavours in itertools.product(
        enumerate([(0, 2.5), (2.5, 5)], start=1),
        itertools.combinations([('uds', {1, 2, 3}), ('b', {5}), ('g', {21})], r=2),
    ):
        bins = []
        for i, flavour_ids in enumerate([flavours[0][1], flavours[1][1]]):
            df_bin = df[
                (np.abs(df.GenJet_eta) >= eta_bin[0])
                & (np.abs(df.GenJet_eta) < eta_bin[1])
                & df.flavour.isin(flavour_ids)
            ]
            bins.append(df_bin.groupby(pd.cut(df_bin.GenJet_pt, binning)))

        plot_median_residual(
            os.path.join(args.outdir, 'residual'),
            bin_centers, (flavours[0][0], flavours[1][0]), bins, eta_bin, ieta
        )
