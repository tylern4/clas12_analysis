import numpy as np
import pandas as pd
from lmfit import Model
from lmfit.models import *
from scipy.special import erfc
import boost_histogram as bh
from pyarrow import csv
import pyarrow as pa
pa.set_cpu_count(8)

ENERGY = 4.81726
EK = 5.499
E16 = 5.75

Q_FULL = 4348.46636E-6  # 03/04/2021
Q_EMPTY = 1211.06E-6

overlapSettings = {"E99-107": {"name": "",
                               "color": 'r',
                               "symbol": '*',
                               "energy": E16},
                   "K. Park 2014": {"name": "",
                                    "color": 'g',
                                    "symbol": 'd',
                                    "energy": EK, },
                   }

which_plot = {
    -1.0: [0, 0],
    -0.8: [0, 1],
    -0.6: [1, 0],
    -0.4: [1, 1],
    -0.2: [2, 0],
    0.0: [2, 1],
    0.2: [3, 0],
    0.4: [3, 1],
    0.6: [4, 0],
    0.8: [4, 1]
}
plot_label = {
    -1.0: "$\cos(\\theta)=[-1.0,-0.8)$",
    -0.8: "$\cos(\\theta)=[-0.8,-0.6)$",
    -0.6: "$\cos(\\theta)=[-0.6,-0.4)$",
    -0.4: "$\cos(\\theta)=[-0.4,-0.2)$",
    -0.2: "$\cos(\\theta)=[-0.2,0.0)$",
    0.0: "$\cos(\\theta)=[0.0,0.2)$",
    0.2: "$\cos(\\theta)=[0.2,0.4)$",
    0.4: "$\cos(\\theta)=[0.4,0.6)$",
    0.6: "$\cos(\\theta)=[0.6,0.8)$",
    0.8: "$\cos(\\theta)=[0.8,1.0)$"
}

# Two sets of W and Q2 bins
w_bins_e99 = np.array([1.1, 1.12, 1.14, 1.16, 1.18, 1.2, 1.22, 1.24, 1.26, 1.28, 1.3,
                       1.32, 1.34, 1.36, 1.38, 1.4, 1.42, 1.44, 1.46, 1.48, 1.5, 1.52,
                       1.54, 1.56, 1.58, 1.6, 1.62, 1.64, 1.66, 1.68, 1.7, 1.72, 1.74, 1.76, 1.78, 1.8, 1.82])

w_bins_k = np.array([1.605, 1.615, 1.625, 1.635, 1.645, 1.655, 1.665, 1.675, 1.685, 1.695,
                     1.705, 1.715, 1.725, 1.735, 1.745, 1.755, 1.765, 1.775, 1.785, 1.795, 1.805])

q2_bins_e99 = np.array([1.1, 1.33, 1.56, 1.87, 2.23, 2.66, 3.5])
q2_bins_k = np.array([1.8,  2.2,  2.6,  3.15, 4.0])

theta_bins = np.array([-1.0, -0.8, -0.6, -0.4, -0.2,
                       0.0, 0.2, 0.4, 0.6, 0.8, 1.0])


def read_csv(file_name: str = None, data: bool = False):
    names = [
        "electron_sector",
        "w",
        "q2",
        "pip_theta",
        "pip_phi",
        "mm2"
    ]
    dtype = {
        "electron_sector": "int8",
        "w": "float32",
        "q2": "float32",
        "pip_theta": "float32",
        "pip_phi": "float32",
        "mm2": "float32",
    }

    # Load file into pyTable before changing to pandas
    pyTable = csv.read_csv(
        file_name,
        read_options=csv.ReadOptions(
            use_threads=True, column_names=names, skip_rows=1),
        convert_options=csv.ConvertOptions(column_types=dtype),
    )
    df = pyTable.to_pandas(strings_to_categorical=True)

    if data:
        return df.dropna()

    mc_rec = df[df.type == "mc_rec"].dropna()
    thrown = df[df.type == "thrown"].dropna()
    del df

    return (
        mc_rec,
        thrown,
    )


def luminosity():
    l = 5  # cm
    rho = 0.0708  # g/cm3
    Avigadro = 6.022E23  # mol^−1
    qe = 1.602E-19  # C
    MH = 1.007  # g/mol
    conv_cm2_to_mubarn = 1E-30  # From wolfram alpha

    return (l*rho*Avigadro)/(qe*MH) * conv_cm2_to_mubarn


def momentum_calc(energy, mass):
    return np.sqrt(energy**2 - mass**2)


def virtual_photon_energy_calc(target_mass, w, q2):
    return (((w**2 + q2) / target_mass) - target_mass) / 2


def costheta_e_calc(beam_energy, target_mass, w, q2):
    electron_mass = 5.109989433549345e-4
    nu = virtual_photon_energy_calc(target_mass, w, q2)
    scattered_energy = ((beam_energy) - (nu))
    beam_momentum = momentum_calc(beam_energy, electron_mass)
    scattered_momentum = momentum_calc(scattered_energy, electron_mass)
    return (beam_energy * scattered_energy - q2 / 2.0 - electron_mass**2) / (beam_momentum * scattered_momentum)


def virtual_photon_epsilon_calc(beam_energy, w, q2, target_mass: float = 0.93827203):
    theta_e = np.arccos(costheta_e_calc(beam_energy, target_mass, w, q2))
    nu = virtual_photon_energy_calc(target_mass, w, q2)
    return np.power(1 + 2 * (1 + nu**2 / q2) * np.power(np.tan(theta_e / 2), 2), -1)


def virtual_photon_flux(w: float, q2: float, beam_energy: float = 4.81726, target_mass: float = 0.93827203) -> float:
    alpha = 0.007297352570866302
    epsilon = virtual_photon_epsilon_calc(beam_energy, w, q2, target_mass)
    return alpha / (4 * np.pi * q2) * w / (beam_energy**2 * target_mass**2) * (w**2 - target_mass**2) / (1 - epsilon)


def degauss(x, A, mu, sigma, lambda1, lambda2):
    mu1 = sigma * sigma * lambda1 + x - mu
    mu2 = -sigma * sigma * lambda2 + x - mu
    ret = A * 0.5 / (1.0 / lambda1 + 1.0 / lambda2) * \
        (np.exp(0.5 * np.power(sigma * lambda1, 2) + lambda1 * (x - mu)) * erfc(mu1 / (sigma * np.sqrt(2.0)))
         + np.exp(0.5 * np.power(sigma * lambda2, 2) - lambda2 * (x - mu)) * erfc(-mu2 / (sigma * np.sqrt(2.0))))

    return ret


def gauss(x, A, mu, sig):
    ret = np.exp(-np.power(x - mu, 2.) / (2 * np.power(sig, 2.)))
    return A*ret


def peak(x, c):
    return np.exp(-np.power(x - c, 2) / 16.0)


def lin_interp(x, y, i, half):
    return x[i] + (x[i+1] - x[i]) * ((half - y[i]) / (y[i+1] - y[i]))


def half_max_x(x, y):
    half = np.max(y)/2.0
    signs = np.sign(np.add(y, -half))
    zero_crossings = (signs[0:-2] != signs[1:-1])
    zero_crossings_i = np.where(zero_crossings)[0]
    return [lin_interp(x, y, zero_crossings_i[0], half),
            lin_interp(x, y, zero_crossings_i[1], half)]


def mm_cut(df: pd.DataFrame, sigma: int = 3):
    data = {}

    for sec in range(1, 7):
        if np.sum(df.electron_sector == sec) == 0:
            continue
        y, x = bh.numpy.histogram(
            df[df.electron_sector == sec].mm2, bins=150, density=True
        )
        x = (x[1:] + x[:-1]) / 2

        peak = PseudoVoigtModel(prefix="peak_")
        pars = peak.guess(y, x=x)
        background = GaussianModel(prefix="back_")
        pars.update(background.make_params())
        model = peak * background
        out = model.fit(y, pars, x=x)
        xs = np.linspace(0.3, 1.5, 1000)
        min_cut = out.params['peak_center'] - \
            sigma*out.params['peak_fwhm'] / 2.355
        max_cut = out.params['peak_center'] + \
            sigma*out.params['peak_fwhm'] / 2.355
        data[sec] = (min_cut, max_cut if max_cut <= 1.0 else 1.0)

    return data


def cut_for_MM(rec):
    sector_cuts = mm_cut(rec, sigma=3)

    cuts = False
    mc_cuts = False
    empty_cuts = False

    for sec, min_max in sector_cuts.items():
        cuts |= (
            (rec.electron_sector == sec)
            & (rec.mm2 >= min_max[0])
            & (rec.mm2 <= min_max[1])
        )

    rec = rec[cuts]

    return rec


def hist_data(data, density=True, bins=10):
    data_y, data_x = bh.numpy.histogram(
        data.pip_phi.to_numpy(), bins=bins, range=(0, 2 * np.pi), density=density, threads=4)
    x = (data_x[1:] + data_x[:-1]) / 2.0
    return data_y, x


def prep_for_ana(dataframe, w_bins, q2_bins, theta_bins):
    dataframe.dropna(inplace=True)
    dataframe["cos_theta"] = np.cos(dataframe.pip_theta).astype(np.float32)
    dataframe["w_bin"] = pd.cut(
        dataframe["w"], bins=w_bins, include_lowest=False)
    dataframe["q2_bin"] = pd.cut(
        dataframe["q2"], bins=q2_bins, include_lowest=False)
    dataframe["theta_bin"] = pd.cut(
        dataframe["cos_theta"], bins=theta_bins, include_lowest=False)
    dataframe.dropna(inplace=True)

    return dataframe


def make_cuts(dataframe, w, q2, theta):
    df1 = dataframe[dataframe.w_bin == w].copy()
    df2 = df1[df1.q2_bin == q2].copy()
    del df1
    df3 = df2[df2.theta_bin == theta].copy()
    del df2
    return df3


def statistical(DN_full, DN_empty, denom):

    error = DN_full/(Q_FULL**2) + DN_empty/(Q_EMPTY**2)
    error = np.sqrt(error) / denom

    return error


def get_error_bars(y, mc_rec_y, thrown_y, stat_error):
    F = mc_rec_y/thrown_y

    dF = ((thrown_y - mc_rec_y) * mc_rec_y) / (thrown_y**3)

    error = np.sqrt(dF)/F

    error_bar = np.sqrt(error**2 + stat_error**2)

    return error_bar


@np.vectorize
def isclose(a, b, rel_tol=1e-4, abs_tol=0.0):
    return np.abs(a-b) <= np.maximum(rel_tol * np.maximum(np.abs(a), np.abs(b)), abs_tol)


def A(M, B, C):
    if (C > 0 and np.abs(B) <= 4*C):
        return M**2 + B**2/(8*C) + C
    else:
        return M**2 + np.abs(B) - C


def model_new(x, M, b, c):
    """
    a => sigma_l + sigma_t
    b => epsilon*sigma_tt
    c => Sqrt(2epsilon(1+epsilon))* sigma_lt
    """
    f = A(M, b, c) + b * np.cos(2*x) + c * np.cos(x)
    return f


def fit_model(ax, func, x, y, xs, color, name):
    # Make model from function
    model = Model(func)
    # Make fit parameters
    params = model.make_params()
    # Make sure to set inital values to 1
    # so fit doesn't fail
    for p in params:
        params[p].set(value=1)

    # Fit the model
    try:
        out = model.fit(y, params, x=x)
    except ValueError as e:
        print(e)
        return None
    except TypeError as e:
        print(e)
        return None

    # Plot the fitted model with output parameters and same x's as model
    ax.plot(xs, out.eval(params=out.params, x=xs),
            linewidth=2.0, c=color, label=f'{name}', alpha=0.2)

    # Get uncertinty and plot between 2 sigmas
    # dely = out.eval_uncertainty(sigma=3, x=xs)
    # ax.fill_between(xs,
    #                 out.eval(params=out.params, x=xs)-dely,
    #                 out.eval(params=out.params, x=xs)+dely,
    #                 color=color, alpha=0.1,
    #                 label='$3 \sigma$')

    return out
