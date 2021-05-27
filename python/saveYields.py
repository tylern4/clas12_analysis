#!/usr/bin/env python
import warnings  # noqa
warnings.filterwarnings("ignore")  # noqa

from calcXsections import *
from tqdm import tqdm
import argparse
import os


def main(rec: pd.DataFrame,
         binning,
         out_folder: str = "plots",
         bins: int = 12,
         csvName: str = "results"):

    results: list = []

    if not os.path.exists(f'{out_folder}/crossSections'):
        os.makedirs(f'{out_folder}/crossSections')

    for w in tqdm(binning["wbins"]):
        # Cut data w bin we're in
        rec_w = rec[rec.w_bin == w].copy()

        for q2 in binning["q2bins"]:
            # Cut data for the q2 bin we're in
            rec_wq2 = rec_w[rec_w.q2_bin == q2].copy()

            for theta in binning["thetabins"]:
                # Cut data/mc for the theta bin we're in
                data = rec_wq2[rec_wq2.theta_bin == theta].copy()

                _data_y, _x = hist_data(data, density=False, bins=bins)

                # Get numbers for error
                N_y = _data_y

                # Remove points with 0 data count
                cut = ~(_data_y == 0)
                x = _x[cut]
                data_y = _data_y[cut]
                N_y = N_y[cut]

                # Get bin widths
                delta_W = (w.right-w.left)
                delta_Q2 = (q2.right-q2.left)
                delta_Theta = np.abs(theta.right-theta.left)
                __phis = np.linspace(0, 2 * np.pi, bins)
                delta_phi = __phis[1] - __phis[0]
                kin_bin_width = delta_W * delta_Q2 * delta_Theta * delta_phi

                # Calculate acceptance and correct data
                flux = virtual_photon_flux(w.mid, q2.mid)
                denom = kin_bin_width * flux

                # Normalize with bin widths
                try:
                    y = data_y / denom
                except ValueError:
                    continue

                error_bar = np.sqrt(N_y)

                for phi, cross, err in zip(x, y, error_bar):
                    results.append({"w_left": w.left,
                                    "w_right": w.right,
                                    "q2_left": q2.left,
                                    "q2_right": q2.right,
                                    "cos_theta": theta.left,
                                    "x": phi,
                                    "y": cross,
                                    "err": err})

    output = pd.DataFrame(results)
    output.to_csv(f"{out_folder}/crossSections/{csvName}.csv")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Make Cross Sections")
    parser.add_argument("--mc", dest="mc_data_file_path",
                        type=str, help="MC csv file", required=False)
    parser.add_argument("--data", dest="rec_data_file_path",
                        type=str, help="Data csv file", required=True)
    parser.add_argument("--empty", dest="empty_file_path",
                        type=str, help="Empty run csv file", required=False)
    parser.add_argument("--radcorr", dest="radcorr", type=str,
                        help="Location of radcorr data csv",
                        required=False, default=None)
    parser.add_argument("--out_folder", dest="out_folder", type=str,
                        help="Location final csv", required=False, default=".")
    parser.add_argument("--phibins", help="Number of phi bins",
                        dest="phibins", required=False, default=10, type=int)
    args = parser.parse_args()

    bins = args.phibins

    w_bins = w_bins_e99
    q2_bins = q2_bins_e99
    csvName = f"full_results_{bins}"
    print(csvName)
    # Start to main

    # print("Start setup")
    start = time.time_ns()
    # Load reconstructed file
    rec = read_csv(args.rec_data_file_path, True)

    # Cut for missing mass
    rec = cut_for_MM(rec)
    end = time.time_ns()

    # Make bins in the dataframes from the bins above
    _rec = prep_for_ana(rec, w_bins, q2_bins, theta_bins)

    # Create dict of sorted bins
    _binning = dict()
    _binning["wbins"] = pd.Index.sort_values(pd.unique(_rec.w_bin))
    _binning["q2bins"] = pd.Index.sort_values(pd.unique(_rec.q2_bin))
    _binning["thetabins"] = pd.Index.sort_values(pd.unique(_rec.theta_bin))
    print(f"Done setup: {(end-start)/1E9:0.2f}Sec")

    main(_rec, _binning,
         bins=bins,
         out_folder=args.out_folder,
         csvName=csvName)

    del _rec
