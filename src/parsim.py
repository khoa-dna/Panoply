import json
from itertools import product
import argparse
import os
from os.path import join, basename, dirname


def arg_parser():
    parser = argparse.ArgumentParser()
    parser.add_argument('config_json', type=str,
                        help="path to config json with params contains list of values want to simulate")
    parser.add_argument('output_folder', type=str,
                        help="output folder contains all simulated config json")
    args = parser.parse_args()
    return args


def save_json(data, path, indent=2):
    with open(path, "w") as f:
        json.dump(data, f, indent=indent)


def main():
    args = arg_parser()
    os.makedirs(args.output_folder, exist_ok=True)
    # load config file
    with open(args.config_json, "r") as f:
        config_parsim = json.load(f)
    # obtain params with list values
    sim_params = {key: val for key, val in config_parsim.items()
                  if isinstance(val, list)}
    header = list(sim_params.keys())
    input_config_name = basename(args.config_json).split(".")[0]
    # iterate all combinations of sim params
    for sim_par in product(*list(sim_params.values())):
        config_name = ""
        # overwrite with single param value
        for i, val in enumerate(sim_par):
            config_parsim[header[i]] = val
            config_name += f"_{header[i]}={val}"
        # save output_file
        output_file = join(args.output_folder,
                           f"{input_config_name}{config_name}.json")
        save_json(config_parsim, output_file)


if __name__ == "__main__":
    main()
