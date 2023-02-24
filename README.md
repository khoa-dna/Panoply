# Panoply

## Overview

Panoply generates complete and optimized panels for flow cytometry experiment given a target list of markers by heuristcally searching a provided database of marker-fluor pairs

Constrained beam search algorithm with heuristic function is used to search for target panels. The algorithm initializes with Auto fluorescence marker/fluor pair and continuously adding one new fluor-marker (f-m) pair to "grow" the panel until a complete panel (if any) is generated. At each growing step, marker with the least fluor partners are chosen first. If there are multiple possible fluors that can be chosen for a marker, heuristic function is used to calculate a score which determine best fluors to be added. 

Heuristic funtion is a function of complexity index (CI), fluors-marker intensity, and fluors brightness-marker abundance correlation, etc... (to be optimized)

Some randomness are added in the process to avoid convergence to one single panel and local maximum. 

## Data Files

- Database: Contains all available marker-fluor pairs
- Spectra: Contains spectrum of (ideally) all fluors in Database
- Marker_Abundance: Contains abundance info of (ideally) all markers in Database
- Fluor_Brightness: Contains brightness info of (ideally) all fluors in Database
- Target: target list of markers

## Understand params

There are 2 sets of params. One is used for the search algorithm and the other is used for synthetic data simulation.

### params for search algorithm

- Intensity weight: part of heuristic function. This is the coef of how much f-m intensity influence the selection of next f-m pair. f-m intensity is defined as the number of standard deviation away from the mean of all f-m instensity. 

- Correlation weight: part of heuristic function. This is the coef of how much correlation between fluors brightness-marker abundance affects the choosing of the next f-m pair.

- Random prob: probability of not choosing the best f-m pair according to heuristic function. For example, if there are 3 possible fluors to match with a marker at step 1, heuristic function will rank the 3 pair. If random prob = 0.0, the best pair will always be chosen at each step. This means that the algorithm will always converge to 1 single panel every time it runs. If random_prob = 0.1, there is 10% chance that non-optimal pair is chosen at each growing step.

- Risk prob: (only applicable when branch num > 1). When choosing the panels to be put into Frontier, if risk_prob = 0.0, the best n panels (n = Frontier size) are chosen. If risk_prob = 0.1, (1 - 0.9^n) chance that at least 1 non optimal panels are chosen.

- Branch num: for each panel in Frontier, how many best panels are generated

- Frontier size: out of all generated panels, down select to number of panels equal to Frontier size

Consider an example to understand the affect of each param on the algorithm. Let say we have the following param: branch_num = 2, frontier_size = 3, random_prob = 0.1, risk_prob = 0.1. Initially, one panel starts with Autofluorescence f-m pair. At the next growing step, say there are 5 possible fluors to match with the next marker (marker with least partners). The score for these 5 size-2 panels (autofluor and current marker) are calculated by heuristic function. Because branch_num is set to 2, 2 out of 5 fluors will be chosen without replacement. However, since random_prob = 0.1, at each time a fluor is chosen, there are 10% chance that it is not the one with highest score. Next, since frontier_size = 3, all 2 panels will be chosen into the Frontier. Each panel will contnue to branch out 2 panels. Now, we have 4 panels, but frontier size is 3, so we need to only choose 3 best panels to continue. However, since risk_prob is 0.1. At each time a panel is chosen, there are 10% chance that it will not be the panel with the highest score.

### params for synthetic data

TBA

## Scripts

### Run Panoply with GUI

```
python run_panoply.py
```

This will open up GUI. You can either load an existing config json to populate the data files and params as this [example config](https://github.com/khoa-yellow/Panoply/blob/master/data/panel_30_0818.json) or input your own settings.


### Run Panoply using command line

```
usage: run_panoply_cmd.py [-h] [--num-panel NUM_PANEL] config_json output_folder tag_name

positional arguments:
  config_json           path to config json with params contains list of values want to simulate
  output_folder         output folder for generated panels and diagnostics plot
  tag_name              name unique to this run
  optional arguments:
  -h, --help            show this help message and exit
  --num-panel NUM_PANEL
                        num panels to generate
```

Example

```
python .\src\run_panoply_cmd.py .\data\panel_30_0818.json .\temp\ example_run --num-panel=10 
```

This will generate in the `temp` folder the following files
- `example_run_cis.csv` : contains CI of each panel
- `example_run_panles.csv` : contains Fluors/Markers of each panel
- `example_run_qc.pdf`: histogram of CIs (to be added more diagnostics plots)

### Run Panoply with multiple config files

You may find yourself wanting to run multiple params settings. Follow 2 steps to achieve this.

#### Step 1: Generate multiple configs

```
usage: parsim.py [-h] config_json output_folder

positional arguments:
  config_json    path to config json with params contains list of values want to simulate
  output_folder  output folder contains all simulated config json

optional arguments:
  -h, --help     show this help message and exit
```

This config json is different from the one used previously when running command line or GUI mode. In this config json, some fields may be specified as a list. All combinations of lists among these fields equals the number of simulated config files. [This](https://github.com/khoa-yellow/Panoply/blob/master/data/panel_30_0818_parsim.json) is an example config file for parsim.py and [this](https://github.com/khoa-yellow/Panoply/tree/master/data/simulated_configs) is an example output of parsim.py

Example
```
python .\src\parsim.py .\data\panel_30_0818_parsim.json .\data\simulated_configs\
```

#### Step 2: Run run_panoply_cmd.py on generated configs in step 1 in parallel

To achieve this, run 
```
source run_panoply_multi_config.sh <config_dir> <num_panel> <output_folder> <tag_name>
```

- `config_dir`: directory output in step 1
- `num_panel`: number of panels want to generate
- `output_folder`: folder contain all outputs from running panoply on each config file
- `tag_name`: unique name for this run

Example
```
source run_panoply_multi_config.sh ../data/simulated_configs 10 ../temp test_multi
```