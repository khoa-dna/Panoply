#!/bin/bash
echo "Start Running Panoply on Multiple Configs"
config_dir=$1
num_panel=$2
output_folder=$3
tag_name=$4
SCRIPT_DIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )

echo $SCRIPT_DIR
for config in $config_dir/*
do 
echo run config: $config
config_name=$(basename "${config%.}")
nohup python $SCRIPT_DIR/run_panoply_cmd.py $config $output_folder "$tag_name"_${config_name%.*} --num-panel $num_panel &
done