# MERFISH_simulator

A simulator for generating MERFISH data.

**Warning:** Alpha version, may change rapidly and unpredictably.

## Dependencies

* [MERlin](https://github.com/ZhuangLab/MERlin)

## Usage

This requires an existing MERlin `merfish-parameters` directory.

1. Create a sub-directory in `merfish-parameters` called `simulation`
and put the your simulation parameters json file in this directory.

2. Place a csv file containing positions in the `merfish-parameters/positions` directory.

3. Run simulator as follows:
```
python ~/path/to/mersim/simulator.py \
   -a simulation_parameters.json \
   -m microscope_parameters.json \
   -o data_organization.csv \
   -c codebook.csv \
   -p positions.csv \
   simulate
```

This will create a directory called `simulate` in your MERlin `data` directory containing the simulated images.
