# A library to work with Gaussian from Python

This library was designed to work with Gaussian and its output files using python scripts and IPython notebooks

Two modules are available: GaussianCalc to start calculations and prepare input files and GaussWork to work with output (log) files

## Disclaimer

This library is for personal usage, so some of its methods may change dramatically. If you are using this code, be carefull with updating

## GaussianCalc

Note, that in order for this module to work correctly, path to Gaussian distribution should be stated in `GAUSS_EXEDIR` and `EXECUTABLE` variables. In future this may be unified.

Default amount of shared processors is stated by `SHARED_PROCS` variable, in future this may be modified to an amount of processors available in system or to some diistinct value (say, 1). However, methods are designed to provide npocs as an argument.

### Creating input files

Creating input file is done by method `write_input_params`. It must be provided with the coordinates file, containing atoms and their coordinates in any format accepted by Gaussian. Information on multiplicity and charge must also be written in this file

This method generates an input file defined by arguments. For now the only purpose of this method is to generate large amount of similar `.gjf` files with similar starting coordinates.

`write_multiple_input_params` can be used to create multiple input files with same method and job out of different coordinates file. In future methods for creating multiple input files with different jobs may be provided.

### Starting calculations

Method `start_calc` allows to start the calculation given a specific `.gjf` file. The python script will wait for Gaussian to finish it's calculations, so it may be used to make subsquent calculations. Environmental variable `GAUSS_SCRDIR` may be set by giving `scr_dir` argument. By default the method won't start calculation if `.out` file already exists, in this case the `CalculationError` will rise. This behaviour can be changed by setting `override` parameter to `True`.

Method `start_multiple_calc` can be used to subsuquently run several Gaussian jobs. It uses `start_calc` function, so parameters are similar. In future it may start several jobs in parallel (and may not).

## GaussWork

This module provides functionality to process Gaussian `.out` (`.log`) files.

### Retrieving optimization info

Using the `get_optimized_coords` the optimized coordinates of structure may be retrieved. By default this method will generate an error, if the calculations didn't converged (if there is no flag *Stationary point found*). This can be changed buy applying `ignore_error` parameter as `True`, in this case the last step of optimization will be used.

The coordinates will be returned as a `Molecule` object. You may probably need to modify the `ELEMENTS` dictionary.

### Writing coordniates info

`Molecule` object can be used to generate coordinates file for `write_input_params` file. Molecule, charge and multiplicity must be provided in order to do this. Default charge and multiplicity is 0, 1.

### Z-matrix

`Molecule` object can be used to generate Z-matrix. If some atoms need to be defined based on other than first three atoms, this can be applied by passing dictionary as a parameter. The dictionary must have a key equal to atom number (indexed by Gaussian) and value equal to a list of indexes of the atom's "basis" atoms. For example, if atom 10 must be defined based on distance to atom 5, angle 10—5—4 and torsion angle 10—5—4—6 (note that the order is important, check Z-matrix documentation to learn more), while all other atoms are defined based on atoms 1—3, the dictionary will be `{10: [5, 4, 6]}`.

### JSON

Atoms and molecules (and their coordinates) can be saved to JSON file and generated out of JSON file. It can be usefull for logging information.

### Retrieving free Gibbs' energy

`get_free_gibbs_energy` can be used to get information of compound's free gibbs energy of formation, if the `freq` job was performed. It returns data given by string *Sum of electronic and thermal Free Energies* as a float number.

### Scanning job

If the `scan` job was used, it's results can be obtained using `get_scan_results` function. It provides a dictionary of lists, defined by the Gaussian results table.
