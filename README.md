# AlexAnalysisPackage
Basic analysis package for SusyNt analyses.

## Dependencies
[SusyNt-Read](https://github.com/susynt/susynt-read)

[Restframes](https://github.com/crogan/RestFrames)

[Superflow](https://github.com/alexarmstrongvi/Superflow)

## Prerequisites:
 - svn access
 - kerberos ticket (run `klist` to check)
 - Run `kinit <USER>@CERN.CH` beforehand if no kerberos ticket established (you will be asked for your CERN password during setup_area otherwise. You may or may not need to capitalize your username to get things setup properly)
 - Access to cvmfs (for setting up an AnalysisRelease)
  

## Directions for initial setup 
1. Clone this repository
2. Run `setup_area.sh`
3. Configure `global_variables.py` in `python` directory
4. Run `setup_fax_samples.sh` - you will need to enter your GRID password

After initial setup, only run `setup_env.sh` and `setup_grid.sh` or `setup_fax.sh`

## Expected file structure
analysis 
- AlexAnalysisPackage
- Superflow - useful interface for making nTuples from susyNts, particularly for determining event weights
- susynt-read - central package to prepare a working area for reading SusyNtuples
- SusyNtuple - defines the objects stored in a root ntuple, as well as the methods used to retrieve the objects
- SUSYTools - tool class important for interfacing between xAOD and SusyNts
- RestFrames - 
- inputs - directory for storing text files of fax links to SusyNt samples
- inputs_local - directory for storing text files of local file path links to SusyNt samples

analysis\_run
- ntuples - root files used for plotting
- outputs - root files created by condor submissions
- logs - log files created by condor submissions
- plots - plot images created by plotting scripts
- yields - text files for yields

## Making nTuples for plotting
1. Run `setup_fax.sh`
2. Test on single sample file (`makeMiniNtuples -f inputs/group/sample.txt -n 1000`)
3. Go to top level directory and run `source make_tarball.sh`
4. Go into `analysis_run` directory and run `source clear_for_submission.sh`
    * remove all files from `logs` and `output` instead if you are not interested in keeping them
5. Go into `outputs` directory and run `python submitMakeMiniNtupleToCondor.py`
    * output root files will be stored in `outputs` directory
    * .log, .err, and .out files for each submitted sample are stored in the `logs` directory
6. Run `condor_q $USER` to check on jobs
7. After jobs complete, identify failed samples with `python check_for_failed_samples.py` and print out causes of failures with `print_failed_sample_stats.py`

## Helpful Scripts
### Plotting
### Yield Table
### Checking DSID lists
### Cutflow Comparison
### Other

## Writing configuration files for plotter.py
