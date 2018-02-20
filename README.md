# AlexAnalysisPackage
Basic analysis package for SusyNt analyses.

## Dependencies
[SusyNt Read](https://github.com/susynt/susynt-read)

[Restframes](https://github.com/crogan/RestFrames)

[Superflow](https://github.com/dantrim/Superflow)

## Expected file structure
analysis 
- AlexAnalysisPackage
- Superflow - useful interface for making nTuples from susyNts, particularly for determining event weights
- susynt-read - central package to prepare a working area for reading SusyNtuples
  - SusyNtuple - defines the objects stored in a root ntuple, as well as the methods used to retrieve the objects
  - SUSYTools - tool class important for interfacing between xAOD and SusyNts
- RestFrames - 
- inputs - directory for storing text files of fax links or local file paths to SusyNt samples

analysis\_run
- ntuples - root files used for plotting
  - data
  - mc
  - logs
- outputs - root files created by condor submissions
- logs - log files created by condor submissions
- plots - plot images created by plotting scripts
  

## Directions for initial setup 
1. Clone this repository
2. Define DSID groups in global_variables.py
3. Setup environment
    1. Kinit
    2. (Optional) rucio and fax
4. Run setup_area.sh
5. (?) Configure global_variables.py assuming it can’t all be done automatically
6. (optional) Run setup_fax_samples.sh if the plan is to use SusyNts stored on fax instead of locally

## Making nTuples for plotting
1. Setup fax and rucio
2. Test on single fax link (`makeMiniNtuples -f fax-prefix/path/to/group:group.phys-susy.ID#._00000#.susyNt.root -n 1`)
3. Go to top level directory and run `source make_tarball.sh`
4. Go into `analysis_run` directory and run `source clear_for_submission.sh`
    * remove all files from `logs` and `output` instead if you are not interested in keeping them
5. Go into `outputs` directory and run `python submitMakeMiniNtupleToCondor.py`
    * output root files will be stored in `outputs` directory
    * .log, .err, and .out files for each submitted sample are stored in the `logs` directory
6. Run `condor_q $USER` to check on jobs

- Optional checks
  - `chain_process.py` - looks through `outputs` and creats TChains for each background and store missing root files
  - `check_for_failed_samples.py` - read log files for confirmation of success and record failed samples
  - `print_failed_sample_stats.py` - read through missing or failed log files to determine causes of failure

## Writing configuration files for plotter.py
