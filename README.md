# AlexAnalysisPackage
Basic analysis package for SusyNt analyses.

## Dependencies
[SusyNt Read](https://github.com/susynt/susynt-read)

[Restframes](https://github.com/crogan/RestFrames)

[Superflow](https://github.com/dantrim/Superflow)

## Expected file structure
analysis 
- AlexAnalysisPackage
- Superflow
- susynt-read
- SusyNtuple
- SUSYTools
- RestFrames
- inputs

analysis\_run
- flat\_ntuples - root files used for plotting
- outputs - root files created by condor submissions
- plots - plot images created by plotting scripts
  

## Directions for setup 
1. Open up clean environment 
2. Create directories `analysis_<susy_tag>` and `analysis_<tag>_run`
3. cd into `analysis_<tag>`
4. Follow Directions on [SusyNt Read](https://github.com/susynt/susynt-read) 
    * Initializer kerberos (e.g. `kinit ALARMSTR@CERN.CH`)
    * Change [--dev] to [—stable] when running setup_area.sh
5. Mv RootCore, SusyNtuple, and SUSYTools out of susynt-read into parent directory
6. Git clone Superflow and AlexAnalysisPackage
7. Git clone RestFrames and follow setup instructions in INSTALL file
8. make `inputs` directory and run `python susynt-read/python/available_datasets` within that directory
    * Make sure rucio is setup
9.  Make directories for all the background groupings used in the analysis and a general MC directory
10. Run `make_condor_lists.py` over txt file in `inputs_LFV` and output them into their corresponding directory
    * Make sure that fax is setup
11. Run `condor_sample_groupings.py` on the MC output directory from (10) and input into the `inputs` directory
    * Edit defined groups to match the desired DSID lists
    * Some of the subprocess commands may be commented out so uncomment those
12. Edit `submitMakeMiniNtupleToCondor.py` to grab the right files in the right directory
13. Go into `analysis_<tag>_run` and create the directories: logs, outputs, plots, and samples
14. Set the relevant global variable values in `scripts/global_variables.py`
15. Run `make_dsid_filelists_from_inputs.py`
17. Add any manually input cross sections from `SUSYtools/data/mc15_13TeV` into the corresponding tag’s file
18. Setup RootCore (`source RootCore/script/setup.sh`)
19. Setup Root (`source susynt-read/bash/setup_root.sh` from within `analysis_<tag>`)
20. `rc clean` and `rc find_packages` and then `rc compile` the analysis directory
    * Remove superflow analyses in the Superflow if any 
21. Test that `makeMiniNtuples_LFV` works by running on a specific fax link with 1 event
22. Go outside the `analysis_<tag>` directories and `tar analysis_<tag>`
    * `tar cfvz analysis_n<tag>.tgz  --exclude=".*" --exclude="*git*" --exclude="*\.root" analysis_n<tag>/`
    * Make sure there is a soft link area.tgz pointing to the tarred analysis directory

Getting miniNtuples and plots
1. If any changes
    * `source susynt-read/bash/setup_root.sh` from `analysis_<tag>`
    * Compile `analysis_<tag>`
2. Remove all files from outputs and log directory
3. run condor submission script from outputs directory
    * Check status with `condor_q $USER`
5. run `chain_process.py`
6. hadd all root files into signal file
7. run plotting scripts
