# AlexAnalysisPackage
Basic analysis package for SusyNt analyses.

## Dependencies
https://github.com/susynt/susynt-read
https://github.com/dantrim/Superflow
https://github.com/crogan/RestFrames

## Expected file structure
analysis\_<tag> 
- AlexAnalysisPackage
- Superflow
- susynt\_read
- SusyNtuple
- SUSYTools
- RestFrames
- inputs

analysis\_<tag>\_run
- flat_ntuples - root files used for plotting
- outputs - root files created by condor submissions
- plots - plot images created by plotting scripts
-  

## Directions for setup 
Once susynt-read is checked out, and compiled, check out Superflow and this package, then simply recompile.

Once successful, simply do:
    makeMiniNtuple -f file_name -n number_of_events
