# bash calls this exit status.
 test if the errorlevel q returns (0 if all ok) are fine
 eg 0 all ok:
 - 1. if shake error
 - 2. if NaN error
 - 4. if interrupted
 - 8. if crashed
 - 16. if hotatoms during fep, etc
 
also important: could not write/read/out of space/storage hangs/ etc pp
