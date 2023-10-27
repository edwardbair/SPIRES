
# modSPIReS
modifications in progress to SPIReS codes

Added option in the fixpeak section, as follows:
- find peaks in the rgvec and fSCA vectors for a single pixel
- If grain size peaks after fSCA peaks, do not let subsequent grain sizes fall below that peak.
- ElseIf fSCA peaks after grain size peaks, keep the grain size at that fSCA peak as the minimum until the end.

To leave Ned's code as is, added NetMethod=false at line 202. Change this to NedMethod=true to run Ned's version.


20231027 Ned uploads a beta copy with changes
