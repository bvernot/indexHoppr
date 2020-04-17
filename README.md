# indexHoppr


### To install:

Go to the directory where you want this to live, and run:

~~~~
git clone https://github.com/bvernot/indexHoppr
cd indexHoppr
~~~~

Start `R` in that directory, and run:

~~~~
install.packages("devtools")
devtools::install('.')
~~~~

Now try loading the library:

~~~~
library(indexHoppr)
~~~~

If it works, exit out of R. If it doesn't.. 


### To run:

From any directory, you should be able to run this and get the help message:

~~~~
Rscript /path/to/github/repos/indexHoppr/findContamCLI.R -h
~~~~

An example command:

~~~~
time Rscript /path/to/github/repos/indexHoppr/findContamCLI.R \
  --splits /mnt/sediments/analyzed_runs/190221_M06210_B24425_MTcapHuman/split/splittingstats.txt \
  --num-plot-libs 30  \
  --plots MTcapHuman_ContamPLOTS_190221_M06210_B24425_A19070_A19073.pdf \
  -libs A19070 A19073 \
  --table MTcapHuman_ContamSCORES_190221_M06210_B24425_A19070_A19073.tsv \
  -nc 40 
  -rcf 20
~~~~
