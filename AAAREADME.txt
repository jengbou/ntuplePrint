go into root
root.exe
.L EMJselect.C
EMJselect("ntuple.root","histos.root",1000,0.2,0.9,0.9,1)
exit root and go into it again to look at the histograms it made
root.exe histos.root
TBrowser a



to see an example of how to properly combine the hists from QCD samples made in the 5 HT bins, look at the QCDcombiner subdirectory of this area and look at the AAAREADME there.
