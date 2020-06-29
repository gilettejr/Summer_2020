# Summer_2020
Python code for working on UKIRT-WFCAM data

data_load is base class, data_read and bookkeeping inherit off this. Contains functions for reading in and processing UKIRT-WFCAM files (filepath will need to be altered accordingly for code to run). Also allows intermediate data to be saved, which can then be read in using data_read. This is much faster than using data_load, as cuts are not repeatedly made each time the code is run.
