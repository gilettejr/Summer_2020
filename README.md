# Summer_2020
Python code for working on UKIRT-WFCAM data

Must be run inside UKIRT-WFCAM directory, on my onedrive. Virtual Python environment for running also available.

Before anything else, run the reload_all_data function from reload_all_data.py. This activates the upload_from_scratch method in bookkeeping.py which processes the raw data, saving intermediate files as it goes.

Once this is complete, you can view cmds, ccs, spatial distributions etc by first creating an object containing the relevant dataset for the galaxy at a certain stage in the analysis process, then calling the relevant graphing function from agb_plotter.py

For example, if we want to look at the k-j cmd for NGC205 after the foreground colour cut has been made:

#create instance of class with chosen galaxy and stage of analysis desired


galaxy_object=data_read(galaxy='ngc205',stage='fore_cut')

#call the static method plot_kj_cmd in the basic_agb_plotter class to plot the data


basic_agb_plotter.plot_kj_cmd(galaxy_object)

The options for the galaxy argument are: cls_cut, fore_cut, agb, agb_crossed, cm
