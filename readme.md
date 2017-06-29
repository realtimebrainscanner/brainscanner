# Real-time EEG source localization
**Using a sparse Bayesian approach**

<!-- # Introduction -->
Real-time EEG source localization is a framework for performing real-time EEG source localization. Developed at Technical University of Denmark (DTU Compute) for recording EEG signals from a Smarting mBrainTrain EEG device. The developed program performs and visualizes source localization in real-time as well as provides real-time task classification.
Currently only Windows is supported.



## System
The used equipment is an EEG headset from Smarting mBrainTrain that connects via Bluetooth.
Device specifics:
- Sample rate: 250 Hz (500 Hz is also possible)
- 24 sensor EEG cap from EasyCap (https://www.mbraintrain.com/wp-content/uploads/2016/01/RBE-24-STD_legacy.pdf)


<!-- **System overview** -->
![](figures/systemoverview.png)

<!--- - [ ] insert figure of headset-->


The headset connects to a computer via Bluetooth using a BlueSoleil Bluetooth dongle that makes the device available as a Serial COM interface. See [mBrainTrain user manual](https://mbraintrain.com/wp-content/uploads/2016/08/SMARTING-User-Manual.pdf) for more details on this.

Using the OpenVIBE Acquisition Server, the data is sampled from the EEG device and made available through lab streaming layer ([LSL](https://github.com/sccn/labstreaminglayer)). By using the [lab streaming layer library for Matlab](https://github.com/sccn/labstreaminglayer/tree/master/LSL/liblsl-Matlab) the data is read into the Matlab application.


### Installation
Install the following programs
 - BlueSoleil Bluetooth driver (can be found on the accompanying USB drive, not the Bluetooth dongle)
 - OpenViBE [http://openvibe.inria.fr/downloads/](http://openvibe.inria.fr/downloads/). Replace the default OpenViBE Acquisition Server config file with the one supplied in this folder. Otherwise make sure that the driver settings are correct every time the OpenViBE Acquisition Server is started. To replace the file move *openvibe-acquisition-server.conf* to *C:\Users\[username]\AppData\Roaming\openvibe*
 <!-- - Lab streaming layer []() -->
 - Clone or download this repository
 - Install PsychoPy (http://www.psychopy.org/about/overview.html), NumPy, SciPy for creating and running stimulation files.


### Setting up the connection
More information can be found in the mBrainTrain user manual (see [user manual](https://mbraintrain.com/wp-content/uploads/2016/08/SMARTING-User-Manual.pdf)).

#### Setting up the Bluetooth connection:
Insert Bluetooth dongle and open "BlueSoleil Space". Connect to the smarting device by right clicking on the question mark (the dongle can only be matched with a specific EEG amplifier) and choose "Connect to Bluetooth seriel port (COMX)". Make a note of the seriel port number (e.g. "COM6" as in the below figure). Disconnect the serial port by first right clicking on the question mark and choose "Disconnect Bluetooth serial port (COM6)".

#### Using SmartingStreamer for impedance measurements 
Open the Smarting Streamer and click connect. Choose the correct seriel port and press the orange connect button. Start by clicking "Measure impedence on ref" and click "start streaming signals". Fill gell into the ground electrode and the reference electrode and see the color of the reference channel change from red to green. Click "stop streaming signals" and choose "Measure impedance" and "Select Channels" followed by "start streaming signals". When the electrodes have changed from red to green you are ready to view the EEG signal (click "show signals"). There will be an artefact from the impedance measurement if you have not changed to "No impendance measurement".
When you are satisfied with the impedances click disconnect in the top.

#### Setting up the OpenViBE connection
Open the "openvibe acquisition server". Click "Driver Properties" and change the port number to the seriel port noted earlier, press Enter before clicking Apply (otherwise the change will not be saved). This sted has to be conducted even though the openvibe-acquisition-server.conf file has been replaced. Click "Connect" and "Play", when this button becomes available.

#### Using the BrainScanner GUI
Open Matlab and go to the brainscanner folder. Type "BrainScanner" in the Command window. Click "Connect headset" and then "Start/stop". 

##### The Visualization window
Preprocessing steps of the visualized data can be activated by clicking (in the "Visualization" square) 
- "Filter": bandpass filters the data to [1,45] Hz)
- "Rereference": changes the montage to have average reference
- "ASR": performs artifact removal by artifact subspace reconstruction (Mullen 2015). ASR can only be activated if a data file has been loaded by clicking "Load ASR data". This button both imports data for ASR and trains the ASR model. Note that the ASR is trained with the preprocessing steps chosen in the "Visualization" square.

"Plot data" will show the EEG signal with the chosen preprocessing settings.
"Plot brain" will perform source localization using teVG (Hansen 2013, 2014) and then shows the localized cortical activitions. If teVG has been trained it will use the sparsity parameter obtained otherwise it will use a predefined value.

"Train VG" is best performed when the headset is not connected. Load an existing EEG file which you wish to use for training the VG parameter by clicking "Load file" in the "EEG data source" window. Choose the desired settings for training the model in the visualization window. Normally, at least the filter and rereference options should be chosen. Click "Train VG". The parameter estimation will take some time, depending on the length of the inputtet EEG data. A mat-file called "gamma" will be saved. Remember to rename of remove this file is you do not wish to use it as input for "Plot brain". 

The forward model used in both "Plot brain" and "Train VG" was constructed in BrainStorm using a the ICBM152 template and was build using OPENMEEG BEM headmodelling (Gramfort, 2010).

##### The Experiment window
Classification is performed after a classification model has been trained. The training will be conducted using the preprocessing options chosen in the experiment window. Normally at least "Filter" should be chosen (if the classification model uses source localization the "Rereference should also be chosen"). If ASR is needed the load a suitable data set for this and the click ASR. When the desired options have been chosen, then click "Train Model" ("Start/stop" should be off). First the user is asked if source localization should be performed. If "1" is typed VG will be train the sparsity parameter regardless if this has already been done in "Train VG".
A stimulation file and an EEG data file will be then be requested and finally the user will be asked whether the EEG data should be delayed with respect to the stimuli file. 

The classification model is build towards classifying open vs. closed eyes. Filtering is e.g. set to [5,15] Hz and the electrodes of interest are O1 and O2 (is source localization is not chosen).

##### Common warnings encountered in the BrainScanner app
The command window will display "input:noData" if the buffer does not contain data. Most of the times this just means that data was requested before new data was available. If you look in the log-file you often see that an empty data block is followed by a full data block (of e.g. 32 samples) and that a block of samples are in the buffer. The buffered data will be used/saved in the next iteration. 

To minimize disruption of data streaming have as few programs active as possible, and use the "Airplane mode" of the pc.

The message "input:shortBlock" will be displayed if less samples than the predetermined blocksize are available. This happens occasionaly and the system handles it by throwing away these samples. This favors realtime applications but will cause small shifts in the saved data with respect to timing and will thus affect offline analysis. This could be avoided by saving all data points.


### Running an open/closed eyes experiment

- Turn on "Airplane mode"
- Connect and disconnect to the correct EEG amplifier using BlueSoleil, open SmartStreaming and fill gel in electrodes while checking impedances. Disconnect SmartStreaming and start streaming through OpenVibe. See above for more details.
- Open BrainScanner in Matlab and connect headset and start streaming.
- Save 1 minute of data while the subject has open eyes and with as few artifacts as possible.
- Run and save an experiment training run.
- Set "Start/stop" to off.
- In the experiment window click "Filter", "Rereference", "Load ASR file" (choose the EEG data containing the 1 minute of data) and then click "ASR".
- Click "Train Model", type "1" in the command window for extracting features from source localization, choose the stimulation and EEG file used in the experiment training run. Delay by 5 seconds.
- When the training is completed turn on "Start/stop" this will now output "Closed eyes" or "Open eyes" in the Command window. A probality of closed eyes is also outputtet.
- To save a test run start saving EEG data at the same time as a stimulation file is run.


### Reading data

**"Data buffer"**

The data is read into a pseudo-buffer with a buffer size of 64 samples equal to two blocks (of each 32 samples) from the OpenViBE Acquisition Server.



## Real-time Matlab application
The Matlab [application](https://github.com/realtimebrainscanner/brainscanner/blob/master/BrainScanner.m)

![](figures/processing.png)





### Brain plotting

![](figures/brainSpin.gif)

## References
Mullen, T., Kothe, C., Chi, M., Ojeda, A., Kerth, T., Makeig, S., … Cauwenberghs, G. (2015). Real-time Neuroimaging and Cognitive Monitoring Using Wearable Dry EEG. IEEE Transactions on Biomedical Engineering, 62(11), 2553–2567. http://doi.org/10.1109/TBME.2015.2481482

Hansen, S. T., Stahlhut, C., & Hansen, L. K. (2013). Expansion of the Variational Garrote to a Multiple Measurement Vectors Model. In Twelfth Scandinavian Conference on Artificial Intelligence (pp. 105–114). Retrieved from http://scholar.google.com/scholar?hl=en&btnG=Search&q=intitle:Expansion+of+the+Variational+Garrote+to+a+Multiple+Measurement+Vectors+Model#0

Hansen, S. T., & Hansen, L. K. (2014). EEG source reconstruction using sparse basis function representations. In 2014 International Workshop on Pattern Recognition in Neuroimaging (pp. 1–4). IEEE. http://doi.org/10.1109/PRNI.2014.6858521

Gramfort, A., Papadopoulo, T., Olivi, E., & Clerc, M. (2010). OpenMEEG: opensource software for quasistatic bioelectromagnetics. Biomedical Engineering Online, 9(1), 45. http://doi.org/10.1186/1475-925X-9-45

## Acknowledgment

This work is supported by

[BaSiCs project by DTU and DRCMR](http://www.drcmr.dk/basics)

**Furthermore, we would like to give credit to**

- Carsten Stahlhult for the initial brain plotting functions and more
- labstreaminglayer
