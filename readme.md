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


### Setting up the connection
#### Setting up the Bluetooth connection:
Insert Bluetooth dongle and open "BlueSoleil Space". Connect to the smarting device by right clicking on the question mark (the dongle can only be matched with a specific EEG amplifier) and choose "Connect to Bluetooth seriel port (COMX)". Make a note of the seriel port number (e.g. "COM6" as in the below figure). Disconnect the serial port by first right clicking on the question mark and choose "Disconnect Bluetooth serial port (COM6)".


#### Using SmartingStreamer for impedance measurements 
Open the Smarting Streamer and click connect. Choose the correct seriel port and press the orange connect button. Start by clicking "Measure impedence on ref" and click "start streaming signals". Fill gell into the ground electrode and the reference electrode and see the color of the reference channel change from red to green. Click "stop streaming signals" and choose "Measure impedance" and "Select Channels" followed by "start streaming signals". When the electrodes have changed from red to green you are ready to view the EEG signal (click "show signals"). There will be an artefact from the impedance measurement if you have not changed to "No impendance measurement".
When you are satisfied with the impedances click disconnect in the top.

#### Setting up the OpenViBE connection
Open the "openvibe acquisition server". Click "Driver Properties" and change the port number to the seriel port noted earlier, press Enter before clicking Apply (otherwise the change will not be saved). This sted has to be conducted even though the openvibe-acquisition-server.conf file has been replaced. Click "Connect" and "Play", when this button becomes available.

#### Using the BrainScanner GUI
Open Matlab and go to the brainscanner folder.



More information can be found in the mBrainTrain user manual (see [user manual](https://mbraintrain.com/wp-content/uploads/2016/08/SMARTING-User-Manual.pdf) page 9).



Select the correct COM-port in the driver settings of OpenViBE Acquisition Server. Always press *Enter* and click the *Ok* button, else the change is not registered by the program.

- [] Insert screenshot of OpenVIBE setup!!!


For more about the possible settings, see []().


### Reading data

**"Data buffer"**

The data is read into a pseudo-buffer with a buffer size of 64 samples equal to two blocks (of each 32 samples) from the OpenViBE Acquisition Server.



## Real-time Matlab application
The Matlab [application](https://github.com/realtimebrainscanner/brainscanner/blob/master/BrainScanner.m)

![](figures/processing.png)





### Brain plotting

![](figures/brainSpin.gif)


## Acknowledgment

This work is supported by

[BaSiCs project by DTU and DRCMR](http://www.drcmr.dk/basics)

**Furthermore, we would like to give credit to**

- Carsten Stahlhult for the initial brain plotting functions and more
- labstreaminglayer
