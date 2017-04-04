# Real-time EEG source localization
**Using a sparse Bayesian approach**

<!-- # Introduction -->
Real-time EEG source localization is a framework for performing real-time EEG source localization. Developed at Technical University of Denmark (DTU Compute) for recording EEG signals from a Smarting mBrainTrain EEG device, performing, and visualizing source localization in real-time. 
Currently only Windows is supported. 



## System
The used equipment is an EEG headset from Smarting mBrainTrain that connects via Bluetooth. 
Device specifics:
- Sample rate: 250 Hz (500 Hz also possible)
- 24 sensor EEG cap from EasyCap (22 channels + 2 mastoid)


<!-- **System overview** -->
![](figures/systemoverview.png)

- [] insert figure of headset


The headset connects to a computer via Bluetooth using BlueSoleil Bluetooth dongle that makes the device available as a Serial COM interface. See [mBrainTrain user manual](https://mbraintrain.com/wp-content/uploads/2016/08/SMARTING-User-Manual.pdf) for more details on this.

Using the OpenVIBE Acquisition Server, the data is sampled from the EEG device and made available through lab streaming layer ([LSL](https://github.com/sccn/labstreaminglayer)). By using the [lab streaming layer library for Matlab](https://github.com/sccn/labstreaminglayer/tree/master/LSL/liblsl-Matlab) the data is read into the Matlab application. 


### Installation 
 - BlueSoleilt [..]()
 - OpenViBE [http://openvibe.inria.fr/downloads/](http://openvibe.inria.fr/downloads/)
 - Lab streaming layer []()
 - This repository
 - Move *openvibeConfig* to *~/...*


### Setting up the connection
Follow the guide by mBrainTrain (see [user manual](https://mbraintrain.com/wp-content/uploads/2016/08/SMARTING-User-Manual.pdf) page 9) to set up the BlueSoleil dongle and connect to the headset. 

Replace the OpenViBE default config file with [this one](). 

Select the correct COM-port in the driver settings of OpenViBE Acquisition Server. Always press *Enter* and click the *Ok* button, else the change is not registered by the program. 

- [] Insert screen shot of OpenVIBE setup!!!

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
