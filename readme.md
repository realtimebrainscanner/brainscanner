# Real-time EEG source localization
**Using a sparse Bayesian approach**

<!-- # Introduction -->
Real-time EEG source localization is a framework for performing real-time EEG source localization. Developed at Technical University of Denmark (DTU Compute) for recording EEG signals from a Smarting mBrainTrain EEG device, performing, and visualizing source localization in real-time.



## System
The used equipment is an EEG headset from Smarting mBrainTrain that connects via Bluetooth. 
Device specifics:
- Sample rate: 250 Hz (500 Hz also possible)
- 24 sensor EEG cap from EasyCap (22 channels + 2 mastoid)


<!-- **System overview** -->
![](figures/systemoverview.png)

- [] insert figure of headset


The headset connects to a computer via Bluetooth using BlueSoleil Bluetooth dongle that makes the device available as a Serial COM interface. See [mBrainTrain user manual](https://mbraintrain.com/wp-content/uploads/2016/08/SMARTING-User-Manual.pdf) for more details on this.

Using the OpenVIBE Acquisition Server, the data is sampled from the EGE device and made available through lab streaming layer ([LSL](https://github.com/sccn/labstreaminglayer)). By using the [lab streaming layer library for Matlab](https://github.com/sccn/labstreaminglayer/tree/master/LSL/liblsl-Matlab) the data is read into the Matlab application. 


### Setting up the connection
Follow the guide by mBrainTrain to set up the BlueSoleil dongle and connect to the headset. 
Replace the OpenViBE default config file with [this one](). 

- [] Insert screen shot of OpenVIBE setup!!!


### Reading data

**"Data buffer"**
The data is read into a pseudo-buffer with a buffer size of 64 samples equal to two blocks (of each 32 samples) from the OpenViBE Acquisition Server. 



## Real-time Matlab application
The Matlab [application]() 

![](figures/processing.png)





### Brain plotting

![](figures/brainSpin.gif)


## Credits

We would like to give credit to
- Carsten Stahlhult for the initial brain plotting functions and more
- labstreaminglayer
