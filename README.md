# Time_Stretcher

Start with a basic Time Stretcher. The second problem is that we are not accounting for changes in the phase of our 
frames as we shift them! 

# Phase_Vocoder

A basic Phase Vocoder is implemented. Implemented a vocoder effect that changes key without changing speed.

# Voice_Changer

A high quality voice changer is implemented here. The basic Phase Vocoder is modified so that it performs uniform pitch shifting of the audio data. 
This is achieved by resampling each of the output frames to a length N*RA/RS. Through a series of signal processing methods, finally implemented changing speed without changing pitch.

# Vocoder_with_time-dependent_variable_parameters

Based on the existing vocoder, some interesting sound effects can be achieved by setting Q to time-dependent variable parameters.
