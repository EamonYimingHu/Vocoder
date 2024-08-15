基于声码器Vocoder进行的音频DSP处理。利用STFT（短时距傅里叶变换）以及插值方法，分部实现了TimeStretcher用于改变音频的速度，PhaseVocoder相位声码器，以及变声器。并以变声器实现逻辑为基础，使用JUCE构建变声器插件

# Time_Stretcher

Start with a basic Time Stretcher. The problem is that we are not accounting for changes in the phase of our 
frames as we shift them! So it's been improved by us to make the vocoder.

# Phase_Vocoder

A basic Phase Vocoder is implemented. We are now changing the audio's speed without changing it's pitch. 
But there is one more improvement to be made, and that is we can make the sound more realistic by interpolation and resampling.

# Voice_Changer

A high quality voice changer is implemented here. The basic Phase Vocoder is modified so that it performs uniform pitch shifting of the audio data. 
This is achieved by resampling each of the output frames to a length N*RA/RS using interpolation methods. 
Through a series of signal processing methods, finally implemented changing pitch without changing speed.

# Vocoder_with_time-dependent_variable_parameters

Based on the existing vocoder, some interesting sound effects can be achieved by setting Q to time-dependent variable parameters.
