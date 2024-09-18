# 基于声码器Vocoder原理进行的音频DSP项目
利用STFT（短时距傅里叶变换）以及多种插值方法，在保证音质的情况下，分部实现了TimeStretcher用于改变音频的速度，PhaseVocoder相位声码器，以及Voice_Changer变声器。并以变声器实现逻辑为基础，使用JUCE构建变声器插件.

从TS改进到相位声码器的过程本质上是一个将相邻帧直接相位缝起来的过程
1. 分帧：语音信号短时平稳性 10-30ms 相邻帧重叠的比例通常是帧长的一半或者75%，以保证连续性，同时减少因分帧带来的边界效应。overlap的目的是避免信号中具有一些跳变信号导致信号分析时的准确度不够,说白点就是提高分辨率；       加窗（矩形窗，海宁窗，汉明窗）：窗函数也叫做加权函数，可以让信号更好地满足FFT处理的周期性要求；截断都有频域能量泄露，而窗函数可以减少截断带来的影响
2. STFT（基2傅里叶变化）:获取波的相位信息，以用于估计瞬时频率
3. 重采样方法进行变速变调。重采样是一种变调变速的音频处理办法,而我们变声器需要的是变调不变速的办法,当然,变声器未必一定需要使用重采样来完成,还有诸如频域变换法,正弦建模法等办法来完成,但前者容易在频域引入额外的分量,后者的计算量也较为庞大
4. 时域压扩：对每个帧进行一系列处理比如拉伸或者压缩，最后在将这些帧重新叠加成合成信号。这个重叠部分需要进行一系列的处理以减少类似于相位不连续、幅度波动造成的影响。 （e.g.比如要把一段音频信号的语速放慢，只需要每一个帧都进行延扩，然后在首位拼接起来就行了，但造成的劣势也是显而易见的，它会造成拼接处的波形不连续（语音信号中的基音断裂））； 这里重叠相加算法是基于相位声码器的原理，因为WSOLA的缺点，在作用于冲击瞬态信号如打击乐等时合成失真明显(原音频信号截取一帧后，通过波形相似匹配下一帧，但两个帧都包含一个瞬态的音频信号)
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
