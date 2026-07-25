# Acoustic 2D Drone Localization Tracking System
## Project Overview
This repository contains an array signal processing pipeline designed to calculate and track the 2D spatial trajectory of a drone. Developed as a core academic engineering project within the joint SICOM (Signal, Image, Communication, Multimedia) branch at Grenoble INP (Phelma / ENSE3). 
Poor recording conditions have heavily impacted the possible paths to detecting the trajectory. I tried TDOA, intensity differences, and DOA from every mic combinations without finding a viable trajectory for the drone.


## Repository Structure
Subject 25-26.pdf               PDF presenting the project
Audio_Track-1..4.wav            Raw recordings from the 4-microphone array
filtrage.m                      Zero-phase bandpass filtering of the raw tracks → Audiofiltre1..4.wav
initialisation.mlx              Loads filtered audio, defines mic geometry, computes all pairwise τ via windowed cross-correlation
DOA 12-34 mic pair.m            Takes pre-calculated time delays, computes the DOA angles for the 12-34 mic pairs per time window, and graphs the resulting 2D trajectory.
Drone trajectory video.mp4      Video of the proposed drone trajectory resulting from the Direction of Arrival (DOA) method
intensity_attempt.m             Attempt at localizing the drone from relative signal intensity/energy across the 4 microphones — not conclusive 
DOA_any_mic_pair_attempt.m      Generalization of the DOA method to all microphone pair combinations (not just 12/34) — not conclusive 
TDOA_attempt.m                  Attempt at TDOA multilateration: pairwise delays converted to hyperbolic constraints, solved via nonlinear least squares — not conclusive


## Usage
Place the 4 raw .wav recordings in the working directory.
Run filtrage.m to produce Audiofiltre1.wav … Audiofiltre4.wav.
Run initialisation.mlx to compute the τ delays for every microphone pair.
Run DOA_12_34_mic_pair.m to estimate DOA angles and visualize the resulting trajectory. 
TDOA_attempt.m, intensity_attempt.m, and DOA_any_mic_pair_attempt.m are the exploratory scripts referenced in Limitations — after initialisation.mlx.


## Limitations
*  Wind noise dominates the 0–150 Hz band and overwhelms the drone's acoustic signature across most of the recording.
*  Some channels show clipping, which distorts the cross-correlation peaks used for both TDOA and DOA.
*  The 0.2 m microphone spacing gives a very narrow physically-valid delay range, making delay estimates highly sensitive to noise.
*  None of the three methods produced a trajectory that seemed plausible.

## Technologies Used
*   **Languages:** MATLAB
*   **Key Libraries:** MATLAB Signal Processing (xcorr, butter, filtfilt);
*   **Domains:** Array Signal Processing, Direction of Arrival (DOA) Estimation, Acoustic Filtering
