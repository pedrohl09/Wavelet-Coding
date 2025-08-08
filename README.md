# Wavelet-Coding

This repository contains a Python implementation of a digital communication system simulator that uses wavelet codes. This project is a Python translation of the original C code developed by Professor Luiz Gonzaga of the Department of Communication Engineering (DCO) at the Federal University of Rio Grande do Norte (UFRN).

## Description

The simulator models a complete communication chain, including:

* **Wavelet Coding:** The information bits are encoded using a wavelet transform.
* **APK Modulation (MPSK):** The encoded symbols are modulated using Amplitude and Phase-Shift Keying (APK), specifically M-ary Phase-Shift Keying (MPSK).
* **Rayleigh Fading Channel:** The modulated signal is transmitted through a channel with Rayleigh fading, which simulates a realistic wireless communication environment.
* **Demodulation and Decoding:** The received signal is demodulated and decoded to recover the original information bits.

The simulation calculates the Bit Error Rate (BER) for different Signal-to-Noise Ratios (SNR) and plots the results.

## How to Run

1.  **Clone the repository:**
    ```bash
    git clone https://github.com/pedrohl09/Wavelet-Coding.git
    ```
2.  **Navigate to the project directory:**
    ```bash
    cd Wavelet-Coding
    ```
3.  **Run the simulation:**
    ```bash
    python3 wavelet_coding.py
    ```

The simulation will generate two text files (`wav_py_oop.txt` and `Hard2x8Ray_py_oop.txt`) with the simulation results and a PNG image (`grafico_ber_vs_snr.png`) with the BER vs. SNR plot.

## Original Author

The original C implementation and the concepts behind this simulator were developed by **Professor Luiz Gonzaga** from the Department of Communication Engineering (DCO) at the Federal University of Rio Grande do Norte (UFRN). This Python version aims to provide a more modern and accessible implementation of his work.
