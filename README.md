# Gaussian UV-Vis Spectrum Generator

A desktop GUI tool (written in Python) designed for analyzing and visualizing electronic spectra (UV-Vis) based on Gaussian output files (TD-DFT).

This application was created for scientific workflows. It allows for quick verification of electronic transitions, monitoring of spin contamination, and generation of publication-ready plots.

## Key Features

* **Data Loading**: Supports standard Gaussian 16 log files (`.log`, `.out`).
* **Visualization**:
    * Interactive spectrum plot (Gaussian convolution) and transition sticks.
    * "Publication-ready" aesthetics (white background, clear axes).
    * Adjustable Full Width at Half Maximum (FWHM/$\sigma$).
* **Data Analysis**:
    * Transition table with filtering by oscillator strength ($f$).
    * **Spin Verification**: Automatic comparison of calculated $<S^2>$ vs. ideal values to flag states with significant spin contamination.
* **Export**:
    * **Excel (.xlsx)**: Exports data to a single file with separate sheets: spectrum data, all transitions, and filtered transitions.
    * **Figures**: Saves plots as PNG, PDF, or SVG (default 600 DPI).

## Requirements

The script requires Python 3 and the following libraries:

* `customtkinter` (GUI)
* `pandas` & `numpy` (data analysis)
* `matplotlib` (plotting)
* `openpyxl` (Excel export support)

## Installation and Usage

1.  Install the required libraries:
    ```bash
    pip install customtkinter pandas matplotlib numpy openpyxl
    ```
2.  Run the script:
    ```bash
    python main_ctk.py
    ```
