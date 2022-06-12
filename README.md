# Welcome to KNSegmentation!

KNSegmentation is a method to segment 12-lead ECG signal into P-wave, QRS complex, and T-wave. For each wave, it identifies the onset, offset, and peak(s) of the wave. For QRS complex, it also identifies individual waves comprising the complex (Q, R, R', S, S').
KNSegmentation is light-weight and based on a recursive divide-and-conqure algorithm called Piecewise Linear Approximation, which reduces the representation of the ECG signal to a sequence of its significant points.
## Future Release
The Matlab section of the code will be migrated to Python shortly and the tool will be expanded to allow users to upload their own ECG.

# Usage

One part of KNSegmentation is implemented using Matlab and the other part using Python (I will migrate everything to Python soon). Before getting into execuation, you first need to load your ECG data.

## Step 0: Input & Output
- **Input:** 12-lead ECGs must be in the **CSV format**, in the dimension of **12 x n**, and stored in **Data/ECG** directory. 
- **Output:** KNSegmentation reads through Data/ECG directory, perfroms segmentation on each ECG, and saves the results in **Data/FinalAnnotation**.

## Step 1: Run Matlab script "PLARunner.mat"

Run PLARunner.mat using Matlab.
> **Requirements:** Matlab 2019 or newer with Signal Processing package installed.


## Step 2: Run Python script "Main.py"

Run Main.py and done!
> **Requirements:** Python 3.8, Numpy, SciPy, Pandas, and MatPlot.

## Demo

I have uploaded three ECGs to demo the software. You just need to run Main.py to see the results.
