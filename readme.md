Projection Two-Photon Lithography Simulation Toolbox for MATLAB (OS-FDM 2D)
v1.0 (Archived, not actively supported)
Rushil Pingali and Sourabh K. Saha
Georgia Institute of Technology

This toolbox can be used to predict projection two-photon lithography process outcomes.
It contains the following files:

- GettingStarted.mlx: MATLAB live script demonstrating how to use the following functions in sequence to perform a simulation experiment
- make_system_profile.m: lets you save the parameters of your printing system as a .mat file for easy access
- optical_dosage.m: function that returns intensity and optical dosage distributions when given a digital micromirror device (DMD) projection
- STEAM_Lab_GT.mat: sample outputs from optical_dosage.m containing parameters for the authors' printing systems at Georgia Institute of Technology 
- fdm_2D.m: function that performs a finite difference reaction-diffusion simulation to predict resulting degree of conversion (DOC) for a given optical dosage. 2D implementation.

Link to accompanying paper/references: Pingali, R., and Saha, S. K. "Rapid Modeling of Photopolymerization in Projection Two-Photon Lithography Via an Operator Splitting Finite Difference Method." ASME. J. Micro Nano Sci. Eng. March 2024; 12(1): 011001. https://doi.org/10.1115/1.4065706

License: MIT License

License text:


Copyright 2024 Rushil Pingali and Sourabh K. Saha 

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
