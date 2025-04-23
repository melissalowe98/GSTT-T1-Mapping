The following MATLAB (R2023B) script was designed to perform T1-fitting on MOLLI data. 

The script has the following steps:

**MASKING**
1. Load the DICOM data: The images are read and sorted based on inversion time, and stored in a 4D array
(rows × cols × slices × inversion times).
2. Create ROI mask: The user selects an ROI from the first slice and first inversion time using an interactive
tool createMask(imellipse). A binary mask is created and expanded across all slices and inversion times to
create a 4D mask.
3. Apply ROI mask: Mask is applied to the image data, removing unwanted background pixels.

**REGISTRATION**
1. Choose reference image: The reference image is chosen to give the best contrast, typically the 5th or 6th
inversion time.
2. Sequential registration: Each subsequent inversion time image is aligned to the previously registered image
using imregister() with a rigid transformation. This is performed forwards and backwards.

**T1 FITTING**
1. Create 2D Matrix: A 2D matrix is created, where each column represents the signal intensities of a single
voxel across all inversion times.
2. Fitting: A nonlinear least squares curve-fitting approach is used, applying the MOLLI-specific equation
S(t) = |A − Be^(−t/T1⋆)| where A and B are fitting coefficients, and T1⋆ is the apparent relaxation time before
Look-Locker correction.
3. Look-locker correction: The final T1 value is calculated using a look-locker correction factor T1 = T1⋆(A/B − 1).
4. Outlier rejection: Remove any outlier values e.g., negative T1, extremely high values.
5. Standard deviation (SD) computation: The residuals between the fitted model and measured signal intensities
are analysed. The SD of the fitted T1 values within the selected ROI is computed to quantify variability and
assess measurement precision. This provides insight into the reliability of the computed T1 values.
6. T1 Maps: The computed voxel-wise T1 values are reshaped into a T1 map using the ‘imagine’ tool This
provides a quantitative visualisation of tissue properties.


