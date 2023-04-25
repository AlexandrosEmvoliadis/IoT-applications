**PROJECT FOR DETECTING PARKINSON DISEASE, USING A SMART GASTRICT WEARABLE BELT**

*AUTh Post-Graduate Programme on Advanced Systems of Information and Computers*

EEG_Data Contains Real Data Gathered by the Gastric Belt, before and after having a meal

- n file-notation states for healthy, while pd states for Parkinson Disease
- data gathered from 7 healthy and 7 parkinson disease patients, hence 14 time-series

**FEATURES EXTRACTION**
*Time-series Processing*:
1. Min Max scaling from [0,65535] to [-1,1]
2. Standard Nosmalization
3. Downsampling, due to very slow waves
4. Band-Pass Butterworth IIR Filter (2nd-4th order, 4th is chosen due to more flat response)

*Features extraction*:
1. Mean and Standard Deviation For Dominant Power and Dominant Frequency (DP & DF)
2. Overall Dominant Power and Frequency (ODP & ODF)
3. Dominant Power Ratio (DPR)
4. Power and Frequency Instability Coefficient (FIC & PIC)
5. Percentage of Normal Gastric Slow Waves
6. Mean and Standard Deviation for PEPD
7. Crest Factor (CF)

Additionally computed:
1. Mean and Standard Deviation for Kurtosis and Skewness
2. Kurtosis and Skewness for Power and Frequency Ratio

**HYPOTHESIS TESTING**

Perform statistical tests for any features between Healthy and Parkinson Disease, using:
1. Wilcoxon rank-sum test
2. Mann-Whitney U test

to reveal whether or not these samples are drawn from the same distribution

**FURTHER PROCESSING FOR CLASSIFICATION**

Very small number of samples, containing outliers. Some options:
1. Delete any feature contains outliers (will result to less features)
2. Delete any sample that contains outliers (will result to even smaller dataset)
3. Replace outliers by *IQR* method (Chosen One)
    - Detect Outliers via IQR
    - Replace these values with a random sample from the same Class 

**CLASSIFICATION**

To increase the number of samples, Generative Adversarial Modeling is utilized. The provided scripts utilize *SDV* library for data augmentation, using CopolaGAN and CT-GAN.

Several classifiers tested, using GridSearch Library for parameter-tuning.

Experiments on Decision Trees, trained on synthetic data and tested on real, reveal:
*F1 = 92.3%, Precision = 100%, Recall = 85.7% and Accuracy = 92.86%*, trained on 2000 samples

Experiments with PCA and Vanilla Auto-Encoder, revealed poor performance.

**RUN**

- Make sure that your environment contains svd library
- Make sure that you are inside './Parkinson_detection' folder
- process_data.py will:
  - Perform time-series process
  - Features extraction (original, and original + additional features)
  - Create folders to store the .csv files
- data_augmentation.py will generate synthetic data using CopulaGAN and CT-GAN, for 100,1000 and 10.000 samples.
- augmented_data_tests_and_ml.py will:
  - Perform statistical analysis on the generated data
  - Perform Classification with an SVM*



*Any traditional ML Classifier can be tested, using GridSearch library
