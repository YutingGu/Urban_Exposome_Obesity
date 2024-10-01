# Linking the Urban Exposome to Obesity via Proteome in the UK Biobank
### Overview
This project explores the relationship between urban environmental exposures (the urban exposome) and obesity, using proteomic data from the UK Biobank. The goal is to investigate how different urban exposures contribute to obesity outcomes through biological mediation via the proteome.

### Key Objectives
- Identify **urban expotypes** using a variety of environmental exposures.
- Characterise these expotypes and their biological embodiment through **proteomic signatures**.
- Assess the impact of these expotypes on **obesity** outcomes.

### Dataset
The study uses data from 51,879 individuals and incorporates 15 urban exposomes such as air pollution, proximity to major roads, and access to green spaces. Proteomic data consists of 1,343 proteins, with 43,572 individuals having complete data on obesity-related outcomes (BMI, LDL cholesterol, systolic blood pressure, and HbA1c).

### Methodology
The analysis workflow is broken down into several steps:

- Expotype Identification: Clustering individuals based on urban exposomes using k-means clustering.
- Exposome Characterisation: Logistic regression is applied to assess the contribution of different urban exposures to expotype classification.
- Proteomic Signature Identification: Stability selection using LASSO identifies proteins associated with specific expotypes.
- Obesity Association: Linear and logistic regression models investigate the association between expotypes and obesity-related outcomes.

### Results
1. The analysis identified four distinct expotypes with varying exposure to pollutants, traffic, and natural environments.
2. Each expotype had a unique set of proteomic signatures, with some proteins showing a significant association with obesity outcomes.
3. Expotype 0, for example, exhibited high levels of particulate matter pollution and was associated with a higher risk of obesity.

<img src="https://github.com/user-attachments/assets/dc6da2a1-65a4-4093-aa6b-a0c2690d5446" alt="CompEpi" width="800"/>

### Relevant Implementations
Relevant code for the analysis can be found in the [Scripts](Scripts) directory. The scripts include:

**STEP 0. K-means clustering** for expotype identification in Python.
**STEP 1. Univariate analysis** of the impact of urban exposomes on expotype assignment (Logistic regression)
**STEP 2. Stability LASSO** for protein selection
**STEP 3. Univariate analysis** of expotype-specific effects on obesity outcomes (Linear and Logistic Regression)

Acknowledgement
This project is a team collaboration as part of the Computational Epidemiology module (Teammates: Amin Moghaddam, Hanh Lan Bui, Lucie Frechin, Riya Nagar).
