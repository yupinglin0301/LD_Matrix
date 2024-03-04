

## 1. Project Overview

This (sub-)repository contains code related to our manuscript " ".

This repository stores data and analysis modules to create linkage disequilibrium (LD) block and conduct LD dimension reduction on genotype data.


## 2. File Structure & Description
~~~~~~~
|-- data
    |-- geno.bed # Binary data file for genotype
    |-- geno.bim # Binary data file for SNP info
    |-- geno.fam # Binary data file for individual information
|-- results
    |-- DisasterResponse.db # database to save clean data
|-- scripts
    |-- classifier.pkl # saved model
    |-- train_classifier.py # machine learning script that creates and trains a classifier, and stores the classifier into a pickle file
|-- Log
|-- README.md
|-- requirements.txt # list of necessary python packages
~~~~~~~

### Reference Data
The reference panel based on the 1000 Genomes Project is available here:
+ https://www.cog-genomics.org/plink/1.9/resources




## 3. Instructions for getting started
### 3.1. Cloning
To run the code locally, create a copy of this GitHub repository by running the following code in terminal:
```sh
git clone https://github.com/timchansdp/Disaster-Response-Pipelines.git
```

### 3.2. Dependencies
The code is developed with Python 3.9.1 and is dependent on python packages listed in `requirements.txt`. To install required packages, run the following command in the project's root directory:
```sh
pip install -r requirements.txt
```
