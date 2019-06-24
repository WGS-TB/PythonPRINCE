# PRINCE - VNTR copy number approximation

PRINCE estimates Variable Number Tandem Repeats (VNTR) copy number from raw next generation sequencing (NGS) data.

## Build status

[![Build Status](https://travis-ci.org/WGS-TB/PythonPRINCE.svg?branch=master)](https://travis-ci.org/WGS-TB/PythonPRINCE)
[![Coverage Status](https://coveralls.io/repos/github/WGS-TB/PythonPRINCE/badge.svg?branch=master)](https://coveralls.io/github/WGS-TB/PythonPRINCE?branch=master)

## Getting Started

### Prerequisites

Python 2.7 or 3.6 with the following packages installed

```
biopython
scipy
numpy
```

### Installing

We recommend installing PRINCE into a conda environment. If you don't already have conda installed, you can download it [here](https://conda.io/miniconda.html). 

```
conda create -n prince biopython scipy
source activate prince
git clone https://github.com/WGS-TB/PythonPRINCE.git
cd PythonPRINCE
pip install -e .
```

You can verify that prince has been installed by typing:

```
prince -h
```

You should see the following output:

```
usage: prince [-h] [-bo BOOST_OUTPUT] [-to TARGET_OUTPUT] [-tmp TEMPLATES]
              [-tf TARGET_FILE] [-bf BOOSTING_FILE] [-k K] [-cn COPYNUMBER]
              [-p PRIMERS] [-np NUM_PROCS]

Prince Options.

optional arguments:
  -h, --help            show this help message and exit
  -bo BOOST_OUTPUT, --boost_output BOOST_OUTPUT
                        output file for training data / training data used to
                        predict copy numbers for queries
  -to TARGET_OUTPUT, --target_output TARGET_OUTPUT
                        output file for query copy number predictions
  -tmp TEMPLATES, --templates TEMPLATES
                        VNTR templates. Default is for M.TB
  -tf TARGET_FILE, --target_file TARGET_FILE
                        target genome names in a text file
  -bf BOOSTING_FILE, --boosting_file BOOSTING_FILE
                        training genome file names in a text file
  -k K, --k K           Kmer size used during read recruitment.
  -p PRIMERS, --primers PRIMERS
                        Flanking sequences used in coverage adjustments
  -np NUM_PROCS, --num_procs NUM_PROCS
                        Number of cores for parallel processing.
                        
```

## Using PRINCE

We recommend using Prince's pre-trained model and settings for querying.

### Querying

To query PRINCE you'll need a target file (eg. samples.txt) with the paths to your fastq files.
Each sample should take up one line. If you are using paired data it is only necessary to specify one file path and name.
Eg. sample_1.fq + sample_2.fq can be specified by sample or sample_.    

samples.txt should look something like this.
```
first_sample
second_sample.fq
third_sample.fastq
sample_folder/fourth_sample_
sample_folder/fifth_sample_1.fastq.gz  sample_folder/fifth_sample_2.fastq.gz
```
Once you have your target file you can run PRINCE.
Specify a target output file (eg. output.txt) with -to. If the file doesn't exist PRINCE will create one. 

```
prince -tf samples.txt -to output.txt
``` 

Each line in output.txt will correspond to the predicted VNTR copy numbers for the corresponding sample in your target file.

### Training

If, for some reason, you'd like to train PRINCE on new data you'll need to do a couple things.

You'll need simulated reads from read simulation software. We used ART.
To use read simulation software you'll need at least one assembled genome of the species of interest.
You'll need to create multiple altered versions of the genome by deleting the VNTR regions of the assembled genome and then inserting various number of copies of the VNTR template.  
In the end you should have at least one altered genome with one copy at each VNTR region, a genome with two copies at each region, a genome with three copies etc.
This should go up to at least 4 but we see no issue with increasing up to 10.
Once you have your altered genomes you can create simulated reads using your preferred software. 

Create a training file (eg. training_samples.txt) with the paths to all your genomes. Each line will correspond to a genome and the number of copies at each VNTR region within that genome, separated by a space.

training_samples.txt should look something like this.
```
genome_w_1_copy 1
genome_w_2_copies 2
sample_folder/genome_w_3_copies 3
another_folder/genome_w_4_copies 4
genome_w_10_copies 10
```

Then to run, specify your training output file with -bo, if the file doesn't exist PRINCE will create one. .
```
prince -bf training_samples.txt -bo training_output.txt
```


To use your new training data on your queries specifiy the training output file.
```
prince -tf samples.txt -to output.txt -bo training_output.txt
```

To accelerate training, specify the number of processes using -np (default is 1).
```
prince -tf samples.txt -to output.txt -bo training_output.txt -np 16
```

## Built With

* [Python 2.7](https://www.python.org/download/releases/2.7/)
* [ART](https://www.niehs.nih.gov/research/resources/software/biostatistics/art/index.cfm)

## Contributing

## Authors

* **Julian Booth**
* **Margaryta Vityaz** 
* **Merhdad Mansouri** 
* **Leonid Chindelevitch** 

## License


## Acknowledgments
