# PRINCE - VNTR copy number approximation

PRINCE estimates Variable Number Tandem Repeats (VNTR) copy number from raw next generation sequencing (NGS) data.

## Getting Started

### Prerequisites

Python 2.7 with the following packages installed

```
biopython
scipy
numpy
```

### Installing

```
git clone https://github.com/WGS-TB/PythonPRINCE.git
```

To check PRINCE is installed properly run 

```
prince -to test_output.txt -tf sample_targets.txt
```
The output file should contain two rows with 24 random real numbers.

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

Create a separate training file for each copy number with the paths to all your genomes with that many copies at each VNTR region.
Specify your training output file.
```
prince -bf training_samples_cn_1.txt -bo training_output.txt -cn 1
prince -bf training_samples_cn_2.txt -bo training_output.txt -cn 2
prince -bf training_samples_cn_3.txt -bo training_output.txt -cn 3
prince -bf training_samples_cn_4.txt -bo training_output.txt -cn 4
```
To use your new training data on your queries specifiy the training output file.
```
prince -tf samples.txt -to output.txt -bo training_output.txt
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
