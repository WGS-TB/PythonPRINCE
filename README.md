# Prince-


## Data Preparation

In order to achieve great results with Prince, it is necessary to provide training data that resembles Your data as closely as possible:
1. "Reference" Genome

   - If available, it is recommended to use a reference genome that is related to Your genome(s), however it is neccessary for it to contain a range of copy numbers within it. 

   - If no reference genome is available, genomes of similar length and with inserted VNTRs can be generated from scratch, or using any other alternative method.

2. Generating reads from the reference

   - There are numerous tools available for generating reads from a reference. We recommend using ART(https://www.niehs.nih.gov/research/resources/software/biostatistics/art/).

   - The coverage and length of reads must match that of Your data.
   
   
 ## Running Prince
 
1. Prince is written in Python 2.7. You'll need that.

2. Downoald this repo. It contains all the files necessary to run the code. Place your training and actual reads into the same directory.

3. Go to Prince.py. There you will need to write down the names of the files containing your training and actual reads, specify the coverage and maybe play with the parameters, if you need to do so.

4. Run the modified Prince.py function.

5. An example of how Prince.py should look like is included above.
   
   





