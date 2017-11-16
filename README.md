# Prince-


## Data Preparation

In order to achieve great results with Prince, it is necessary to provide training data that resembles Your data as closely as possible:
1. "Reference" Genome

   - If available, it is recommended to use a reference genome that is related to Your genome(s). The reference genome should be modified to include a range of copy numbers. For example, one could replace all the occurances of the templates with the same template, but with a copy number of one, then create an instance with a copy nuber of two, and so on. You can do this using inject_repeats.py.

   - If no reference genome is available, genomes of similar length and with inserted VNTRs can be generated from scratch, or using any other alternative method.

2. Generating reads from the reference

   - There are numerous tools available for generating reads from a reference. We recommend using ART(https://www.niehs.nih.gov/research/resources/software/biostatistics/art/).

   - The coverage and length of reads must match that of Your data.
   
   
 ## Running Prince
 
1. Prince is written in Python 2.7 and uses the biopythin library.

2. Downoald this repo. It contains all the files necessary to run the code. Place your training and actual reads into the same directory.

3. Go to Prince.py. You will need to write down the names of the files containing your training and actual reads, specify the coverage and modify the parameters, if you need to do so.

4. Run the modified Prince.py function.

5. An example of how Prince.py should look like is included above.
   
   





