import os
from prince.predict import get_data, get_equations, get_copy_number
from prince.match_score import compute_match_score
import time
         
import multiprocessing as mp

def test_target(opts, template_obj, primers):
    # Get the query paths
    with open(opts.target_file) as file:
            queries = [line.rstrip("\n") for line in file]
    
    pool = mp.Pool(processes=opts.num_procs)
    # Run analyses in multiple processes 
    results = [pool.apply_async(compute_match_score,(query, template_obj, opts.k, primers)) 
                   for query in queries]
    match_scores = [result.get() for result in results]
    
    # Write results
    data = get_data(opts.boost_output)
    equations = get_equations(data)
    with open(opts.target_output,'a+') as file:
        if os.path.getsize(opts.target_output) == 0:
            file.write("Templates,")
            file.write(",".join(template_obj["Names"]))
            file.write('\n') 
        for (targetMatchScore, targetFileName) in match_scores:
            file.write(targetFileName)
            for t, ms in enumerate(targetMatchScore):
                slope, intercept = equations[t]
                y_predict = get_copy_number(ms, slope, intercept)
                file.write("," + "{:.2f}".format(y_predict))
            file.write('\n')
