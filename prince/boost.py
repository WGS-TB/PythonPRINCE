from prince.match_score import compute_match_score

import multiprocessing as mp

def run_boosts(opts, template_obj, primers):
    output = ""
    
    with open(opts.boosting_file) as file:
        samples = [line.rstrip("\n").split() for line in file]
        
    pool = mp.Pool(processes=opts.num_procs)
    # Run analyses in multiple processes 
    results = [pool.apply_async(compute_match_score,(sample, template_obj, opts.k, primers, cn)) 
                   for sample, cn in samples]

    match_score_vectors = [result.get() for result in results]
    
    
    for line_number,(matchscore_vector, file_name, cn) in enumerate(match_score_vectors):
        for loci_num, matchscore in enumerate(matchscore_vector): 
            output += file_name + "," + template_obj["Names"][loci_num] + "," + str(loci_num) + "," + str(cn) + "," + str(matchscore) + "\n"
    return output