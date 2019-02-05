from prince.match_score import compute_match_score

import multiprocessing as mp

def run_boosts(opts,template_obj, primers):
    with open(opts.boosting_file) as file:
        samples = [line.rstrip("\n") for line in file]
        
    pool = mp.Pool(processes=opts.num_procs)
    # Run analyses in multiple processes 
    results = [pool.apply_async(compute_match_score,(sample, template_obj, opts.k, primers)) 
                   for sample in samples]
    match_score_vectors = [result.get() for result in results]
    
    with open(opts.boost_output, "a") as f:
        for matchscore_vector, file_name in match_score_vectors:
            for loci_num, matchscore in enumerate(matchscore_vector): 
                f.write(file_name + "," + template_obj["Names"][loci_num] + "," + str(loci_num) + "," + str(opts.copynumber) + "," + str(matchscore) + "\n")
