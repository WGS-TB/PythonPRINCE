from math import sqrt
from prince.predict import get_data, get_equations, get_copy_number
from prince.match_score import compute_match_score
import time
import multiprocessing as mp

def partition(lst, n):
    division = len(lst) / float(n)
    return [ lst[int(round(division * i)): int(round(division * (i + 1)))] for i in range(n) ]

def multiple_targetMatchScore(opts, queries, templates, templateKmers):
    match_scores = []
    for query in queries:
        targetFileName = query.split("/")[-1]
        print("Querying %s" % targetFileName)
        start_time = time.time()
        match_score = compute_match_score(query, templates, templateKmers, opts.k) 
        match_scores.append( (targetFileName,match_score) )
        print("Done with %s in time %s" % (targetFileName,str(time.time()-start_time)))
    return match_scores

def test_target(opts, templates,templateNames, templateKmers):
    # Get the query paths
    queries = []
    with open(opts.target_file) as file:
        for line in file:
            queries.append( line.rstrip("\n") )
    # Find match scores
    if opts.num_procs > 1:
        queries_partition = partition(queries,opts.num_procs)
        pool = mp.Pool(processes=opts.num_procs)
        # Run analyses in multiple processes 
        results = [pool.apply_async(multiple_targetMatchScore,(opts,queries,templates,templateKmers)) 
                   for queries in queries_partition]
        match_scores_list = [p.get() for p in results]
        match_scores = []
        for lst in match_scores_list:
            match_scores += lst 
    else:
        match_scores = multiple_targetMatchScore(opts,queries,templates,templateKmers)
    # Write results
    data = get_data(opts.boost_output)
    equations = get_equations(data)
    with open(opts.target_output,'w') as file:
        file.write("Templates,")
        file.write(",".join(templateNames))
        file.write('\n') 
        for (targetFileName, targetMatchScore) in match_scores:
            file.write(targetFileName)
            for t, ms in enumerate(targetMatchScore):
                slope, intercept = equations[t]
                y_predict = get_copy_number(ms, slope, intercept)
                file.write("," + "{:.2f}".format(y_predict))
        file.write('\n')
