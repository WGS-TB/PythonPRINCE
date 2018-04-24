from match_score import compute_match_score

def run_boosts(opts,templates,templateNames,templateKmers):
    import time
    start_time = time.time()
    with open(opts.boosting_file) as file:
        sample = file.readline().strip("\n")
        while sample:

            sampleFileName = sample.split("/")[-1]  
	    print("\nTraining with %s" % sampleFileName)                
            boostingMatchScore = compute_match_score(sample, templates, templateKmers, opts.coarse, opts.fine, opts.screen)
       	    print(time.time() - start_time)
            with open(opts.boost_output, "a") as f:
                for t_num, matchscore in enumerate(boostingMatchScore):
                    f.write(sampleFileName + "," + templateNames[t_num] + "," + str(t_num) + "," + str(opts.copynumber) + "," + str(
                        matchscore) + "\n")

            sample = file.readline().strip("\n")
