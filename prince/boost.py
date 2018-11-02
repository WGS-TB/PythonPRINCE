from prince.match_score import compute_match_score

def run_boosts(opts,templates,templateNames,templateKmers):
    with open(opts.boosting_file) as file:
        sample = file.readline().strip("\n")
        while sample:
            sampleFileName = sample.split("\t")[0].split("/")[-1] 
            sample_forward = sample.split("\t")[0]
            try:
                sample_reverse = sample.split("\t")[1]
            except:
                sample_reverse = ""
            print("\nTraining with %s" % sampleFileName)
            boostingMatchScore = compute_match_score(sample_forward, sample_reverse,templates, templateKmers, opts.k)
            with open(opts.boost_output, "a") as f:
                for t_num, matchscore in enumerate(boostingMatchScore):
                    f.write(sampleFileName + "," + templateNames[t_num] + "," + str(t_num) + "," + str(opts.copynumber) + "," + str(
                        matchscore) + "\n")
            print("Done with %s" % sampleFileName)
            sample = file.readline().strip("\n")
