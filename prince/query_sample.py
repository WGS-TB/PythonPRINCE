import os
from prince.predict import get_data, get_equations, get_copy_number
from prince.match_score import compute_match_score
import time

def test_target(opts, template_obj, primers):
    with open(opts.target_file) as file:
        query = file.readline().strip("\n")
        while query:
            start_time = time.time()
            targetFileName = query.split("\t")[0].split("/")[-1] 
            target_forward = query.split("\t")[0]
            try:
                target_reverse = query.split("\t")[1]
            except:
                target_reverse = ""
            print("\nQuerying %s" % targetFileName)
            
            targetMatchScore = compute_match_score(target_forward, target_reverse, template_obj, opts.k, primers)
            
            data = get_data(opts.boost_output)
            equations = get_equations(data)
            predictions = []
            
            # Write target predictions to text file
            with open(opts.target_output, "a+") as f:
                if os.path.getsize(opts.target_output) == 0:
                    f.write("Templates,")
                    f.write(",".join(template_obj["Names"]))
                    f.write("\n")
                f.write(targetFileName)
                for t, ms in enumerate(targetMatchScore):
                    slope, intercept = equations[t]
                    y_predict = get_copy_number(ms, slope, intercept)
                    f.write("," + "{:.2f}".format(y_predict))
                
                f.write("\n")
            print("Done with %s" % targetFileName)
            print(time.time()-start_time)
            query = file.readline().strip("\n")
