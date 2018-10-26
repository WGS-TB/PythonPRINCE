from math import sqrt
from prince.predict import get_data, get_equations, get_copy_number
from prince.match_score import compute_match_score
import time

def test_target(opts, templates,templateNames, templateKmers):
    with open(opts.target_file) as file:
        query = file.readline().strip("\n")
        while query:
            start_time = time.time()
            #add split on tab to deal with second column scenarios
            targetFileName = query.split("\t")[0].split("/")[-1] #CHANGE
            print("\nQuerying %s" % targetFileName)
            query1 = query.split("\t")[0]
            try:
                query2 = query.split("\t")[1]
            except:
                query2 = ""
            targetMatchScore = compute_match_score(query1, query2, templates, templateKmers, opts.k)
            
            data = get_data(opts.boost_output)
            equations = get_equations(data)
            predictions = []
            
            # Write target predictions to text file
            with open(opts.target_output, "a+") as f:
                if f.readline() == "":
                    f.write("Templates,")
                    f.write(",".join(templateNames))
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
