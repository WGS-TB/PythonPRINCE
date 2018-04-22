from math import sqrt
from predict import get_data, get_equations, get_copy_number
from match_score import compute_match_score

def test_target(opts, templates, templateKmers):
    with open(opts.target_file) as file:
	query = file.readline().strip("\n")
	while query:

	    targetFileName = query.split("/")[-1] #CHANGE
    	    targetMatchScore = compute_match_score(query, templates, templateKmers,
                                             opts.coarse,opts.fine,opts.screen)

    	    data = get_data(opts.boost_output)
    	    equations = get_equations(data)
    	    predictions = []

    	    # Write target predictions to text file
 
    	    with open(opts.target_output, "a") as f:
        	f.write(targetFileName)
        	for t, ms in enumerate(targetMatchScore):
            	    slope, intercept = equations[t]
            	    y_predict = get_copy_number(ms, slope, intercept)
            	    f.write("," + str(y_predict))

        	f.write("\n")
	    query = file.readline().strip("\n")
