from math import sqrt
from predict import get_data, get_equations, get_copy_number
from match_score import compute_match_score

def test_target(target,target_output_file, boost_file, templates, templateKmers,filteringKmerLength, matchingKmerLength, f0):

    targetMatchScore = compute_match_score(target, templates, templateKmers,
                                             filteringKmerLength, matchingKmerLength, f0)[0]

    data = get_data(boost_file)
    equations = get_equations(data)
    predictions = []

    # Write target predictions to text file
 
    with open(target_output_file, "a") as f:
        f.write(target[0])
        for t, ms in enumerate(targetMatchScore):
            slope, intercept = equations[t]
            y_predict = get_copy_number(ms, slope, intercept)
            f.write("," + str(y_predict))

        f.write("\n")
