"""
Calculate a set of metrics, then calculate the weighted sum for the specified
metrics for a given input map.
"""


def ws_3_metrics(v1, v2, v3, w1, w2, w3, r1, r2, r3, s1, s2, s3):
    # Determine the scaled values (range of 0 - 1)
    scaled_value_1 = (v1 - min(r1))/(max(r1) - min(r1))
    scaled_value_2 = (v2 - min(r2)) / (max(r2) - min(r2))
    scaled_value_3 = (v3 - min(r3)) / (max(r3) - min(r3))
    # Initialise the total weighted sum.
    total = 0
    # If want to maximise the first objective, add to total.
    if s1 == "maximise":
        total = total + w1 * scaled_value_1
    # Otherwise, subtract from total.
    elif s1 == "minimise":
        total = total - w1 * scaled_value_1
    # If want to maximise the second objective, add to total.
    if s2 == "maximise":
        total = total + w2 * scaled_value_2
    # Otherwise, subtract from total.
    elif s2 == "minimise":
        total = total - w2 * scaled_value_2
    # If want to maximise the third objective, add to total.
    if s3 == "maximise":
        total = total + w3 * scaled_value_3
    # Otherwise, subtract from total.
    elif s3 == "minimise":
        total = total - w3 * scaled_value_3
    # Return to weighted total sum.
    return total


def ws_2_metrics(v1, v2, w1, w2, r1, r2, s1, s2):
    # Determine the scaled values (range of 0 - 1)
    scaled_value_1 = (v1 - min(r1)) / (max(r1) - min(r1))
    scaled_value_2 = (v2 - min(r2)) / (max(r2) - min(r2))
    # Initialise the total weighted sum.
    total = 0
    # If want to maximise the first objective, add to total.
    if s1 == "maximise":
        total = total + w1 * scaled_value_1
    # Otherwise, subtract from total.
    elif s1 == "minimise":
        total = total - w1 * scaled_value_1
    # If want to maximise the second objective, add to total.
    if s2 == "maximise":
        total = total + w2 * scaled_value_2
    # Otherwise, subtract from total.
    elif s2 == "minimise":
        total = total - w2 * scaled_value_2
    # Return to weighted total sum.
    return total


def ws(metrics, no_metrics):
    # Determine the weighted sum for a set of input metrics.
    # Metrics is a dictionary that stores a list of information corresponding
    # to each metrics used. The form is:
    # metrics[index #][value, weight, min_value, max_value, form (0 = min, 1 = max)]
    # Initialise the total weighted sum.
    total = 0
    for i in range(0, no_metrics):
        scaled_value = (metrics[i][0] - metrics[i][2]) / (metrics[i][3] -
                                                          metrics[i][2])
        if metrics[i][4] == "maximise":
            total = total + scaled_value * metrics[i][1]
        elif metrics[i][4] == "minimise":
            total = total - scaled_value * metrics[i][1]
    return total
# Module test.
"""
metrics = {
    0: [0.89, 0.33, 0.85, 1, "maximise"],
    1: [0.05, 0.33, 0, 0.1, "maximise"],
    2: [0.003, 0.5, 0, 0.004, "minimise"]
}
no_metrics = 3

weighted_sum = ws(metrics, no_metrics)
print(weighted_sum)
"""