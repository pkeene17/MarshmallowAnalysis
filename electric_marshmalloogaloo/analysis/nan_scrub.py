#Helper that scrubs nan values in a numpy array
#input: data-assumes its a one dimensional numpy array. Otherwise, will return a flattened array
#output: data-but without nans
import numpy as np

def nan_scrub(data):
	return data[~np.isnan(data)]