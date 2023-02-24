import pickle
import numpy as np
import pandas as pd

def row_sum_n_col(array_, n_col):
    return np.sum(np.reshape(array_, (-1, int(np.shape(array_)[1] / n_col), n_col)), axis=2)

def generate_sample_MRIOT():
    with open('processed_gtap/regions_list.pkl', 'rb') as f:
        regions = pickle.load(f)
    with open('processed_gtap/activities_list.pkl', 'rb') as f:
        activities = pickle.load(f)

    # Z 141*65 141*65
    # y 141*65 141
    # x 141*65 1
    # v 141*65 1

    Z = np.random.rand(141 * 65, 141 * 65)
    v = np.random.rand(141 * 65)
    v_by_country = row_sum_n_col(np.reshape(v, (1, -1)), 65).flatten()
    y = np.tile([1/(141*65)] * 141, (141*65, 1))
    y = y * v_by_country[None, :]
    MRIOT = np.concatenate([np.concatenate([Z, y, np.reshape(x, (len(x), 1))], axis=1),
                            np.concatenate([np.reshape(v, (1, len(v))),
                                            np.ones((1, 141 + 1)) * np.nan], axis=1),
                            np.concatenate([np.reshape(x, (1, len(x))),
                                            np.ones((1, 141 + 1)) * np.nan], axis=1)],
                           axis=0)
    df = pd.DataFrame(MRIOT)

    indexes = []
    for region in regions:
        for act in activities:
            indexes.append(region + '_' + act)
    df.set_axis(indexes + ['VA', 'X'], axis=0, inplace=True)
    df.set_axis(indexes + ['y_' + region for region in regions] + ['X'], axis=1, inplace=True)
    with open('processed_gtap/adjusted_MRIOT_array.pkl', 'wb') as f:
        pickle.dump(df, f)
