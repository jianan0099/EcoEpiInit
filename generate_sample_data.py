import pickle
import numpy as np
import pandas as pd


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
    y = np.random.randint(10, 20, size=(141 * 65, 141))
    x = np.sum(Z, axis=1) + np.sum(y, axis=1)
    v = x - np.sum(Z, axis=0)
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
