import numpy as np

def mst_ransac_tdoa_matrix_to_tdoa_vector(tdoam, n_iterations=100, inlier_threshold=0.1):
    n = tdoam.shape[0]
    v = np.zeros(n)

    best_score = 0

    for ransac_iter in range(n_iterations):
        perm = np.random.permutation(n)
        v[perm[0]] = 0
        for i in range(n-1):
            v[perm[i+1]] = v[perm[i]] + tdoam[perm[i], perm[i+1]]
        prediction = np.expand_dims(v, 0) - np.expand_dims(v, 1)

        res = prediction - tdoam

        score = np.sum(np.abs(res.flatten()) < inlier_threshold)
        if score > best_score:
            best_score = score
            tdoav = v.copy()

    prediction = np.expand_dims(tdoav, 0) - np.expand_dims(tdoav, 1)
    res = prediction - tdoam

    # use least squares to re-estimate based on inliers
    inliers = np.where(np.triu(np.abs(res) < inlier_threshold, 1))

    M = np.zeros((len(inliers[0]), n-1))
    b = np.zeros((len(inliers[0]), 1))
    counter = 0
    for i in range(len(inliers[0])):
        if inliers[0][i] >= inliers[1][i]:
            continue

        if inliers[0][i] == 0:
            M[counter, inliers[1][i]-1] = 1
            b[counter] = tdoam[0, inliers[1][i]]
        else:
            M[counter, inliers[0][i]-1] = -1
            M[counter, inliers[1][i]-1] = 1
            b[counter] = tdoam[inliers[0][i], inliers[1][i]]
        counter += 1
    tdoav[0] = 0
    tdoav[1:] = np.linalg.lstsq(M, b, rcond=None)[0].squeeze(1)

    prediction = np.expand_dims(tdoav, 0) - np.expand_dims(tdoav, 1)
    res = prediction - tdoam

    score = sum((np.abs(res) < inlier_threshold).flatten())

    return tdoav, score