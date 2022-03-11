#!/usr/bin/python3

import numpy as np

def block_average(prop, av_prop):

    sigmas_mean = []
    errs_sigma_mean = []
    ks = []

    k=1
    nb = int((len(prop))/k)

    sum_a = np.zeros(nb)

    for i in range(nb):
        tmp = 0
        for j in range(k):
            tmp += prop[k*i+j]
        sum_a[i] = tmp

    new_a_0 = (1/k)*sum_a
    a_0 = new_a_0

    sum_sigma = 0

    for i in range(0, len(a_0), 1):
        sum_sigma += (new_a_0[i]-av_prop)**2

    sigma2_a = (1/(nb))*sum_sigma

    sigma2_mean = sigma2_a/(nb-1)
    sigma_mean = np.sqrt(sigma2_a)/np.sqrt((nb-1))
    err_sigma = 1/(np.sqrt(2*(nb-1)))

    ks.append(np.sqrt(k))
    sigmas_mean.append(sigma_mean)
    errs_sigma_mean.append(err_sigma*sigma_mean)

    k = k*2

    while (k<len(prop)+1):

        nb = int((len(prop))/k)
        sum_a = np.zeros(nb)

        for i in range(0, int(len(a_0))-1, 2):
            index = int(np.floor(i/2))
            sum_a[index] += a_0[i]
            sum_a[index] += a_0[i+1]
            sum_a[index] = sum_a[index]/2

        a_0 = sum_a

        sum_sigma = 0

        for i in range(0, len(a_0), 1):
            sum_sigma += (a_0[i]-av_prop)**2

        if (nb>1):
            sigma2_a = (1/(nb))*sum_sigma

            sigma2_mean = sigma2_a/(nb-1)
            sigma_mean = np.sqrt(sigma2_a)/np.sqrt((nb-1))
            err_sigma = 1/(np.sqrt(2*(nb-1)))

        else:

            sigma_mean = 0
            err_sigma = 0

        ks.append(np.sqrt(k))
        sigmas_mean.append(sigma_mean)
        errs_sigma_mean.append(err_sigma*sigma_mean)

        k = k*2

    return sigmas_mean, errs_sigma_mean

def get_error(prop, av_prop, k):

    nb = int((len(prop))/k)

    sum_a = np.zeros(nb)

    for i in range(nb):
        tmp = 0
        for j in range(k):
            tmp += prop[k*i+j]
        sum_a[i] = tmp

    new_a_0 = (1/k)*sum_a

    sum_sigma = 0

    for i in range(0, len(new_a_0), 1):
        sum_sigma += (new_a_0[i]-av_prop)**2

    sigma2_a = (1/(nb))*sum_sigma

    sigma2_mean = sigma2_a/(nb-1)
    sigma_mean = np.sqrt(sigma2_a)/np.sqrt((nb-1))

    return sigma_mean
