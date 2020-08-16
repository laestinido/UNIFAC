import numpy as np

def calc_ln_gamma_c(x, D):
    
    print('Intermediate Values for ln_\u03B3_c calculation')
    print('-----------------------------------------------')
    #let N be the number of species in the mixture
    N = np.shape(D)[1] - 4

    def calc_r(D):
        r = np.zeros(N)
        for i in range(N):
            r[i] = sum(D[:, 2] * D[:, 4+i])
        return r.round(4)

    def calc_q(D):
        q = np.zeros(N)
        for i in range(N):
            q[i] = sum(D[:, 3] * D[:, 4+i])
        return q.round(4)

    def calc_phi(x, r):
        phi = np.zeros(N)
        for i in range(N):
            phi[i] = (r[i] * x[i]) / np.sum(x * r)
        return phi.round(4)

    def calc_theta(x, q):
        theta = np.zeros(N)
        for i in range(N):
            theta[i] = (q[i] * x[i]) / np.sum(x * q)
        return theta.round(4)

    def calc_l(r, q):
        l = np.zeros(N)
        for i in range(N):
            l[i] = (5 * (r[i] - q[i]) ) - (r[i]-1)
        return l.round(4)
    
    
    r = calc_r(D)
    print('r: {}'.format(r))
    
    q = calc_q(D)
    print('q: {}'.format(q))
    
    phi = calc_phi(x, r)
    print('\u03A6: {}'.format(100*phi))
    
    theta = calc_theta(x, q)
    print('\u03F4: {}'.format(100*theta))
    
    l = calc_l(r, q)
    print('l: {}'.format(l))
    
    ln_gamma_c = np.zeros(N)
    for i in range(N):
        ln_gamma_c[i] = np.log(phi[i] / x[i]) + (5 * q[i] * np.log(theta[i] / phi[i])) + l[i] - ((phi[i] / x[i]) * np.sum(x*l))
    
    print('\u03B3_c: {}\n'.format(ln_gamma_c.round(4)))
    return ln_gamma_c.round(4)
  

#Functions for Solving the Residual Contribution
def calc_ln_gamma_r(x, D, a_mn, T):
    
    
    def calc_ln_F(x, D, a_mn, T, i):
    #use i=0 if relative to the whole mixture
    
        #for superscripting variables during printing
        sup = ['', '\u00B9', '\u00B2', '\u00B3', '\u2074', '\u2075']
        if i==0:
            print('\nln_\u0393 ' + sup[i] + ' calculation values')
        else:
            print('\nln_\u0393\u207D' + sup[i] + '\u207E calculation values')
        
        #K is the number of functional groups
        K = len(D)
    
        def calc_psi(a_mn, T):
            return np.exp(-a_mn/T).round(4)

        def calc_X(x, D, i):
            X = np.zeros(K)
            for k in range(K):
                if i == 0:
                    X[k] = np.sum(x*D[k, 4:]) / np.sum(x * D[:, 4:])
                else:
                    X[k] = np.sum(D[k, i+3]) / np.sum(D[:, i+3])
            return X.round(4)

        def calc_theta(X, D):
            theta = np.zeros(K)
            print(X)
            for k in range(K):
                theta[k] = (D[k, 3] * X[k]) / np.sum(D[:, 3] * X)
            return theta.round(4)

        psi = calc_psi(a_mn, T)
        #print('\u03A8: {}\n'.format(psi))
    
        X = calc_X(x, D, i)
        print('X: {}'.format(X))
    
        theta = calc_theta(X, D)
        print('\u03F4: {}'.format(theta))
    
        ln_F = np.zeros(K)
    
        for k in range(K):
        
            theta_psi_m = 0
        
            for m in range(K):
                theta_psi_n = 0
                for n in range(K):
                    theta_psi_n += theta[n]*psi[n, m]
                theta_psi_m += (theta[m]*psi[k, m]) / theta_psi_n
        
            ln_F[k] = D[k, 3] * (1 - (np.log(np.sum(theta*psi[:, k]))) - theta_psi_m)
        
            #removes value of ln_F[k] in the list if k group is not present in molecule i
            if i != 0 and D[k, i+3] == 0:
                ln_F[k] = 0
    
        if i==0:
            print('ln\u0393: {}'.format(ln_F.round(6)))
        else:
            print('ln_\u0393\u207D' + sup[i] + '\u207E: {}'.format(ln_F.round(6)))
    
        return ln_F.round(6)

    print('Intermediate Values for ln_\u03B3_r calculation')
    print('-----------------------------------------------')
    #let N be the number of species in the mixture
    N = np.shape(D)[1] - 4
    
    ln_gamma_r = np.zeros(N)
    
    ln_F = calc_ln_F(x, D, a_mn, T, 0);
    for i in range(N):
        ln_F_i = calc_ln_F(x, D, a_mn, T, i+1)
        ln_gamma_r[i] = np.sum(D[:, i+4] * (ln_F - ln_F_i) )
    
    print('\u03B3_r: {}\n'.format(ln_gamma_r.round(5)))
    return ln_gamma_r.round(5)


def calc_gamma(x, D, a_mn, T):
    
    N = np.shape(D)[1] - 4
    ln_gamma_c = calc_ln_gamma_c(x, D)
    ln_gamma_r = calc_ln_gamma_r(x, D, a_mn, T)
    
    gamma = np.exp(ln_gamma_c + ln_gamma_r)
    
    print('\nSolved Activity Coefficients')
    print('------------------------------')
    
    sub = ['\u2081', '\u2082', '\u2083', '\u2084', '\u2085', '\u2086', '\u2087', '\u2088', '\u2089']
    for i in range(len(gamma)):
        print('\u03B3' + sub[i] + ': {}'.format(gamma[i].round(3)))
    return gamma.round(4)