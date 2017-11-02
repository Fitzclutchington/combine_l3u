import numpy as np

def remap(int[:,:] ii, int[:,:] jj, float[:,:] H, float[:,:] T, float[:,:] S,
          float[:,:] t_utc, float[:,:] sza_interp ):
    cdef int i, j
    cdef int N = ii.shape[0]
    cdef int M = ii.shape[1]

    for i in range(N):
        for j in range(M):
            H[ii[i,j],jj[i,j]] = H[ii[i,j],jj[i,j]] + 1
            T[ii[i,j],jj[i,j]] = T[ii[i,j],jj[i,j]] + t_utc[i,j]
            S[ii[i,j],jj[i,j]] = S[ii[i,j],jj[i,j]] + sza_interp[i,j]

def smoothing(char[:,:] QL, float[:,:] S, float[:,:] T, float[:,:] H, 
              double[:,:] l3u_sza, float[:,:] l3u_time, double[:,:] f, double f_sum):

    cdef int w = 5
    cdef int n, m, i, j, ind_y, ind_x
    cdef int N = S.shape[0]
    cdef int M = S.shape[1]
    cdef int h_sum = 0
    cdef int win_size =  (w*2 + 1)**2
    cdef int win_length = 2*w + 1
    cdef int count = 0
    cdef double numerator_s =0
    cdef double numerator_t = 0

    for n in range(N):
        for m in range(M):
            if QL[n,m] >= 0:
                
                # count h
                for i in range(-w,w+1):
                    for j in range(-w,w+1):
                        ind_y = n+i
                        ind_x = m+j
                        if ind_x >= 0 and ind_x < M and ind_y >= 0 and ind_y < N:
                            if H[ind_y,ind_x] > 0:
                                h_sum += 1

                if h_sum == win_size:
                    #print h_sum
                    for i in range(0,win_length):
                        for j in range(0,win_length):
                            # print i,j, n-w+i, m-w+j
                            ind_y = n-w+i
                            ind_x = m-w+j
                            if ind_x >= 0 and ind_x < M and ind_y >= 0 and ind_y < N:
                                numerator_s += (f[i,j]*S[ind_y, ind_x]/H[ind_y, ind_x])
                                numerator_t += (f[i,j]*T[ind_y, ind_x]/H[ind_y, ind_x])
                    l3u_sza[n,m] = numerator_s/f_sum
                    l3u_time[n,m] = numerator_t/f_sum
                    numerator_t = 0
                    numerator_s = 0
                    h_sum = 0

                elif h_sum > 0:
                    
                    for i in range(-w,w+1):
                        for j in range(-w,w+1):                            
                            ind_y = n+i
                            ind_x = m+j
                            if ind_x >= 0 and ind_x < M and ind_y >= 0 and ind_y < N:
                                if H[n+i, m+j] > 0:
                                    numerator_s += (S[ind_y, ind_x]/H[ind_y, ind_x])
                                    numerator_t += (T[ind_y, ind_x]/H[ind_y, ind_x])
                                    count +=1
                    l3u_sza[n,m] = numerator_s/count
                    l3u_time[n,m] = numerator_t/count
                    numerator_t = 0
                    numerator_s = 0
                    h_sum = 0
                    count = 0    
                else:
                    numerator_t = 0
                    numerator_s = 0
                    h_sum = 0   


