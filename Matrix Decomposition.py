import numpy as np

def PLU_decomposition(A,verbose=False):
    m=A.shape[0]
    n=A.shape[1]
    P=np.eye(n)
    L=np.zeros([n,n])
    U=np.array(A).copy()
    if m!=n:
        print "This matrix cannot realize PLU_decomposition because it is not a square matrix!"
        return
    else:
        for i in range(n-1):
            for j in range(i+1,n):
                for k in range(n-1,i,-1):
                    if U[i,i] < U[k,i]:
                        L[i, :], L[k, :] = np.array([L[k, :], L[i, :]]).copy()
                        P[i, :], P[k, :] = np.array([P[k, :], P[i, :]]).copy()
                        U[i, :], U[k, :] = np.array([U[k, :], U[i, :]]).copy()
                if U[i,i]==0:
                    L[j,i]=0
                else:
                    L[j,i]=U[j,i]/float(U[i,i])
                    U[j,:]=U[j,:]-L[j,i]*U[i,:]

        for i in range(n): L[i,i]=1

        if verbose :
            print"The PLU_decomposition decomposition of A is:"
            print "The P matrix is:"
            print (P)
            print "The L matrix is"

            print (L)
            print "The U matrix is"
            print (U)

        return P,L,U


def Gram_Schimidt(A,verbose=False):
    m = A.shape[0]
    n = A.shape[1]
    R=np.zeros([m,n])
    Q=np.zeros([m,m])
    R[0,0]=np.linalg.norm(A[:,0])
    if R[0, 0] == 0:
        print"cannot realiaze Gram_Schimidt decomposition"
    else:
        Q[:,0]=A[:,0]/R[0,0]
    for i in range(1,n):
        Q[:,i]=A[:,i]
        for j in range(i):
            R[j,i]=np.dot(Q[:,j].T,Q[:,i])
            Q[:,i]-=R[j,i]*Q[:,j]
        R[i,i]=np.linalg.norm(Q[:,i])
        if R[i,i]==0:
            print"cannot realiaze Gram_Schimidt decomposition"
        else:
            Q[:, i] = Q[:, i] / R[i, i]
    if verbose:
        print"The Gram_Schimidt decomposition of A is:"
        print ("The Q matrix is:")
        print (Q)
        print ("The R matrix is:")
        print (R)
    return  Q,R

def householder_reduction(A,verbose=False):
    m = A.shape[0]
    n = A.shape[1]
    E=np.eye(m)
    T=np.array(A).copy()
    P=np.eye(m)
    for i in range(n-1):
        u=T[i:,i]-np.linalg.norm(T[i:n,i])*E[i:,i]
        u=np.array([u])
        u=u.T
        I = np.eye(m)
        I[i:,i:]=np.eye(n-i)-2*np.dot(u,u.T)/(np.linalg.norm(u)*np.linalg.norm(u))
        T=np.dot(I,T)
        P=np.dot(I,P)
        Q=P.T
    if verbose:
        print "The householder_reduction of A is:"
        print "The Q matrix is:"
        print Q
        print "The R matrix is:"
        print T
    return  Q,T

def Givens_reduction(A,verbose=False):
    m = A.shape[0]
    n = A.shape[1]
    P=np.eye(m)
    R=np.array(A).copy()
    for k in range(m):
        for j in range(n-1,k,-1):
            tmpa = R[k,k]
            tmpb = R[j,k]
            mag = np.sqrt(tmpa*tmpa+tmpb*tmpb)
            c= tmpa / mag
            s= tmpb / mag
            I=np.eye(m)
            I[k,k]=c
            I[k,j]=s
            I[j,j]=c
            I[j,k]=-s
            P=np.dot(I,P)
            R=np.dot(I,R)
            Q=P.T
    if verbose:
        print"The Givens_reduction of A is:"
        print ("The Q matrix is:")
        print (Q)
        print ("The R matrix is:")
        print (R)
    return  Q,R
# A=np.array([[1 , 2 , -3 , 4],[4 , 8 , 12 , -8],[2 , 3 , 2 , 1 ],[-3 , -1 , 1 , -4 ]])
# A = np.array([[1, 19, -34], [-2, -5, 20], [2, 8, 37]])
#A = np.array([[0, -20, -14], [3, 27, -4], [4, 11, -2]])
# A = np.array([[3, 2, 1], [2, -3, 4], [5, 1, -1], [7, 4, 2]])
# A=np.array([[0,0],[0,0]])
if __name__=='__main__':
    #A = np.random.randn(9,9)
    A = np.random.randint(0,10000,[4,4])
    print"the A matrix is:"
    print A
    print "please choose the way to decompose A . "
    print  "Press '1' to realize the PLU_decomposition of A. "
    print  "Press '2' to realize the Gram_Schimidt of A."
    print  "Press '3' to realize the householder_reduction of A."
    print  "Press '4' to realize the Givens_reduction of A."
    Q = np.empty(A.shape)
    R = np.empty(A.shape)
    choose_mode=input("please choose the mode:")
    if choose_mode==1:
        P, L, U=PLU_decomposition(A)
        print  np.sum(np.linalg.norm(np.dot(P,A) - np.dot(L,U)))
        print "The P matrix is:"
        print (P)
        print "The L matrix is"

        print (L)
        print "The U matrix is"
        print (U)

    elif choose_mode==2:
        Q, R  = Gram_Schimidt(A)
        print  np.sum(np.linalg.norm(A - np.dot(Q, R)))

    elif choose_mode==3:
        Q, R = householder_reduction(A)
        print  np.sum(np.linalg.norm(A - np.dot(Q, R)))
    elif choose_mode==4:
        Q,R = Givens_reduction(A)
        print  np.sum(np.linalg.norm(A - np.dot(Q, R)))

    if choose_mode is not 1:
        print ("The Q matrix is:")
        print (Q)
        print ("The R matrix is:")
        print (R)

        print "QR decomposition in numpy is \n", np.linalg.qr(A)