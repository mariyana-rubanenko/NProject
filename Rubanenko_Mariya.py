import numpy as np
import math 
import matplotlib.pyplot as plt
p0 = np.array([6,1])
p1 = np.array([2,6])
p2 = np.array([-4,3])
p3 = np.array([-3,-1])
p4 = np.array([2,-2])
list_poly = [p0, p1, p2, p3, p4]

def Length (A,B):
    return np.sqrt((B[0]-A[0])*(B[0]-A[0])+(B[1]-A[1])*(B[1]-A[1]))
    
def Determinante(A,B,C,D):
    return A*D-B*C
    
def AreaTriangle (A,B,C):
    return np.fabs(1/2*((B[0]-A[0])*(C[1]-A[1])-(B[1]-A[1])*(C[0]-A[0])))

def Intersect(A1,A2,B1,B2):
    l1 = A1[1] - A2[1]
    b1 = A2[0] - A1[0]
    g1 = Determinante(A1[0],A1[1],A2[0],A2[1])
    l2 = B1[1] - B2[1]
    b2 = B2[0] - B1[0]
    g2 = Determinante(B1[0],B1[1],B2[0],B2[1])
    
    flag = 0
    x = A1[0]
    y = A1[1]
    if (Determinante(l1,l2,b1,b2)!= 0):
       x = - Determinante(g1,g2,b1,b2)/Determinante(l1,l2,b1,b2)
       y = - Determinante(l1,l2,g1,g2)/Determinante(l1,l2,b1,b2)
       flag = 1
    return [x,y,flag]

def Indexes(poly, i, k):
    if (i >= len(poly)):
        ii = i % len(poly)
    else:
        ii = i
    if (ii == 0):
        m_i = len(poly) - 1
    else:
        m_i = ii - 1
            
    if (k >= len(poly)):
        kk = k % len(poly)
    else:
        kk = k
    if (kk == (len(poly) - 1)):
        s_k = 0
    else:
        s_k = kk + 1
    return [m_i, ii, kk, s_k]

def AreaSigma(poly, i, k):
    sum_area = 0
    mi = Indexes(poly, i, k)[0]
    ii = Indexes(poly, i, k)[1]
    kk = Indexes(poly, i, k)[2]
    sk = Indexes(poly, i, k)[3]
    I = [Intersect(poly[mi],poly[ii],poly[kk],poly[sk])[0], Intersect(poly[mi],poly[ii],poly[kk],poly[sk])[1]]
    if (kk < ii):
        for count in range(ii, len(poly) - 1):
            sum_area = sum_area + AreaTriangle (poly[count],poly[count + 1],I)
        sum_area = sum_area + AreaTriangle (poly[len(poly) - 1],poly[0],I)
        for count in range(0, kk):
            sum_area = sum_area + AreaTriangle (poly[count],poly[count + 1],I)
    else:
        if (ii != kk):    
            for count in range(ii, kk):
                sum_area = sum_area + AreaTriangle (poly[count],poly[count + 1],I)
        else:
            sum_area = 0
    return sum_area
      
def EndPts (poly, i, k, sigma, lamba): 
    num = len(poly) - 1
    if ((i == 0) and (k != num)):
        I = Intersect(poly[0],poly[num],poly[k],poly[k + 1])
        IP_mi = Length(I, poly[num])
        IP_sk = Length(I, poly[k + 1])
        area_mi_k = AreaTriangle(I, poly[num], poly[k])
        area_i_sk = AreaTriangle(I, poly[i], poly[k + 1])

    if ((i != 0) and (k == num)) :
        I = Intersect(poly[i],poly[i-1],poly[num],poly[0])
        IP_mi = Length(I, poly[i-1])
        IP_sk = Length(I, poly[0])
        area_mi_k = AreaTriangle(I, poly[i-1], poly[k])
        area_i_sk = AreaTriangle(I, poly[i], poly[0])

    if ((i != 0) and (k != num)) :
        I = Intersect(poly[i],poly[i - 1],poly[k],poly[k + 1])
        IP_mi = Length(I, poly[i - 1])
        IP_sk = Length(I, poly[k + 1])
        area_mi_k = AreaTriangle(I, poly[i - 1], poly[k])
        area_i_sk = AreaTriangle(I, poly[i], poly[k + 1])
        
    T = [1,1]
    if ((i == 0) and (k == num)):
        T[0] = 0
        T[1] = 0
    else:
        sigma_ik = AreaSigma(poly, i, k)
        IP_i = Length(I, poly[i])
        IP_k = Length(I, poly[k])

    if (area_mi_k > sigma + sigma_ik):
        t1 = math.log((IP_k)/(2*lamba*IP_sk))
    else:
        t1 = math.log(2*lamba)
    if (area_i_sk > sigma+sigma_ik):
        t2 = -math.log((IP_i)/(2*lamba*IP_mi))
    else:
        t2 = -math.log(2*lamba)
    return [t1,t2]

def Straight(A, B, t):
    x_t = A[0] + t * (B[0] - A[0])
    y_t = A[1] + t * (B[1] - A[1])
    return x_t, y_t

def Polygon(poly):
    for count in range(0, len(poly) - 1):
        plt.plot(*Straight(poly[count], poly[count + 1], np.linspace(0, 1, 100)))
    return plt.plot(*Straight(poly[len(poly) - 1], poly[0], np.linspace(0, 1, 100)))

def Hyperbola(lamba, O, A, B, t):
    x_t = O[0] + lamba * (np.exp(t) * (A[0] - O[0]) + np.exp(-t) * (B[0] - O[0]))
    y_t = O[1] + lamba * (np.exp(t) * (A[1] - O[1]) + np.exp(-t) * (B[1] - O[1]))
    return [x_t, y_t]
   
def DiscrNumber(sigma, t1, t2, eps):
    return round(math.fabs((sigma/(12 * eps))**(1/3) * (t2 - t1)))
    
def BuildPoly(poly, lamba, I, sk, mi, t):
    P0 = Hyperbola(lamba, I, poly[sk], poly[mi], t)[0]
    P1 = Hyperbola(lamba, I, poly[sk], poly[mi], t)[1]
    P = np.array([P0, P1])
    return P

def Erosion(sigma, eps, iterat):
    #Print start polygon
    Polygon(list_poly)
    
    #Create array of iter+1
    poly_cn = list_poly
    
    #Create array of iter
    poly_c = []
    
    for it in range(0, iterat):
        
        #array of iter+1(poly_cn) write in array poly_c
        poly_c.clear()
        poly_c.extend(poly_cn)
        poly_cn.clear()   
        numb = len(poly_c)
        
        for i in range(0, numb):
            for k in range(i, numb + i):
                
                #Find true indexes
                [m_i, ii, kk, s_k] = Indexes(poly_c, i, k)
                #Point intersect
                I = [Intersect(poly_c[ii], poly_c[m_i], poly_c[kk], poly_c[s_k])[0], Intersect(poly_c[ii], poly_c[m_i], poly_c[kk], poly_c[s_k])[1]]
                #if P_i,P_i-1 and P_kP_k+1 not intersect and P_i != P_k then k++
                if ((Intersect(poly_c[ii], poly_c[m_i], poly_c[kk], poly_c[s_k])[2] == 0) and (ii != kk)):
                    continue
                else:
                    #Area of Sigma_i_k
                    ik_sigma = AreaSigma(poly_c, ii, kk) 
                    ssigma = sigma + ik_sigma 
  
                    if ((AreaTriangle(I, poly_c[kk], poly_c[ii]) <= ssigma) and (AreaTriangle(I, poly_c[s_k], poly_c[m_i]) >= ssigma)):
                        lamba = math.sqrt(ssigma/(4 * (AreaTriangle(I, poly_c[s_k], poly_c[m_i]))))
                        # Hyperbola M(t) t1<=t<=t2
                        [t1, t2] = EndPts(poly_c, ii, kk, sigma, lamba)
                        #Build Hyperbola
                        plt.plot(*Hyperbola(lamba, I, poly_c[s_k], poly_c[m_i], np.linspace(t1, t2, 100)))

                        #Point of next polygon
                        numbnext = DiscrNumber(ssigma, t1, t2, eps)
                        #Find point next polygon and push in new array 
                        for s in range(0, numbnext):
                            ts = (1 - s/numbnext) * t1 + (s/numbnext) * t2
                            poly_cn += [BuildPoly(poly_c, lamba, I, s_k, m_i, ts)]

    plt.show()

#Print function, which build sigma-erosion with sigma=0.5 or sigma=..., eps=0.01 and use iteration=5
#Please change one of print(Erosion)
print(Erosion(0.5, 0.01, 5))
#print(Erosion(0.1, 0.01, 5)) 
#print(Erosion(1, 0.01, 5))
#print(Erosion(1.5, 0.01, 5))  