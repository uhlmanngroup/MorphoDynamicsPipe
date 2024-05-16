import numpy as np

######################################################################
########### function cpda #####################################
######################################################################

## The following code is based on an implementation on mathwork by Awrangjeb (Awrangjeb, Robust Image Corner Detection based on the Chord-to-Point Distance Accumulation Technique).

# L is the chord length
def cpda(curve, L, closed=True):
    xs = curve[0,:]
    ys = curve[1,:]
    curveLength = len(xs)  
    if curveLength <= L:
        L = curveLength-1

    H_L = np.zeros((curveLength))
    for k in range(0,curveLength):
        xk = xs[k] # (x1,y1) = point at which distance will be accumulated
        yk = ys[k]      
                        
        for i in range(k-L,k):
            if i<0 and not closed:
                continue
            else:
                x1 = xs[i] # (leftx,lefty) = current left point for which distance will be accumulated
                y1 = ys[i] 

            if i+L<curveLength:
                x2 = xs[i+L] # (rightx,righty) = current right point for which distance will be accumulated
                y2 = ys[i+L]
            else:
                if closed:
                    x2 = xs[i+L-curveLength] # (rightx,righty) = current right point for which distance will be accumulated
                    y2 = ys[i+L-curveLength]
                else:
                    break

            a = y2-y1 # coefficients of st. line through points (x1,y1) and (x2,y2)
            b = x1-x2
            c = x2*y1 - x1*y2
            if a==0 and b == 0:
                dist = 0
            else:
                dist = (a*xk + b*yk + c)/np.sqrt(a*a+b*b)
            H_L[k] = H_L[k] + dist

    return H_L
