from matplotlib import pyplot as plt
import numpy as np
def pl(fname, tit):
    with open(fname, 'r') as tmp:
        tmp = tmp.readlines()
        n = 0
        for line in tmp:
            if line.split()[0] == '#':
                break
            else:
                n += 1
        arr = np.zeros((n, n))
        i = 0
        l = 10
        for  _, line in enumerate(tmp):
            


            #print(i, line.split()[0])
            if line.split()[0] == '#':
                plt.figure(tit)
                plt.title(np.sum(arr))
                plt.imshow(arr.transpose(), origin='lower', extent=(0,1,0,1))
                plt.colorbar()
                plt.savefig(tit +str(l)+ '.png')
                plt.clf()
                #plt.show()
                l += 1

                i = 0
            else:
                arr[i ,:] = list(map(float, line.split()))
                i += 1

def xypl(fname, tit):
    tmp1 = open('xy1.out', 'r') 
    tmp1 = tmp1.readlines()
    tmp2 = open('xy2.out', 'r') 
    tmp2 = tmp2.readlines()

    k1 = 0
    k2 = 0

    t = 0
    print(len(tmp1))
    print(len(tmp2))
    while (k1 < len(tmp1)-1) and (k2 < len(tmp2) - 1):
        plt.figure()
        xs = []
        ys = []
        for i, line in enumerate(tmp1[k1:]):
            if (line.split()[0] == '#'):
                k1 += i+1
                break
            x, y = map(float, line.split())
            xs.append(x)
            ys.append(y)
        plt.plot(xs, ys, 'or', markersize=1, alpha=0.1)
    
        xs = []
        ys = []
        for j, line in enumerate(tmp2[k2:]):
            if (line.split()[0] == '#'):
                k2 += j+1
                break
            x, y = map(float, line.split())
            xs.append(x)
            ys.append(y)

        plt.plot(xs, ys, 'ob', alpha=0.1, markersize=1)
        for i in range(10):
            plt.plot([i*0.1, i*0.1], [0, 1], '--k', alpha=0.1)
            plt.plot([0, 1], [i*0.1, i*0.1], '--k', alpha=0.1)
        plt.xlim(0, 1)
        plt.ylim(0, 1)
        plt.savefig('xy%s.png' %t)
        plt.cla()
        plt.close()
        t += 1
        print(k1, k2)

#pl('N.out', 'N')
#pl('de.out', 'dE')
pl('e0.out', 'E_b')
#pl('e1.out', 'E_u')
#pl('u.out', 'U')
#pl('v.out', 'V')
#pl('m.out', 'M')
#pl('w.out', 'W')
#pl('dw.out', 'dW')
#pl('p.out', 'P')
#xypl('xy2.out', 'xy')
