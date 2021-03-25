import ROOT
import numpy as np
import matplotlib.pyplot as plt

pmtId = "Pmt3"
inFile = ROOT.TFile.Open ( "zooTraces100binsPMT3.root" ," READ ")

hist = inFile.Get ("stTraces")
#histNo = inFile.Get ("stTracesNo")

arr0 = []
arr1 = []

for e in range(0, hist.GetNbinsX()):
    if hist.GetBinContent(e, 1) > 0:
        arr0 = []
        for i in range( 0, hist.GetNbinsY()+1 ):
            arr0.append( hist.GetBinContent(e, i) )
        arr1.append( arr0 )

'''
for i in range(0, len(arr1)):
    plt.plot(arr1[i])

plt.xlim(1,2100)
#plt.xlabel("Bins (from 0 to 99, first bins; from 100 to 200, last bins)")
plt.xlabel("Time / (8.3 ns)")

#plt.ylim(140,270)
plt.ylabel("FADC")

plt.title( "Traces for Station 1851" )
#plt.title( "Traces for Station 1851\n("+str(len(arr1))+"/608)" )
plt.savefig("../../plots/traces1851"+pmtId+".png", dpi=300)
plt.show()
'''
# ======================================
# ======================================

print("Doing first and last bins")

arr2 = []
arr3 = []
tmp = True

for e in range(0, hist.GetNbinsX()):
    if hist.GetBinContent(e, 1) > 0:
        arr2 = []
        tmp = True
        for i in range( 1, 101 ):
            if hist.GetBinContent(e, i) < 235 or hist.GetBinContent(e, i) > 255:
                tmp = False
        for i in range( 1947, 2046 ):
            if hist.GetBinContent(e, i) < 235 or hist.GetBinContent(e, i) > 255:
                tmp = False
        if tmp:
            for i in range( 0, hist.GetNbinsY()+1 ):
                arr2.append( hist.GetBinContent(e, i) )
            arr3.append( arr2 )

nb = 100

xf = np.linspace(0, nb, nb+1)
xl = np.linspace(nb+1, 2*nb+1, nb+1)
tmpfb = []
tmplb = []
fb = []
lb = []

for e in range(0, len(arr1)):
    tmpfb = []
    tmplb = []
    for i in range(1, nb+2):
        tmpfb.append( arr1[e][i] )
        tmplb.append( arr1[e][-i] )
    fb.append(tmpfb)
    lb.append(tmplb)

full = []
tmp = []

for e in range(0, len(arr1)):
    tmp = []
    for i in range(1, 2*nb):
        if i<nb+1:
            tmp.append( fb[e][i] )
        else:
            tmp.append( lb[e][i-nb] )
    full.append(tmp)
full = np.array(full)
for i in range(0, len(arr1)):
    #plt.plot(xf, fb[i])
    #plt.plot(xl, lb[i])
    plt.plot(full[i])

plt.xlim(0,2*nb)
plt.xlabel("Bins (from 0 to 99, first bins; from 100 to 200, last bins)")
plt.xlabel("Time / (8.3 ns)")
plt.axvline(x=nb, c='k', lw=4)
plt.ylabel("FADC")

plt.title( "First and Last 100 bins for station 1851's traces" )
#plt.title( "Traces for Station 1851\n("+str(len(arr1))+"/608)" )
plt.savefig("../../plots/traces1851FL"+pmtId+".pdf", dpi=300)
plt.show()

#nb = 100
#xf = np.linspace(0, nb, nb + 1)
#xl = np.linspace(nb+1, 2*nb+1, nb + 1)
#tmpfb = []
#tmplb = []
#fb = []
#lb = []
#
#for e in range(0, len(arr3)):
#    tmpfb = []
#    tmplb = []
#    for i in range(1, nb+2):
#        tmpfb.append( arr3[e][i] )
#        tmplb.append( arr3[e][-i] )
#    fb.append(tmpfb)
#    lb.append(tmplb)
#
#for i in range(0, len(arr3)):
#    plt.plot(xf, fb[i])
#    plt.plot(xl, lb[i])
#
#
#plt.xlim(0,nb*2)
#plt.xlabel("Bins (from 0 to 99, first bins; from 100 to 200, last bins)")
#plt.xlabel("Time / (8.3 ns)")
#plt.axvline(x=100, c='k', lw=4)
#plt.ylim(235,255)
#plt.ylabel("FADC")
#
#plt.title( "First and Last 100 bins for station 1851's traces" )
##plt.title( "Traces for Station 1851\n("+str(len(arr1))+"/608)" )
#plt.grid(axis='y')
#plt.savefig("../../plots/traces1851FLzoom"+pmtId+".pdf", dpi=300)
#plt.show()


#arr0 = []
#arr1 = []
#
#for e in range(0, histNo.GetNbinsX()):
#    if histNo.GetBinContent(e, 1) > 0:
#        arr0 = []
#        for i in range( 0, histNo.GetNbinsY()+1 ):
#            arr0.append( histNo.GetBinContent(e, i) )
#        arr1.append( arr0 )
#
#
#xf = np.linspace(0, 100, 101)
#xl = np.linspace(101, 201, 101)
#tmpfb = []
#tmplb = []
#fb = []
#lb = []
#
#for e in range(0, len(arr1)):
#    tmpfb = []
#    tmplb = []
#    for i in range( 1, 102):
#        tmpfb.append( arr1[e][i] )
#        tmplb.append( arr1[e][-i] )
#    fb.append(tmpfb)
#    lb.append(tmplb)
#
#
#
#for i in range(0, len(arr1)):
#    plt.plot(xf, fb[i])
#    plt.plot(xl, lb[i])
#
#
#plt.xlim(1,200)
##plt.xlabel("Bins (from 0 to 99, first bins; from 100 to 200, last bins)")
#plt.xlabel("bin / (8.3 ns)")
#
#plt.ylim(200,280)
#plt.ylabel("FADC")
#
#plt.title( "Traces for Station 1851\n("+str(len(arr1))+"/608)" )
#plt.grid("on")
#plt.savefig("../../plots/traces1851NoFirstLast.pdf", dpi=300)
