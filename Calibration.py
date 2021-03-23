from readspc import spectrum2D
import xlrd

sheet = xlrd.open_workbook('spectral_radiance_Wi_17G_73161PTB20.xlsx')
sheet = sheet.sheet_by_index(0)
L = []
for i in range(7,sheet.nrows):
    row = [sheet.cell_value(i, 0), sheet.cell_value(i, 1), sheet.cell_value(i, 2)]
    L.append(row)

frame = spectrum2D('FS_x2_500-710/FS_x2_  32.spc')
dtmeas = frame.gate_time
frame = spectrum2D('FS_x2_500-710/FS_x2_bg_  32.spc')
dtcalib = frame.gate_time
wavelength = frame.wavelength

intensity = []
for i in range(1,10):
    if i < 10:
        num = '0' + str(i)
    else:
        num = str(i)
    frame = spectrum2D('FS_x2_500-710/FS_x2_  33-SG-700-910-Step-1-raw-Frame-'+num+'.spc')
    intensity.append(frame.intensity)
BGintensity = []
for i in range(1,5):
    num = str(i)
    frame = spectrum2D('FS_x2_500-710/FS_x2_bg_  33-SG-700-910-Step-1-raw-Frame-'+num+'.spc')
    BGintensity.append(frame.intensity)

Smeas = []
maxx = len(intensity[0])
maxy = len(intensity[0][0])
nframe = len(intensity)
nBGframe = len(BGintensity)
for x in range(maxx):
    col = []
    for y in range(maxy):
        meanU = 0
        for i in range(nframe):
            meanU = meanU + intensity[i][x][y]
        meanUBG = 0
        for i in range(nBGframe):
            meanUBG = meanUBG + BGintensity[i][x][y]
        col.append(meanU/(nframe)-meanUBG/(nBGframe))
    Smeas.append(col)
    print(x)


Lcalib = []
for wl in wavelength:
    Ldown = L[0][0]
    i = 0
    while Ldown < wl.value:
         i = i + 1
         Ldown = L[i][0]
    Lup = L[i-1][0]
    ratio = (wl.value-Ldown)/(Lup-Ldown)
    Lcalib.append(L[i][1]*(1-ratio) + L[i-1][1]*ratio)
#print(Smeas)
#print(dtmeas)
#print(Scalib)
#print(dtcalib)
#print(Lcalib)
#Lmeas = (Smeas/dtmeas)/(Scalib/dtcalib)*Lcalib
#print(Lmeas)