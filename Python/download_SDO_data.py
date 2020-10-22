"""
  Download SDO data with drms
"""

import drms

# -- Select series
series = 'hmi.M_45s'
#series = 'aia.lev1_euv_12s'

# -- Select data intervals
tsel = '2014.02.04_03:58:00_TAI-2014.02.04_04:02:00_TAI'
#tsel = '2012.07.04_13:03:00_TAI-2012.07.04_13:55:00_TAI'
#tsel = '2012.07.04_15:20:00_TAI-2012.07.04_16:10:00_TAI'

# -- Wavelenght
wave1 = 4500
#wave2 = 171

# -- Make string
qstr = '%s[%s]' % (series, tsel)#, wave2)
#qstr = '%s[%s][%d]' % (series, tsel, wave1)#, wave2)

keys = ['T_REC', 'T_OBS', 'DATAMIN', 'DATAMAX', 'DATAMEAN', 'DATARMS',\
    'DATASKEW', 'DATAKURT', 'QUALITY']
print('Data export query:\n  %s\n' % qstr)

email = 'andmartinezci@unal.edu.co'
c = drms.Client(verbose=True)
r = c.export(qstr, method='url', protocol='fits', email=email)

print('\nRequest URL: %s' % r.request_url)
print('%d file(s) available for download.\n' % len(r.urls))

out_dir='/home/angel/'
r.download(out_dir)
print('Download finished.')

exit()





c = drms.Client()
print(c.series(r'aia.lev1'))

c = drms.Client(email='andmartinezci@unal.edu.co',verbose=True)
r = c.export('aia.lev1[2014-02-04T03:59:00Z/1m][335]',protocol='fits')

print(r.data.filename)
print(r.urls.url[0])
r.download(".")

exit()



c = drms.Client(email='andmartinezci@unal.edu.co',verbose=True)
r = c.export('hmi.v_45s[2014.02.04_03:59:00_TAI-2014.02.04_04:01:00_TAI]{continuum}',protocol='fits')

print(r.data.filename)
print(r.urls.url[0])
r.download(".")#/run/media/angel/KINGSTON/data/")#'/run/media/angel/VERBATIM16G/', None, verbose=True)


