pro tsfmt

nmax=14L
tsary=findgen(nmax)*100.
print, tsary
print, tsary, format='("nmf2/fac_ne=[",5(",",f12.7),"$")'
print, "]"

print, 'tsfmt finished'
end ; pro tsFmt
