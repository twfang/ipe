pro get_v20,time_string $
,tec4,nmf24,hmf24 $
,tec5,nmf25,hmf25 $
,tec6,nmf26,hmf26

;136ut
;v21
restore ,'tmp20150730_ipe_80_1453556049_'+time_string+'.sav'
tec4sv=tec0
nmf24sv=nmf20
hmf24sv=hmf2


;v22
restore ,'tmp20150730_ipe_80_1453558480_'+time_string+'.sav'
tec5sv=tec0
nmf25sv=nmf20
hmf25sv=hmf2

;v23
restore ,'tmp20150730_ipe_80_1453558551_'+time_string+'.sav' 
tec6sv=tec0
nmf26sv=nmf20
hmf26sv=hmf2


tec4=tec4sv
tec5=tec5sv
tec6=tec6sv
nmf24=nmf24sv
nmf25=nmf25sv
nmf26=nmf26sv
hmf24=hmf24sv
hmf25=hmf25sv
hmf26=hmf26sv
end


