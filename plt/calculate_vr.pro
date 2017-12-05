pro calculate_vr, vr,ve,mp,lp,VEXBTH,VEXBE, z_km,mlat_deg,jmin_in,jmax_is

  earth_radius = 6.3712E+06     ;!..Earth radius [meter]
  ht90  = 90.0E+03              ;!reference height in meter
  r90 = earth_radius + ht90     ;![m]
  time_step=10.                 ;[sec]


   theta_t1 = -mlat_deg[JMIN_IN[lp]]*!pi/180.0+!pi*0.50 ;plasma_grid_GL NH
   theta_t0 = theta_t1 - ( VEXBth[lp,mp] * time_step ) / r90
   coslambda_m  = COS ( !pi*0.50 - theta_t0 )
   r0_apex = ( earth_radius + ht90 ) *  coslambda_m *  coslambda_m ;[meter]
   midpoint = JMIN_IN[lp] + ( JMAX_IS[lp] - JMIN_IN[lp] )/2
   ;z_km=plasma_grid_Z*1.0E-3
   plasma_grid_Z = z_km[midpoint]*1.0E+3
   r1 = plasma_grid_Z + earth_radius ;[meter]
   Vr = (r1-r0_apex)/time_Step ;meter/second
;Veast
   Ve = VEXBe[lp,mp]*(r1)/(r90)
end;pro calculate_vr
