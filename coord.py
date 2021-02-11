def CalcCoord2(fl, sat, yr, day, wkdy, cor, lps, day_2, month):



	#######################################################################################################
	# This function was developed by Loram Siqueira and João Francisco Galera Monico - 2020 ###############
	# If any extra information is needed please send an e-mail to loram.siqueira@unesp.br #################
    #######################################################################################################


	######################################################################################################	
	# fl is a file read through file.readlines() e must contain the broadcast ephemeris without header####
	

    for j in range(len(fl)):
        sc = 0
        if fl[j].startswith(sat):
		###############################################################
		# This section calculates the coordinates for GLONASS #########
		###############################################################
		
            if 'R' in sat:
                #Constants #########################################
                mu = 398600.44
                c20 = -1082.6257e-6
                ae = 6378.136
                we = 0.7292115e-4
                Tu = jdn(yr,month, day_2)/36525
                H0 = (24110.54841 + 8640184.812866*Tu + 0.093104*(Tu ** 2) - 6.2e-6 *(Tu ** 3))#*(2*pi)/86400
                w3 = (we + 4.3e-15*Tu)*86400/(2*pi)


                ####################################################
                ## Elements coming from the the nav file          ##
                ####################################################
                TauN = float(fl[j][23:42])
                GammaN = float(fl[j][42:61])
                tb = int(float(fl[j][61:80]))
                X0 = float(fl[j + 1][0:23])
                Xdot = float(fl [j+1][23:42])
                xdotdot = float(fl[j + 1][42:61])
                Y0 = float(fl[j + 2][0:23])
                Ydot = float(fl[j + 2][23:42])
                ydotdot = float(fl[j + 2][42:61])
                Z0 = float(fl[j + 3][0:23])
                Zdot = float(fl[j + 3][23:42])
                zdotdot = float(fl[j + 3][42:61])

                if X0 == 0:
                    break
				######################################################
                ## Start of calculations #############################
                ######################################################

                for i in range (tb-15*60, tb+15*60):
                    if (i % 900 == 0):
                        Ste = (H0 + (we*(i - 3*60*60))*86400/(2*pi))%86400
                        Ste = Ste*(2*pi)/86400
                        Xi = X0*math.cos(Ste) - Y0*math.sin(Ste)
                        Yi = X0*math.sin(Ste) + Y0*math.cos(Ste)
                        Zi = Z0
                        Xdoti = Xdot*math.cos(Ste) - Ydot*math.sin(Ste) - we*Y0
                        Ydoti = Xdot*math.sin(Ste) + Ydot*math.cos(Ste) + we*X0
                        Zdoti = Zdot

					####################################################
					##### Initial conditions for integration ###########
					####################################################
                        w10 = X0
                        w20 = Y0
                        w30 = Z0
                        w40 = Xdot
                        w50 = Ydot
                        w60 = Zdot
                        deltat = i - tb -lps
                        h = sign(deltat)

                        for j in range (abs(deltat)):
                            r = math.sqrt((w10**2) + (w20**2) + (w30**2));

                            k11 = h * (w40);
                            k12 = h * (w50);
                            k13 = h * (w60);


                            k14 = h * ((-mu*w10) / (r**3) + (((3.0 / 2.0) * c20 * mu * (ae ** 2)*w10) / (r **5) ) * (1 - 5*(w30**2) / (r ** 2) ) + (we**2) * w10 + 2 * we * w50 + xdotdot)

                            k15 = h * (-mu*w20 / (r**3) + (((3.0/2.0) * c20 * mu * (ae**2)*w20) / (r**5)) * (1 - (w30**2)*5 / (r**2) ) + (we**2) * w20 - 2 * we * w40 + ydotdot)

                            k16 = h * (-mu*w30 / (r**3) + (((3 / 2) * c20 * mu * (ae**2)*w30) / (r**5) ) * (3 - 5*(w30**2) / (r**2) ) + zdotdot)

                            
                            k21 = h * (w40 + 0.5 * k14);
                            k22 = h * (w50 + 0.5 * k15);
                            k23 = h * (w60 + 0.5 * k16);

                            k24 = h * ((-mu * (w10 + 0.5 * k11))/ (r**3) + (3 / 2 * c20 * mu * (ae**2) * (w10 + 0.5 * k11) )/ (r**5) *(1 - 5 * (w30 + 0.5 * k13)**2/ (r**2) ) + (we**2) * (w10 + 0.5 * k11) + 2 * we * (w50 + 0.5 * k15) + xdotdot)
                            k25 = h * (-mu / (r**3) * (w20 + 0.5 * k12) + 3 / 2 * c20 * mu * (ae**2) / (r**5) * (w20 + 0.5 * k12) *(1 - 5 / (r**2) * (w30 + 0.5 * k13)**2) + (we**2) * (w20 + 0.5 * k12) - 2 * we * (w40 + 0.5 * k14) + ydotdot)

                            k26 = h * (-mu / (r**3) * (w30 + 0.5 * k13) + 3 / 2 * c20 * mu * (ae**2) / (r**5) * (w30 + 0.5 * k13) *(3 - 5 / (r**2) * ((w30 + 0.5 * k13)**2)) + zdotdot)
                            k31 = h * (w40 + 0.5 * k24);
                            k32 = h * (w50 + 0.5 * k25);
                            k33 = h * (w60 + 0.5 * k26);

                            k34 = h * (-mu*(w10 + 0.5 * k21) / (r**3)  + 3.0 / 2.0 * c20 * mu * (ae**2)* (w10 + 0.5 * k21) / (r**5) *(1 - 5 * (w30 + 0.5 * k23)**2)/ (r**2)  + (we**2) * (w10 + 0.5 * k21) + 2 * we * (w50 + 0.5 * k25) + xdotdot)
                            k35 = h * (-mu / (r**3) * (w20 + 0.5 * k22) + 3 / 2 * c20 * mu * (ae**2) / (r**5) * (w20 + 0.5 * k22)*(1 - 5 / (r**2) * (w30 + 0.5 * k23)**2) + (we**2) * (w20 + 0.5 * k22)-2 * we * (w40 + 0.5 * k24) + ydotdot);
                            k36 = h * (-mu / (r**3) * (w30 + 0.5 * k23) +3 / 2 * c20 * mu * (ae**2) / (r**5) * (w30 + 0.5 * k23) *(3 - 5 / (r**2) * ((w30 + 0.5 * k23)**2)) + zdotdot);

                            k41 = h * (w40 + k34);
                            k42 = h * (w50 + k35);
                            k43 = h * (w60 + k36);

                            k44 = h * (-mu / (r**3) * (w10 + k31) + 3 / 2 * c20 * mu * (ae**2) / (r**5) * (w10 + k31) * (1 - 5 / (r**2) * ((w30 + k33)**2)) + (we**2) * (w10 + k31) + 2 * we * (w50 + k35) + xdotdot)
                            k45 = h * (-mu / (r**3) * (w20 + k32) + 3 / 2 * c20 * mu * (ae**2) / (r**5) * (w20 + k32)*(1 - 5 / (r**2) * ((w30 + k33)**2)) + (we**2) * (w20 + k32) - 2 * we * (w40 + k34) + ydotdot)
                            k46 = h * (-mu / (r**3) * (w30 + k33) + 3 / 2 * c20 * mu * (ae**2) / (r**5) * (w30 + k33) * (3 - 5 / (r**2) * ((w30 + k33)**2)) + zdotdot);


                            w11 = w10 +  (k11 + 2 * k21 + 2 * k31 + k41)/6;
                            w21 = w20 +  (k12 + 2 * k22 + 2 * k32 + k42)/6;
                            w31 = w30 +  (k13 + 2 * k23 + 2 * k33 + k43)/6;
                            w41 = w40 +  (k14 + 2 * k24 + 2 * k34 + k44)/6;
                            w51 = w50 +  (k15 + 2 * k25 + 2 * k35 + k45)/6;
                            w61 = w60 +  (k16 + 2 * k26 + 2 * k36 + k46)/6;

                            w10 = w11;
                            w20 = w21;
                            w30 = w31;
                            w40 = w41;
                            w50 = w51;
                            w60 = w61;
						############################################################
						#### x, y and z position and velocity ######################
						### de position is referred to WGS84 #######################
						############################################################
						
                        vx = '{:14.6f}'.format(w40*1000)
                        vy = '{:14.6f}'.format(w50*1000)
                        vz = '{:14.6f}'.format(w60*1000)
                        x = '{:16.6f}'.format(w10*1000-0.36)
                        y = '{:16.6f}'.format(w20*1000-0.08)
                        z = '{:16.6f}'.format(w30*1000-0.18)
                        cor.write('P' + sat + ' , ' + str(i) + ',' + str(x) + ',' + str(y) + ',' + str(z) + ',' + str(vx) + ',' + str(vy) + ',' + str(vz) + '\n')

			################################################################################
		    ##### This section calculates the coordinates for GPS, Galileo and BDS #########
		    ################################################################################
            if sat[0] in ['E', 'G', 'C']:
			################################################################################
			##### information from the navigation file #####################################
			################################################################################
                a0 = float(fl[j][23:42])
                a1 = float(fl[j][42:61])
                a2 = float(fl[j][61:80])
                iode = float(fl[j + 1][4:23])
                crs = float(fl[j + 1][23:42])
                deltan = float(fl[j + 1][42:61])
                mm0 = float(fl[j + 1][61:80])
                cuc = float(fl[j + 2][4:23])
                e = float(fl[j + 2][23:42])
                cus = float(fl[j + 2][42:61])
                a = (float(fl[j + 2][61:80]))**2
                toe = float(fl[j + 3][4:23])
                cic = float(fl[j + 3][23:42])
                omega0 = float(fl[j + 3][42:61])
                cis = float(fl[j + 3][61:80])
                i0 = float(fl[j + 4][4:23])
                trans = float(fl[j + 7][4:23])
                crc = float(fl[j + 4][23:42])
                w = float(fl[j + 4][42:61])
                omegap = float(fl[j + 4][61:80])
                i = float(fl[j + 5][4:23])
                if 'C' in sat:
                    sc = 14  ### correction for the BDS time reference 
                if 'G' in sat:
                    start = int(toe - 3600)
                    end = int(toe + 3600)
                if 'E' in sat:
                    start = int(toe - 10*60)
                    end = int(toe + 180*60)
                if 'C' in sat:
                    start = int(toe - 60*60)
                    end = int(toe + 60*60)
                for j in range(start, end):
                    count = 0
                    if (j % 900 == 0):
                        t = j
                        if 'G' in sat:
                            GM = 3.986005e14
                            we = 7.2921151467e-5
                        if 'C' in sat:
                            GM = 3.986004418e14
                            we = 7.2921150e-5
                        if 'E' in sat:
                            GM = 3.986004418e14
                            we = 7.2921151467e-5
                        else:
                            GM = 3.986004418e14
                        n0 = math.sqrt(GM / (a ** 3))
                        delta_t = t - toe -sc
                        n = n0 + deltan
                        mk = mm0 + n * delta_t
                        mkdot = n
                        aux = mk
                        for k in range(4):
                            Ek = mk + (e * math.sin(aux))
                            aux = Ek

                        Ek = aux
                        Ekdot = mkdot/(1 - e*math.cos(Ek))

                        dts = ((a0 + (float(j) - float(toe)) * a1 + a2 * (float(j) - float(toe)-14) ** 2))
                        cos_vk = (math.cos(Ek) - e) / (1 - (e * math.cos(Ek)))
                        sen_vk = (math.sqrt(1 - e * e) * math.sin(Ek)) / (1 - (e * math.cos(Ek)))
                        tan_vk = sen_vk / cos_vk

                        vk = arctan2(sen_vk, cos_vk)
                        vkdot = (math.sin(Ek)*Ekdot*(1+e*math.cos(vk)))/(math.sin(vk)*(1 - e*math.cos(Ek)))

                        phik = vk + w
                        deltauk = cuc * math.cos(2 * phik) + cus * math.sin(2 * phik)
                        uk = phik + deltauk

                        deltark = crc * math.cos(2 * phik) + crs * math.sin(2 * phik)
                        rk = a * (1 - e * math.cos(Ek)) + deltark

                        deltaik = cic * math.cos(2 * phik) + cis * math.sin(2 * phik)
                        ik = i0 + i * delta_t + deltaik

                        ukdot = vkdot + 2*(cus*math.cos((2*uk))-cuc*math.sin((2*uk)))*vkdot
                        rkdot = (a*e*math.sin(Ek)*n)/(1-e*math.cos(Ek)) + 2*(crs*math.cos(2*uk)-crc*math.sin(2*uk))*vkdot
                        ikdot = i + (cis*math.cos(2*uk)-cic*math.sin(2*uk))*2*vkdot

                        xp = rk * math.cos(uk)
                        yp = rk * math.sin(uk)

                        omegakdot = (omegap - we)


                        xpkdot = rkdot*math.cos(uk)-yp*ukdot
                        ypkdot = rkdot*math.sin(uk) + xp*ukdot

                        omegak = omega0 + (omegap-we) * delta_t - we * toe

                        x = (xp * math.cos(omegak) - yp * math.cos(ik) * math.sin(omegak))#/1000
                        y = (xp * math.sin(omegak) + yp * math.cos(ik) * math.cos(omegak))#/ 1000
                        z = yp * math.sin(ik)#/1000


                        xkdot = (xpkdot - yp*math.cos(ik)*omegakdot)*math.cos(omegak)- \
                                (xp*omegakdot + ypkdot*math.cos(ik) - yp*math.sin(ik)*ikdot)*math.sin(omegak)

                        ykdot = (xpkdot-yp*math.cos(ik)*omegakdot)*math.sin(omegak)+ \
                                (xp*omegakdot+ypkdot*math.cos(ik) - yp*math.sin(ik)*ikdot)*math.cos(omegak)

                        zkdot = ypkdot*math.sin(ik) + yp*math.cos(ik)*ikdot

					##############################################################
					## Final velocities and position #############################
					##############################################################
                        vx = '{:14.6f}'.format(xkdot)
                        vy = '{:14.6f}'.format(ykdot)
                        vz = '{:14.6f}'.format(zkdot)
                        x ='{:16.6f}'.format(x)
                        y = '{:16.6f}'.format( y )
                        z = '{:16.6f}'.format( z)


                        cor.write('P'+sat + ' , ' + str(j) + ',' + str(x) + ',' + str(y) +',' + str(z) + ',' + str(vx) + ',' + str(vy) +',' + str(vz) +'\n') 

def comparason (final, bdr, cor):
	import math
	from numpy import *
	import datetime
	import lib as lb
	from numpy.linalg import inv

	c_or = open(cor, 'r')
	c_or = c_or.readlines()
	
	year = 2019
	
	#################################################################################
	#### range here is for the days and goes from 1 to 366 depending on the year ####
	#################################################################################
	for i in range (1,2):

    ## Calculating yaw angle##
    data = str(year) +' '+ str(i)
    x = datetime.datetime.strptime(data, '%Y %j')
    day = float(x.strftime('%d'))
    month = float(x.strftime('%m'))
    T = lb.jdn(year,month,day)/36525
    lambdaSun = 280.4606184 + 3600.77005361*T
    Msun = (357.277233+35999.05034*T)*pi/180
    lambdaEcl = (lambdaSun + 1.914666471*math.sin(Msun) + 0.019994643*math.sin(2*Msun))*pi/180
    Obliq = (23.439291 - 0.0130042*T)*pi/180
    r0 = (149.619 - 2.499*math.cos(Msun) - 0.021*math.cos(2*Msun))*10**6*(1000)
    rsun = [math.cos(lambdaEcl), math.cos(Obliq)*math.sin(lambdaEcl), math.sin(Obliq)*math.sin(lambdaEcl)]
    r0sun = ([r0*math.cos(lambdaEcl), r0*math.cos(Obliq)*math.sin(lambdaEcl), r0*math.sin(Obliq)*math.sin(lambdaEcl)])

	###################################################################
	### creating the results files with and without PCO corrections ###
	###################################################################
	
    d = str(i)
    day = d.zfill(3)
    res = open (day, 'w')
    res_cor = open (day + '_cor', 'w')
	
    ############################################################
    ### Abrir arquivo de efemerides precisar e coletar #########
    ### todos os satelites disponiveis                 #########
    ############################################################


    fl = open (final, 'r')
    fl = fl.readlines()
    gl = open (brdc, 'r')
    gl = gl.readlines()
    sat = []
	
	################################################################
	### Creates a list with all available satellites in the file ###
	################################################################
	
    for lines in fl:
        sat.append(lines[11:15])
    sat = list(sorted(set((sat))))


    for sats in sat:
        for i in range(len(c_or)):
            if sats[1:4] in c_or[i]:
                aux = c_or[i + 1].split(',')
                x_cor = float(aux[0])/1000
                y_cor = float(aux[1])/1000
                z_cor = float(aux[2])/1000
        for lines in fl:
            if sats in lines:
                TimeS = int(lines[1:11])
                line1 = lines
                for l in gl:
                    l = l.split(",")
                    if sats in l[0]:
                        TimeB = int(l[1])
                        if 'R' in sats:
                            TimeB = TimeB + 15*60
                        if TimeS == (TimeB):
                            s = (l[0])
                            t = '{:10.2f}'.format(float(l[1]))
                            x = float(l[2])
                            y = float(l[3])
                            z = float(l[4])
                            vx = float(l[5])
                            vy = float(l[6])
                            vz = float(l[7])
                            xp = float(line1[15:29])*1000
                            yp = float(line1[29:43])*1000
                            zp = float(line1[44:56])*1000
                            V3 = [vx, vy, vz]
                            V1 = [x, y, z]
                            V1m = linalg.norm(V1)
                            V2 = [xp, yp, zp]
                            V2m = linalg.norm(V2)

                            ru = V2 / V2m # #Radial component
                            nu = cross(V2, V3) / linalg.norm(cross(V2, V3)) # Along track
                            tu = cross(nu, ru)  # Cross-track

                            V2 = [[xp], [yp], [zp]] #Precise coordinates
                            V1 = [[x], [y], [z]] #broadcast coordinate
                            Rct = matrix([ru,tu,nu]) #rotation matrix for RTN

                            Diff = Rct*V2 - Rct*V1  # Results without correction
                            dradial = float(Diff[0])
                            dalong = float(Diff[1])
                            dcross = float(Diff[2])
                            mod = '{:15.2f}'.format(math.sqrt(dradial ** 2 + dalong ** 2 + dcross ** 2))

                            dradial = '{:15.2f}'.format(dradial)
                            dalong = '{:15.2f}'.format(dalong)
                            dcross = '{:15.2f}'.format(dcross)


                            n = cross(nu, rsun)
                            rsunproj = cross(n,nu)
                            rsat = matrix([xp, yp, zp])
                            k = matrix(-ru)

                            aux = r0sun- transpose(V2)
                            e = (aux / linalg.norm(aux))
                            j = cross(k,e)
                            i = cross(j,k)

                            R= matrix(concatenate((i, j, k )))
                            R = linalg.inv(R)
                            delta = matrix([[x_cor],[y_cor],[z_cor]])
                            rsat2 = transpose(rsat)+R*delta
                            V4 = rsat2

                            Diff2 = Rct*V4 - Rct*V1 #Results with correction
                            dradial2 = float(Diff2[0])
                            dalong2 = float(Diff2[1])
                            dcross2 = float(Diff2[2])
                            mod2 = '{:15.2f}'.format(math.sqrt(dradial2 ** 2 + dalong2 ** 2 + dcross2 ** 2))
                            dradial2 = '{:15.2f}'.format(dradial2)
                            dalong2 = '{:15.2f}'.format(dalong2)
                            dcross2 = '{:15.2f}'.format(dcross2)
                            aux =  abs(float(mod))
                            #print dradial2

                            res.write (t + ','+ s + ',' +  str(dradial) +','+  str(dalong) +','+ str(dcross) +',' + str(mod)+ '\n')
                            res_cor.write(t + ','+ s + ',' +  str(dradial2) +','+  str(dalong2) +','+ str(dcross2) +','+ str(mod2)+ '\n')