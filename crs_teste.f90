PROGRAM crs

	IMPLICIT REAL*8(a-h,o-z)
	REAL(8) conj_sol(3,21), norma(21)
	REAL(8) xqmin,xqmax,zqmin,zqma,mmin,mmax, Rrms
	REAL(8) P(3), finalmod(3)
	REAL*8 perfil(10000),anom(10000),anomc(10000)
	REAL*8 mi,md,mom
	INTEGER(4) modelos, iteracoes
	
	OPEN(2,file='anomalia.txt') 
	OPEN(3,file='grid.txt') 
	OPEN(4,file='phi.txt')


	ig=1
	DO WHILE (.true.)
		READ(2,*,end=111) a1,a2
		perfil(ig)=a1
		anom(ig)=a2
		ig=ig+1
	END DO
	111	CONTINUE
	CLOSE(2)

	np=ig-1

	zp=-1d0
	yq=0d0


!!!Inicio do CRS

	xqmin = 0d0
	xqmax = 100d0
	zqmin = 0d0
	zqmax = 100d0
!	mmin = 1400d0
!	mmax = 1900d0

!!! o valor de modelos iniciais geralmente eh 7* n de parametros
!! pode ser diferente, vamos testar 7 primeiro
	modelos = 21

!!!! calculos dos modelos iniciais

	DO i=1,modelos
    	conj_sol(1,i) = rand()*(xqmax-xqmin) + xqmin
    	conj_sol(2,i) = rand()*(zqmax-zqmin) + zqmin
    	conj_sol(3,i) = 1943.3354689399671
	END DO

!!! Começo da etapa iterativa
!!! Obtemos n=3 modelos do conjunto aleatoriamente
!!! E fazemos uma média deles

	S = rand()*modelos
	T = rand()*modelos
	U = rand()*modelos
	k = IDNINT(S)
	l = IDNINT(T)
	m = IDNINT(U)

	print*, k, l, m


!!! construindo a média de n modelos escolhidos aleatoriamente
	DO i = 1,3
    	P(i) = (conj_sol(i,k) + conj_sol(i,l) + conj_sol(i,m))/3
	END DO

	V = rand()*modelos
	n = IDNINT(V)

!!! média de P + o modelo n+1

	DO i = 1,3
    	finalmod(i) = 2*P(i) - conj_sol(i,n)
	END DO


!!  caracteristicas do campo geomagnético local
	mi=-34d0
	md=0d0
	F=23500d0

	xq=finalmod(1)
	zq=finalmod(2)
	mom=finalmod(3)


!    calculo das componentes do campo regional
	
	CALL dircos(mi,md,amx,amy,amz)

	fx=F*amx
	fy=F*amy
	fz=F*amz


!	 ----cálculo da amomalia de campo total através das subrotinas dipole e dircos

	DO i=1,np
		xp=perfil(i)

		CALL dipole(xq,yq,zq,mi,md,mom,xp,yp,zp,bx,by,bz)
		bxx=bx+fx
		byy=by+fy
		bzz=bz+fz
		bt=sqrt(bxx**2+byy**2+bzz**2)

		anomc(i)=bt-F

	END DO

!  cálculo do Rrms()norma
!  para o modelo T será rms, todos os outros serão no vetor norma

	soma=0d0
	DO i=1,np
		soma=soma+(anomc(i)-anom(i))**2
	END DO

	Rrms=dsqrt(soma/dfloat(np))

	DO j=1,modelos
    	xq=conj_sol(1,j)
    	zq=conj_sol(2,j)
    	mom=conj_sol(3,j)

    	DO i=1,np
	    	xp=perfil(i)

	    	CALL dipole(xq,yq,zq,mi,md,mom,xp,yp,zp,bx,by,bz)
	    	bxx=bx+fx
	    	byy=by+fy
	    	bzz=bz+fz
	    	bt=sqrt(bxx**2+byy**2+bzz**2)

	    	anomc(i)=bt-F

    	END DO

    	soma=0d0
    	DO i=1,np
	    	soma=soma+(anomc(i)-anom(i))**2
    	END DO

    	norma(j)=dsqrt(soma/dfloat(np))

	END DO

!!! Fazendo o grid

	DO i=0,100
		DO j = 0,100
			xq = xqmin + i
			zq = zqmin + j
			mom = conj_sol(3,1)
			DO k = 1,np
				xp=perfil(k)

				CALL dipole(xq,yq,zq,mi,md,mom,xp,yp,zp,bx,by,bz)
				bxx=bx+fx
				byy=by+fy
				bzz=bz+fz
				bt=sqrt(bxx**2+byy**2+bzz**2)
		
				anomc(k)=bt-F
		
			END DO
		
			soma=0d0
			DO k=1,np
				soma=soma+(anomc(k)-anom(k))**2
			END DO
		
			WRITE(3,*) xq, zq, dsqrt(soma/dfloat(np))

		END DO
	END DO

	iteracoes = 5000
	cont = 0
!!!!!!!!!!!verificar se a maior norma do vetor normal é menor que 
!!!!!!! Inicio das iteracoes
	DO WHILE (cont.NE.iteracoes)
		cont = cont + 1

    	vmax = MAXVAL(norma)
    	DO i=1,modelos
			IF (norma(i).EQ.vmax) THEN
				imax = i
			END IF
    	END DO

		IF (vmax.LT.Rrms) THEN
			S = rand()*modelos
        	T = rand()*modelos
        	U = rand()*modelos
        	k = IDNINT(S)
        	l = IDNINT(T)
        	m = IDNINT(U)

!!! construindo a média de n modelos escolhidos aleatoriamente
        	DO i = 1,3
            	P(i) = (conj_sol(i,k) + conj_sol(i,l) + conj_sol(i,m))/3
        	END DO

			V = rand()*modelos
        	n = IDNINT(V)
        
!!! média de P + o modelo n+1

        	DO i = 1,3
            	finalmod(i) = 2*P(i) - conj_sol(i,n)
        	END DO
        

        	xq=finalmod(1)
	    	zq=finalmod(2)
	    	mom=finalmod(3)

			DO i=1,np
				xp=perfil(i)
				CALL dipole(xq,yq,zq,mi,md,mom,xp,yp,zp,bx,by,bz)
	        	bxx=bx+fx
	        	byy=by+fy
	        	bzz=bz+fz
	        	bt=sqrt(bxx**2+byy**2+bzz**2)

	        	anomc(i)=bt-F

	    	END DO

!  cálculo do Rrms()norma
!  para o modelo T será rms, todos os outros serão no vetor norma

	    	soma=0d0
	    	DO i=1,np
	        	soma=soma+(anomc(i)-anom(i))**2
	    	END DO

	    	Rrms=dsqrt(soma/dfloat(np))

        	DO j=1,modelos
            	xq=conj_sol(1,j)
            	zq=conj_sol(2,j)
            	mom=conj_sol(3,j)

	        	DO i=1,np
	            	xp=perfil(i)

	            	CALL dipole(xq,yq,zq,mi,md,mom,xp,yp,zp,bx,by,bz)
	            	bxx=bx+fx
	            	byy=by+fy
	            	bzz=bz+fz
	            	bt=sqrt(bxx**2+byy**2+bzz**2)

	            	anomc(i)=bt-F

	        	END DO

	        	soma=0d0
	        	DO i=1,np
	        	soma=soma+(anomc(i)-anom(i))**2
	        	END DO

	        	norma(j)=dsqrt(soma/dfloat(np))

        	END DO

		ELSE
			conj_sol(1,imax)=finalmod(1)
			conj_sol(2,imax)=finalmod(2)
			conj_sol(3,imax)=finalmod(3)

			S = rand()*modelos
        	T = rand()*modelos
        	U = rand()*modelos
			k = IDNINT(S)
        	l = IDNINT(T)
        	m = IDNINT(U)

!!! construindo a média de n modelos escolhidos aleatoriamente
        	DO i = 1,3
            	P(i) = (conj_sol(i,k) + conj_sol(i,l) + conj_sol(i,m))/3
        	END DO

			V = rand()*modelos
        	n = IDNINT(V)
        
!!! média de P + o modelo n+1

        	DO i = 1,3
            	finalmod(i) = 2*P(i) - conj_sol(i,n)
        	END DO
        

			xq=finalmod(1)
			xq=finalmod(2)
	    	mom=finalmod(3)

			DO i=1,np
				xp=perfil(i)

				CALL dipole(xq,yq,zq,mi,md,mom,xp,yp,zp,bx,by,bz)
				bxx=bx+fx
				byy=by+fy
				bzz=bz+fz
				bt=sqrt(bxx**2+byy**2+bzz**2)
			
				anomc(i)=bt-F

	    	END DO

!  cálculo do Rrms()norma
!  para o modelo T será rms, todos os outros serão no vetor norma

	    	soma=0d0
	    	DO i=1,np
	        	soma=soma+(anomc(i)-anom(i))**2
	    	END DO

	    	Rrms=dsqrt(soma/dfloat(np))

			DO j=1,modelos
				xq=conj_sol(1,j)
				zq=conj_sol(2,j)
				mom=conj_sol(3,j)
			
				DO i=1,np
					xp=perfil(i)

					CALL dipole(xq,yq,zq,mi,md,mom,xp,yp,zp,bx,by,bz)
					bxx=bx+fx
					byy=by+fy
					bzz=bz+fz
					bt=sqrt(bxx**2+byy**2+bzz**2)

					anomc(i)=bt-F

				END DO

				soma=0d0
				DO i=1,np
					soma=soma+(anomc(i)-anom(i))**2
				END DO

		    	norma(j)=dsqrt(soma/dfloat(np))

			END DO

		END IF

		vmin = MINVAL(norma)
    	DO i=1,modelos
			IF (norma(i).EQ.vmin) THEN
				imin = i
			END IF
    	END DO
		WRITE(4,*) conj_sol(1,imin), conj_sol(2,imin), vmin
	
	END DO !!!! fim das iteracoes


	vmin = MINVAL(norma)
    	DO i=1,modelos
			IF (norma(i).EQ.vmin) THEN
				imin = i
			END IF
    	END DO

	xq=conj_sol(1,imin)
	zq=conj_sol(2,imin)
	mom=conj_sol(3,imin)
	
	DO i = 1,np

		xp=perfil(i)

		CALL dipole(xq,yq,zq,mi,md,mom,xp,yp,zp,bx,by,bz)
		bxx=bx+fx
		byy=by+fy
		bzz=bz+fz
		bt=sqrt(bxx**2+byy**2+bzz**2)

		anomc(i)=bt-F

	!	WRITE(3,*) perfil(i),anomc(i),anom(i)
		WRITE(6,*) perfil(i),anomc(i),anom(i)
	END DO
	
	WRITE(6,*) 'solucao:' 
	WRITE(6,*) 'xq=', conj_sol(1,imin)
	WRITE(6,*) 'zq=', conj_sol(2,imin)
	WRITE(6,*) 'mom=', conj_sol(3,imin)
	WRITE(6,*) 'rms', norma(imin) 

END PROGRAM

!!!!!!! SUBROTINAS  !!!!!!!!!!

	subroutine dipole(xq,yq,zq,mi,md,moment,xp,yp,zp,bx,by,bz)

!     VARIAVEIS UTILIZADAS
!	xp,yp e zp = COORDENADAS DO PONTO DE OBSERVAÇÃO (metros)
!	xq,yq,zq = coordenadas do centro da esfera (metros)
!	mi = inclinação da magnetização (graus, positivo abaixo da horizontal) 
!	md = declinação da magnetização (graus, positiva para Leste, a partir do Norte verdadeiro)
!	moment = momento de dipolo magnetico (A.m²)
!	bx,by,bz,t = elementos do campo B gerados delo dipolo - dado de saida (nT)
!	eixo x, na direção Norte
!	eixu y, na direção Leste
!	eixo z, para baixo 

	implicit real*8 (a-h,o-z)
    real*8 mi,md,mx,my,mz,moment
	data t2nt/1.d9/,cm/1.d-7/
	pi=2.d0*acos(0.d0)
	call dircos(mi,md,mx,my,mz)
	rx=xp-xq
	ry=yp-yq
	rz=zp-zq
	r2=rx**2+ry**2+rz**2
	r=sqrt(r2)
	if(r.eq.0.)pause 'DIPOLE: bad argument detected.'
	r5=r**5
	dot=rx*mx+ry*my+rz*mz
	bx=cm*moment*(3.d0*dot*rx-r2*mx)/r5
	by=cm*moment*(3.d0*dot*ry-r2*my)/r5
	bz=cm*moment*(3.d0*dot*rz-r2*mz)/r5
	bx=bx*t2nt
	by=by*t2nt
	bz=bz*t2nt
	return
	end	
!!!!!!!!!!!!!!!!!!

    subroutine dircos(incl,decl,a,b,c)
	implicit real*8 (a-h,o-z) 
    real*8 incl
	pi=2.d0*acos(0.d0)
	d2rad=pi/180.d0
	xincl=incl*d2rad
	xdecl=decl*d2rad
	a=cos(xincl)*cos(xdecl)
	b=cos(xincl)*sin(xdecl)
	c=sin(xincl)
	return
	end		
		

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    


