!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!			ELVIS FRANK DOMINGUEZ VIDAL						!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!					contenedor de metodos					!
module subrutinas


	contains
	Subroutine normal

		implicit none
		real,dimension(100,100)::S
		integer:: m,n,i,j,respt
		real:: sx,sy,sxy,sx2,a,b,min,max,olgura
		open(unit=1,file='data1.f90')
		print*, ''
		print*, '|----------------------------------------------|'
		print*,'ingreso de la matriz a(i,j)'
		read(1,*)m,n
		do i=1,m
			read(1,*)S(i,1:n)
			print*,S(i,1:n)
		end do
		close(1)
		sx=0
		sy=0
		sxy=0
		sx2=0
		do i=1,m
			sx=sx + S(i,1)
			sy=sy + S(i,2)
			sxy=sxy + S(i,1)*S(i,2)
			sx2=sx2 + S(i,1)**2
		end do
		a=(m*sxy - sx*sy)/(m*sx2 - sx**2)
		b=(sx2*sy - sxy*sx)/(m*sx2 - sx**2)
		print*, ''
		print*, '|----------------------------------------------|'
		print*,' sx=(',sx,')'
		print*,' sy=(',sy,')'
		print*,' sxy=(',sxy,')'
		print*,' sx2=(',sx2,')'
		print*, ''
		print*, '|----------------------------------------------|'
		print*,'el valor de : a=(',a,')'
		print*,'(F12.3)',a
		print*,'el valor de : b=(',b,')'
		print*,'(F12.3)',b
		print*, ''
		print*, '|----------------------------------------------|'
		print*,'La recta es: y=(',a,')','X + (',b,')'
		
		open(unit=1,file='salida.wxm')
		write(1,200)'/* [wxMaxima batch file version 1] [ DO NOT EDIT BY HAND! ]*/'
		write(1,200)'/* [ Created with wxMaxima version 16.04.2 ] */'
		write(1,200)
		write(1,200)'/* [wxMaxima: input   start ] */'
		write(1,*) 'wxdraw2d(color=blue,'
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		do i=1,m
			if(i==1)then
				min=S(i,1)
				max=S(i,1)
			else
				if(max<S(i,1))then
					max=S(i,1)
				endif
				if(min>S(i,1))then
					min=S(i,1)
				endif
			endif			
		end do
		olgura=max-min
		max=max+abs(0.1*olgura)
		min=min-abs(0.1*olgura)
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		write(1,*) 'explicit(',b,'+',a,'*x ,x,',min,',',max,'),'

		write(1,*) 'color=orange-red,'
		write(1,*) 'point_type=5,'
		write(1,*) 'point_size=1,'

		write(1,*) 'points(['
		do i=1,m
			write(1,*) '['
			do j=1,n
				if(j==n)then
					if(i==m)then
						write(1,*),S(i,j),']'
					else
						write(1,*),S(i,j),'],'
					endif
					else
					write(1,*),S(i,j),','
				endif				
			enddo
		end do
		write(1,*) ']),'

		write(1,*) 'xlabel="X",'
		write(1,*) 'ylabel="Y",'
		write(1,*) 'title="Aproximacion");'
		write(1,200)'/* [wxMaxima: input   end   ] */'
		write(1,200)
		write(1,200)'/* Maxima cant load/batch files which end with a comment! */'
		write(1,200)'Created with wxMaxima"$'
		close(1)

		200 FORMAT(T1,A)
		
	End Subroutine normal

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	Subroutine grado(opc)
		implicit none
		integer,intent(inout)::opc
		real,dimension(100,100)::S
		real,dimension(10,10)::mat
		real,dimension(6)::suba
		integer:: m,n,i,j,k,expo,x
		real:: cte,olgura,min,max
		open(unit=1,file='datap1.f90')
		print*,''
		print*, '|----------------------------------------------|'
		print*,'ingreso de la matriz DATOS a(i,j)'
		read(1,*)m,n
		do i=1,m
			read(1,*)S(i,1:n)
			print 100,S(i,1:n)
		end do
		close(1)
		do i=1,10
			do j=1,10
				mat(i,j)=0;
			enddo
		enddo
		
		!hallando sumatorias
		print*, ''
		print*, '|----------------------------------------------|'
		print*,'	Matriz'
		do i=1,(opc+1)
			expo=i-1;
			do j=1,(opc+1)
				do k=1,m
					mat(i,j)=mat(i,j)+S(k,1)**expo;				
				enddo
				expo=expo+1;
			enddo
			expo=i-1;
			do k=1,m
				mat(i,(opc+2))=mat(i,(opc+2))+S(k,2)*S(k,1)**expo
			enddo
			print 100,mat(i,1:(opc+2))
		enddo

		!triangulo ceros
		do i=1,(opc)
			expo=i-1;
			do j=(i+1),(opc+2)
				cte=(-mat(j,i))/(mat(i,i))
				do k=i,(opc+2)
					mat(j,k)=((mat(i,k))*cte)+mat(j,k)
					!para hacer cero a valores muy pequeños
					if(mat(j,k)<0.000001)then
						mat(j,k)=0;
					endif		
				enddo
			enddo
		enddo

		!Mostrando
		print*, ''
		print*, '|----------------------------------------------|'
		print*,'	Matriz Reducida'
		do i=1,(opc+1)
			print 100,mat(i,1:(opc+2))
		enddo
		
		!obteniendo a0,a1,..
		x=(opc+2)

		do i=1,(x-1)
			
			if(i>1)then
				suba(7-i)=0
				k=0
				do j=(x-i+1),(x-1)
					k=k+1
					suba(7-i)=suba(7-i)-(suba(7-k)*mat((x-i),j))				
				enddo
				suba(7-i)=mat((x-i),x)+suba(7-i)
			else
			suba(7-i)=mat((x-i),x)
			endif
			suba(7-i)=suba(7-i)/mat((x-i),(x-i))
		enddo

		print*, ''
		print*, '|----------------------------------------------|'
		print*,'	Valores'
		print*,'a0=       ,a1=       ,a2=       ...'
		print 100,suba((8-x):6)	


		open(unit=1,file='salida.wxm')
		write(1,200)'/* [wxMaxima batch file version 1] [ DO NOT EDIT BY HAND! ]*/'
		write(1,200)'/* [ Created with wxMaxima version 16.04.2 ] */'
		write(1,200)
		write(1,200)'/* [wxMaxima: input   start ] */'

		write(1,*) 'wxdraw2d(color=blue,'

		write(1,*) 'explicit('
		j=0;
		do i=(8-x),6			
			if(i==6)then
				write(1,*) suba(i),'*x^',j
			else
				write(1,*) suba(i),'*x^',j,'+'
			endif
			j=j+1
		enddo

		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		do i=1,m
			if(i==1)then
				min=S(i,1)
				max=S(i,1)
			else
				if(max<S(i,1))then
					max=S(i,1)
				endif
				if(min>S(i,1))then
					min=S(i,1)
				endif
			endif			
		end do
		olgura=max-min
		max=max+abs(0.1*olgura)
		min=min-abs(0.1*olgura)
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		write(1,*) ',x,',min,',',max,'),'

		write(1,*) 'color=orange-red,'
		write(1,*) 'point_type=5,'
		write(1,*) 'point_size=1,'

		write(1,*) 'points(['
		do i=1,m
			write(1,*) '['
			do j=1,n
				if(j==n)then
					if(i==m)then
						write(1,*),S(i,j),']'
					else
						write(1,*),S(i,j),'],'
					endif
					else
					write(1,*),S(i,j),','
				endif				
			enddo
		end do
		write(1,*) ']),'

		write(1,*) 'xlabel="X",'
		write(1,*) 'ylabel="Y",'
		write(1,*) 'title="Aproximacion");'
		write(1,200)'/* [wxMaxima: input   end   ] */'
		write(1,200)
		write(1,200)'/* Maxima cant load/batch files which end with a comment! */'
		write(1,200)'Created with wxMaxima"$'
		close(1)

		100 FORMAT(F8.5,2X,F8.5,2X,F8.5,2X,F8.5,2X,F8.5,2X,F8.5,2X,F8.5)
		200 FORMAT(T1,A)
		
	End Subroutine grado
	
end module subrutinas

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!					contenedor de funciones					!


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!						menu principal						!
program mNumerico
	use subrutinas

	implicit none
	integer:: opc
	character:: cont

	156 continue
	print*, '|----------------------------------------------|'
	print*, '	TEORIA DE APROXIMACIONES'
	print*, '|----------------------------------------------|'
	print*, '	1 -APROXIMACION NORMAL'
	print*, '	2 -APROXIMACION GRADO 2'
	print*, '	3 -APROXIMACION GRADO 3'
	print*, '	4 -APROXIMACION GRADO 4'
	print*, '	5 -APROXIMACION GRADO 5'
	print*, '|----------------------------------------------|'
	read*, opc

	select case(opc)
		case(1)
			call normal
		case(2)
			call grado(opc)
		case(3)
			call grado(opc)
		case(4)
			call grado(opc)
		case(5)
			call grado(opc)
		case default
		print*, '¡La opcion que ingreso NO ES VALIDA!!!'

		goto 156
	End select
	print*, '|----------------------------------------------|'
	print*, '	DESEA OTRA VEZ EJECUTAR AL PROGRAMA S/N'
	Print*, ''
	read*, cont
	if (cont=='s')then
		goto 156
	End if

end program mNumerico