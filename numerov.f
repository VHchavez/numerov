!____________________________________________________________________________________!
!						                                     !
!  		    			  	                                     !
!                                Física Computacional                                !
!                                   Noviembre 2014                                   !
!                                  V.H.G.Chávez                                      !
!                                                                                    !
!                           Dada la ecuación de Schrödinger                          !
!                    Resuelve el problema del Oscilador Armónico                     !
!                                                                                    !
! 						                                     !                     
!                                  d_x[si]+k^2[si]=0   		                     !
!____________________________________________________________________________________!

       


        program Quantum Oscillator

        integer m, r
        double precision h,h2,E,y,k,x,k_0,k_1,psi_0,psi_1,psi,b,q,t
        double precision energias, error, contador
        integer, parameter ::   N=500     !Iteraciones
        double precision, dimension(-N:N) :: si,xi


        OPEN(10, FIlE='psi.txt')
        OPEN(11, FIlE='psi2.txt')



!____________________________________________________________________________________!
!				Condiciones Iniciales                                !  
!____________________________________________________________________________________!


        energias  =           7                    !Cantidad de Soluciones a Encontrar
        contador  =           0
   
        h         =         .01                    !Incrementos
        !El producto de N*h resulta en el intervalo de integración de psi.

        error     =          .00000000000002       !Incertidumbre aceptada para el valor de Energía. 
        E         =           0.0                    !En dónde comienza la búsqueda de energías. 
        ep        =           0.1                  !Incremento de la Búsqueda
        t         =           0.0
        q         =           0


!____________________________________________________________________________________!
!	                  	Método de Numerov                                    !  
!____________________________________________________________________________________!

        do while (contador .lt. energias) 
         
            !Condiciones iniciales.
            y      =         0.0                    
            k      =         0.0                    
            x      =  -1*(N+2)*h


            !Definimos los parámetros para el algoritmo de Numerov.
            k_0    =           E + x-2*h
            k_1    =           E + x-h
            psi_0  =           0
            psi_1  =           0.00000000000000001
             
        
            q = 0
            do m=-1*N+2, N-2, 1
               x     = x+h 
               xi(m) = x
               k     = 2*E - x**2
               b     = h**2/12
               psi   = (2*(1-5*b*k_1)*psi_1-(1+b*k_0)*psi_0)/(1+b*k)
               si(m) = psi

               if (q .lt. si(m)) then
                   q = si(m)
               end if

               psi_0 = psi_1
               psi_1 = psi
               k_0   = k_1
               k_1   = k

               
            end do
 
!____________________________________________________________________________________!
!	                  	  Búsqueda de Energías                               !  
!____________________________________________________________________________________!
     
             if (abs(si(N-2)) .lt. error) then
                 contador = contador +1
                 do r=-1*N+2, N-2,1


                      WRITE(10,*) xi(r), si(r)/q + t + 0.5
                      WRITE(11,*) xi(r), (si(r)/q)**2 + t + 0.5

                     
                 end do

                      WRITE(10,*) 
                      WRITE(10,*) 
                      WRITE(11,*) 
                      WRITE(11,*) 

                      WRITE(*,*) "Su valor de Energía es: ", E
2000   format(f5.3)

                      t=t+2
                      
             end if
            E = E+ep
            q = q+1


        end do


        return
        STOP
        END
