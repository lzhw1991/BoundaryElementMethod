        subroutine gasave(ISTEP,ID,GA1,GA2,N1,N2)
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c program to save the GAMtrix to hard drive
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        complex GA1(N1,N2),GA2(N2,N2)
        character(77) fname1,fname2
           if(id.lt.10) then
             write(fname1,'(A5,I1,A4)')'gatdL',id,'.dat'
             write(fname2,'(A5,I1,A4)')'gaddL',id,'.dat'
           else
             write(fname1,'(A5,I2,A4)')'gatdL',id,'.dat'
             write(fname2,'(A5,I2,A4)')'gaddL',id,'.dat'
           endif
         open(50+id,file=fname1,form='unformatted')
        open(30+id,file=fname2,form='unformatted')
        WRITE(30+ID) GA1
        WRITE(50+ID) GA2
        close(30+ID)
        close(50+ID)
        END         
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        subroutine gaload(ISTEP,ID,GA,N)
        complex ga(N,N)
        character(77) fname1
           if(id.lt.10) then
             write(fname1,'(A5,I1,A4)')'gatdL',id,'.dat'
           else
             write(fname1,'(A5,I2,A4)')'gatdL',id,'.dat'
           endif
       open(50+id,file=fname1,form='unformatted')
        read(50+ID) GA
        close(50+ID)
        end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        subroutine gaload1(ISTEP,id,ga,N1,N2)
        complex ga(N1,n2)
        character(77) fname2

           if(id.lt.10) then
             write(fname2,'(A5,I1,A4)')'gaddL',id,'.dat'
           else
             write(fname2,'(A5,I2,A4)')'gaddL',id,'.dat'
           endif
        open(30+id,file=fname2,form='unformatted')
        read(30+ID) ga
        close(30+ID)
        end
