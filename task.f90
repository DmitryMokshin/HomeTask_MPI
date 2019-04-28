module Task
use mpi
  contains

  subroutine GetMaxCoordinates(A,x1,y1,x2,y2)
    implicit none
      real(8),dimension(:,:), intent(in) :: A
      integer(4),intent(out) :: x1,y1,x2,y2
      integer(4) :: n,left,right,up,down,m,tr,minpos,i
      real(8),allocatable :: curcol(:),B(:,:)
      real(8) :: cursum,maxsum,ressum
      logical :: transpos
      
  
       m=size(A,dim=1)
       n=size(A,dim=2)
       transpos=.FALSE.
        if (m<n) then
        
           transpos=.TRUE.
           B=transpose(A)
           m=size(B,dim=1)
           n=size(B,dim=2)
           
        else
        
           B=A
           
        endif
        
       allocate(curcol(m))
       	
       maxsum=B(1,1)
       x1=1; y1=1; x2=1; y2=1;

        do left=1, n
        
           curcol=B(:,left)
           
             do right=left,n

                 if (right>left) then
                    curcol=curcol+B(:,right)
                 endif

                cursum=curcol(1)
                up=1; down=1; ressum=0; minpos=0;
                
                  do i=1, size(curcol)                  
                     ressum=ressum+curcol(i)
                     
                      if (ressum > cursum) then
                         cursum=ressum
                         up=minpos+1
                         down=i
                      endif

                      if (ressum < 0) then
                         ressum=0 
                         minpos=i
                      endif
                      
                  enddo
               
          if (cursum > maxsum) then
          
             maxsum=cursum
              x1=up
              x2=down
              y1=left
              y2=right
              
          end if
             enddo
             
        enddo



       deallocate(curcol)

        if (transpos) then
           tr=x1
           x1=y1
           y1=tr

           tr=y2
           y2=x2
           x2=tr
           
        endif

  
  end subroutine GetMaxCoordinates


end module Task
