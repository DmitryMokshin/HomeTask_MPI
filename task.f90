module Task
use :: mpi
  contains

  subroutine GetMaxCoordinates(A,x1,y1,x2,y2)
    implicit none
      real(8),dimension(:,:), intent(in) :: A
      integer(4),intent(out) :: x1,y1,x2,y2
      integer(4) :: n,left,right,up,down,m,tr,minpos,i,k
      integer(4) :: mpiErr, mpiSize, mpiRank
      real(8),allocatable :: curcol(:),B(:,:)
      real(8) :: maxsum
      real(8) :: cursum,ressum,resmaxsum
      logical :: transpos

       call mpi_init(mpiErr)
       call mpi_comm_size(MPI_COMM_WORLD, mpiSize, mpiErr)
       call mpi_comm_rank(MPI_COMM_WORLD, mpiRank, mpiErr)
       
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
       x1=1;x2=1;y1=1;y2=1; maxsum=B(1,1)
        do left=1+mpiRank,n,mpiSize

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

       call mpi_finalize(mpiErr)

  end subroutine GetMaxCoordinates


end module Task
