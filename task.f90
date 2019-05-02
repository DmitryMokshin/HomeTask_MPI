module Task
use::mpi
  contains

  subroutine GetMaxCoordinates(A,x1,y1,x2,y2)
   implicit none
    real(8), dimension(:,:),intent(in) :: A
    integer(4), intent(out) :: x1,y1,x2,y2
    integer(4) :: n,left,right,up,down,m,tr,min_pos,i
    integer(4) :: mpiErr,mpiSize,mpiRank,local_coord(4)
    real(8),allocatable :: current_column(:),B(:,:)
    real(8),allocatable :: local_max_sum_vector(:),global_max_sum(:)
    real(8) :: current_sum,inter_sum,local_max_sum
    logical :: transpos

     call mpi_init(mpiErr)
     call mpi_comm_size(MPI_COMM_WORLD,mpiSize,mpiErr)
     call mpi_comm_rank(MPI_COMM_WORLD,mpiRank,mpiErr)

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

     allocate(current_column(m))
     local_coord=1; local_max_sum=B(1,1)

         do left=1+mpiRank,n,mpiSize
           current_column=B(:,left)

             do right=left,n

                 if (right>left) then
                   current_column=current_column+B(:,right)
                 endif

               current_sum=current_column(1)
               up=1; down=1; inter_sum=0; min_pos=0;

                 do i=1,size(current_column)
                   inter_sum=inter_sum+current_column(i)

                     if (inter_sum>current_sum) then
                       current_sum=inter_sum
                       up=min_pos+1
                       down=i
                     endif

                     if (inter_sum<0) then
                       inter_sum=0
                       min_pos=i
                     endif
                 enddo

                 if (current_sum>local_max_sum) then
                   local_max_sum=current_sum
                   local_coord(1)=up
                   local_coord(2)=down
                   local_coord(3)=left
                   local_coord(4)=right
                 endif
             enddo
         enddo

     allocate(local_max_sum_vector(0:mpiSize-1),global_max_sum(0:mpiSize-1))

     local_max_sum_vector(mpiRank)=local_max_sum
     call mpi_allreduce(local_max_sum_vector,global_max_sum,mpiSize,mpi_real8,MPI_MAX,MPI_COMM_WORLD,mpiErr)
     call mpi_bcast(local_coord,4,mpi_integer,maxloc(global_max_sum,dim=1)-1,MPI_COMM_WORLD,mpiErr)

     deallocate(local_max_sum_vector,global_max_sum)
     deallocate(current_column)

     x1=local_coord(1)
     x2=local_coord(2)
     y1=local_coord(3)
     y2=local_coord(4)

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
