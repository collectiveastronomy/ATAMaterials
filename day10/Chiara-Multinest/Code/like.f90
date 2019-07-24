module like
	
use params
use utils1
implicit none
      
contains
      
      
!=======================================================================

subroutine slikelihood(Cube,slhood)
         
	implicit none
      
	double precision Cube(nest_nPar),slhood
	integer i,j
	
	slhood = 0.0d0
	
	do i = 1, sdim
		Cube(i) = spriorran(i, 1) + ( spriorran(i, 2) - spriorran(i, 1) ) * Cube(i)
	enddo
	
	call nrmhd(slhood,Cube,sdim)

        slhood = log(slhood)

end subroutine slikelihood
      
!=======================================================================

end module like
