module typic_pattern 
    
    !USE - the command use, or use, only must be before implicit none
    !the use command can be used to share entities or data

    implicit none
    
    !you can declare any data that you want to be 'globally' 
    !visible to all the procedures defined later


contains !marks the end of the main code of the program and 
         !indicates the start of the definition of functions and subroutines
   
    !add comments to the top of a block of expressions 
    !explaining how the following code relates to the physical problem.

    !functions, subroutines...

        !if you will use loops add a comment for what you will use it

        !check if there is any division by zero


!function name(arg1, arg2, ....) 
!    [declarations, including those for the arguments]
!    [executable statements] 
!end function [name]

!subroutine name(arg1, arg2, ....) 
!    [declarations, including those for the arguments]
!    [executable statements] 
!end subroutine [name]

end module teste