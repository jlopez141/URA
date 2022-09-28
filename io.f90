module inout


contains

subroutine start_message()
implicit NONE
write(*,fmt='(a)') ""
write(*,fmt='(a)') ""
write(*,fmt='(a)') "       **************************************************"
write(*,fmt='(a)') "       **                                              **"
write(*,fmt='(a)') "       **         Ultimate Rigidity Analyzer           **"
write(*,fmt='(a)') "       **                    URA                      **"
write(*,fmt='(a)') "       **                   v 0.1                      **"
write(*,fmt='(a)') "       **                                              **"
write(*,fmt='(a)') "       **                  created by                  **"
write(*,fmt='(a)') "       **    J. Lopez Zorrilla & X. M. Aretxabaleta    **"
write(*,fmt='(a)') "       **                                              **"
write(*,fmt='(a)') "       **************************************************"
write(*,fmt='(a)') ""
write(*,fmt='(a)') ""

end subroutine
end module
