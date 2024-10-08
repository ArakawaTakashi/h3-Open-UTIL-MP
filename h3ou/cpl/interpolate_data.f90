subroutine interpolate_data(recv_model_name, send_model_name, mapping_tag, sn1, sn2, send_data, &
                            rn1, rn2, recv_data, num_of_data, tn, exchange_tag)
  use h3ou_api, only : h3ou_interpolation
  implicit none
  character(len=*), intent(IN) :: recv_model_name, send_model_name
  integer, intent(IN) :: mapping_tag
  integer, intent(IN) :: sn1, sn2
  real(kind=8), intent(IN) :: send_data(sn1,sn2)
  integer, intent(IN) :: rn1, rn2
  real(kind=8), intent(INOUT) :: recv_data(rn1,rn2)
  integer, intent(IN) :: num_of_data
  integer, intent(IN) :: tn
  integer, intent(IN) :: exchange_tag(tn)

  call h3ou_interpolation(recv_model_name, send_model_name, mapping_tag,&
       send_data, recv_data, num_of_data, exchange_tag(1))

end subroutine interpolate_data
