module colours
  use types, only: WP
  implicit none
  !!! RGB colours in [0-1, 0-1, 0-1] format
  ! Tableau10
  real(kind=WP), dimension(3), parameter :: col_blue = (/ 0.3411764705882353_WP, 0.47058823529411764_WP, 0.6431372549019608_WP/)
  real(kind=WP), dimension(3), parameter :: col_orange = (/ 0.8941176470588236_WP, 0.5803921568627451_WP, 0.26666666666666666_WP/)
  real(kind=WP), dimension(3), parameter :: col_red = (/ 0.8196078431372549_WP, 0.3803921568627451_WP, 0.36470588235294116_WP/)
  real(kind=WP), dimension(3), parameter :: col_teal = (/ 0.5215686274509804_WP, 0.7137254901960784_WP, 0.6980392156862745_WP/)
  real(kind=WP), dimension(3), parameter :: col_green = (/ 0.41568627450980394_WP, 0.6235294117647059_WP, 0.34509803921568627_WP/)
  real(kind=WP), dimension(3), parameter :: col_yellow = (/ 0.9058823529411765_WP, 0.792156862745098_WP, 0.3764705882352941_WP/)
  real(kind=WP), dimension(3), parameter :: col_purple = (/ 0.6588235294117647_WP, 0.48627450980392156_WP, 0.6235294117647059_WP/)
  real(kind=WP), dimension(3), parameter :: col_pink = (/ 0.9450980392156862_WP, 0.6352941176470588_WP, 0.6627450980392157_WP/)
  real(kind=WP), dimension(3), parameter :: col_brown = (/ 0.5882352941176471_WP, 0.4627450980392157_WP, 0.3843137254901961_WP/)
  real(kind=WP), dimension(3), parameter :: col_grey = (/ 0.7215686274509804_WP, 0.6901960784313725_WP, 0.6745098039215687_WP/)
  ! Other
  real(kind=WP), dimension(3), parameter :: col_black = (/ 0.0_WP, 0.0_WP, 0.0_WP/)

end module colours
