module colours
  use types, only: WP
  implicit none
  !!! RGB colours in [0-1, 0-1, 0-1] format
  ! Tableau10
  real(kind=WP), dimension(3), parameter :: col_blue = (/ 0.3411764705882353, 0.47058823529411764, 0.6431372549019608/)
  real(kind=WP), dimension(3), parameter :: col_orange = (/ 0.8941176470588236, 0.5803921568627451, 0.26666666666666666/)
  real(kind=WP), dimension(3), parameter :: col_red = (/ 0.8196078431372549, 0.3803921568627451, 0.36470588235294116/)
  real(kind=WP), dimension(3), parameter :: col_teal = (/ 0.5215686274509804, 0.7137254901960784, 0.6980392156862745/)
  real(kind=WP), dimension(3), parameter :: col_green = (/ 0.41568627450980394, 0.6235294117647059, 0.34509803921568627/)
  real(kind=WP), dimension(3), parameter :: col_yellow = (/ 0.9058823529411765, 0.792156862745098, 0.3764705882352941/)
  real(kind=WP), dimension(3), parameter :: col_purple = (/ 0.6588235294117647, 0.48627450980392156, 0.6235294117647059/)
  real(kind=WP), dimension(3), parameter :: col_pink = (/ 0.9450980392156862, 0.6352941176470588, 0.6627450980392157/)
  real(kind=WP), dimension(3), parameter :: col_brown = (/ 0.5882352941176471, 0.4627450980392157, 0.3843137254901961/)
  real(kind=WP), dimension(3), parameter :: col_grey = (/ 0.7215686274509804, 0.6901960784313725, 0.6745098039215687/)
  ! Other
  real(kind=WP), dimension(3), parameter :: col_black = (/ 0.0_WP, 0.0_WP, 0.0_WP/)

end module colours
