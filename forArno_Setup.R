#For Arno, getting the libraries and other files setup

# 1. libraries
# Installing packages can sometimes generate odd errors. 
# I recommend running one at a time. 
# You'll know it worked when the console displays no errors, followed by a new line starting with ">"
# If you get errors, it may work to simply run that line again, or to add "dependencies = T", such as...
# install.packages('R.matlab', dependencies = T)

install.packages('R.matlab')
install.packages('dplyr')
install.packages('tidyr')
install.packages('mgcv')
install.packages('gratia')
install.packages('ggplot2')

# 2. other files
# These 3 files need to be in the working directory
# 'utils.R', 'trance_area_power_fixed.csv', 'area_coordinates.csv',  