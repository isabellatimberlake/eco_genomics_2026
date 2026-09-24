# Playing in R ####
# 4 hashtags adds a section! appears at bottom of script (or cntrl+shift+R)

x <- 5
students <- data.frame(
  name = c("A", "B", "C", "D", "E", "F"),
  height = c(62, 68, 72, 64, 69, 68)
)

head(students)

# next section ####

class(students) # tells what the variable is (str, data.frame, etc)
str(students) # description of the data frame
students$name[2] # indexing into the data frame column name
students[1, 2] # indexing into the whole data frame (row 1, col 2)

mean(students$height) # 67.16667
