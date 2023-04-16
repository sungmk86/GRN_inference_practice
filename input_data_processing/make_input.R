# import functions
setwd('/home/seongwonhwang/Desktop/projects/git/GRN_inference_practice/input_data_processing')
source("make_input_utils.R")
TEST_ID = 'TEST3'

path_input <- "/home/seongwonhwang/Desktop/projects/mogrify/Statistical\ Consulting/"
path_tf_and_reqdgenes <- "/home/seongwonhwang/Desktop/projects/GRN_in_general/PyG/data/Anonymized_tables/Anonymized_tables/"
path_output <- "/home/seongwonhwang/Desktop/projects/git/GRN_inference_practice/input_data_processing/data"
#############################################
MI <- MakeInput$new(TEST_ID, path_input, path_tf_and_reqdgenes, path_output)
MI$write_files("training")
MI$write_files("predicting")
