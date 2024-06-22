project (speech-detection)

include_directories(${PROJECT_SOURCE_DIR}/Inc)

aux_source_directory(${PROJECT_SOURCE_DIR}/Src SRC_LIST)

add_executable(mfcc ${SRC_LIST})