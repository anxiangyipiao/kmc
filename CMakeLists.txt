cmake_minimum_required(VERSION 3.10)

# 设置项目名称和语言
project(kmcplus CXX)

# 指定 C++ 标准
#set(CMAKE_CXX_STANDARD 17)
#set(CMAKE_CXX_STANDARD_REQUIRED ON)
#set(CMAKE_CXX_EXTENSIONS OFF)

# 编译器选项
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -march=native")

# 第三方库头文件和库文件的路径
set(BLITZ_INCLUDE_DIR "/home/liaosi/local/blitz/include")
set(YAML_CPP_INCLUDE_DIR "/home/liaosi/local/yaml-cpp/include")
set(YAML_CPP_LIBRARY "/home/liaosi/local/yaml-cpp/lib/libyaml-cpp.so.0.6")

# 添加头文件搜索路径
include_directories(${BLITZ_INCLUDE_DIR} ${YAML_CPP_INCLUDE_DIR})

# 添加源代码
file(GLOB_RECURSE SOURCES "src/*.cpp")

# 生成可执行文件
add_executable(kmcplus.x ${SOURCES})

# 链接第三方库
target_link_libraries(kmcplus.x ${YAML_CPP_LIBRARY})

# 安装目标
install(TARGETS kmcplus.x DESTINATION ${CMAKE_SOURCE_DIR}/bin)
