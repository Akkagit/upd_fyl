#Minimum requirement of CMake version : 3.0.0
cmake_minimum_required(VERSION 3.0.0)

#Project name and version number
project(${libVis} VERSION 0.1.0)

# set include path
include_directories(${libMatSrc})
include_directories(${libBaseSrc})
include_directories(${libRTKSrc})
include_directories(${libINSSrc})
include_directories(${libAmbSrc})
include_directories(${libLSQSrc})
include_directories(${libVisSrc})
if(USE_VIS)
	if(CMAKE_SYSTEM_NAME MATCHES "Windows")
		include_directories(${opencv_ROOT}/include)
		include_directories(${GLFW_ROOT}/include)
	else()
		include_directories(${GLFW_INCLUDE_DIRS})
		include_directories(${OpenCV_INCLUDE_DIRS})
	endif()
endif()
if(USE_VIS AND USE_CERES)
        if(CMAKE_SYSTEM_NAME MATCHES "Windows")
		include_directories(${ceres_ROOT}/ceres-solver/include)
		include_directories(${ceres_ROOT}/win/include)
		include_directories(${ceres_ROOT}/glog/src/windows)
	else()
		include_directories(${ceres_ROOT}/glog/src)
		include_directories(${ceres_ROOT}/ceres-solver/include)
	endif()
endif()

# for all the header files and source_files
file(GLOB_RECURSE header_files *.h)
file(GLOB_RECURSE source_files *.cpp)

if(CMAKE_SYSTEM_NAME MATCHES "Windows")
	# find files and put them into different groups
	file(GLOB gall_header_files 	gall/*.h)
	file(GLOB gcoders_header_files 	gcoders/*.h)
	file(GLOB gdata_header_files 	gdata/*.h)
	file(GLOB gio_header_files 		gio/*.h)
	file(GLOB gmodels_header_files 	gmodels/*.h)
	file(GLOB gproc_header_files 	gproc/*.h)
	file(GLOB gprod_header_files 	gprod/*.h)
	file(GLOB gproj_header_files 	gproj/*.h)
	file(GLOB gset_header_files 	gset/*.h)
	file(GLOB gutils_header_files 	gutils/*.h)
	file(GLOB gzstream_header_files gzstream/*.h)
	file(GLOB pugixml_header_files 	pugixml/*.hpp)
	file(GLOB md5_header_files 		md5/*.h)
	file(GLOB gins_header_files 	gins/*.h)
	file(GLOB gvis_header_files     gvis/*.h)

	file(GLOB gall_source_files 	gall/*.cpp)
	file(GLOB gcoders_source_files 	gcoders/*.cpp)
	file(GLOB gdata_source_files 	gdata/*.cpp)
	file(GLOB gio_source_files 		gio/*.cpp)
	file(GLOB gmodels_source_files 	gmodels/*.cpp)
	file(GLOB gproc_source_files 	gproc/*.cpp)
	file(GLOB gprod_source_files 	gprod/*.cpp)
	file(GLOB gproj_source_files 	gproj/*.cpp)
	file(GLOB gset_source_files 	gset/*.cpp)
	file(GLOB gutils_source_files 	gutils/*.cpp)
	file(GLOB gzstream_source_files gzstream/*.cpp)
	file(GLOB pugixml_source_files 	pugixml/*.cpp)
	file(GLOB md5_source_files	 	md5/*.cpp)
	file(GLOB gins_source_files 	gins/*.cpp)
	file(GLOB gvis_source_files     gvis/*.cpp)

	# add files to project
	source_group("CMake Files"  FILES CMakeLists.txt)
	source_group("Header Files" FILES header_files)
	source_group("Soruce Files" FILES source_files)
	
	source_group("Header Files\\gall"		FILES ${gall_header_files})
	source_group("Header Files\\gcoders"	FILES ${gcoders_header_files})
	source_group("Header Files\\gdata"		FILES ${gdata_header_files})
	source_group("Header Files\\gio"		FILES ${gio_header_files})
	source_group("Header Files\\gmodels"	FILES ${gmodels_header_files})
	source_group("Header Files\\gproc"		FILES ${gproc_header_files})
	source_group("Header Files\\gprod"		FILES ${gprod_header_files})
	source_group("Header Files\\gproj"		FILES ${gproj_header_files})
	source_group("Header Files\\gset"		FILES ${gset_header_files})
	source_group("Header Files\\gins"		FILES ${gins_header_files})
	source_group("Header Files\\gutils"		FILES ${gutils_header_files})
	source_group("Header Files\\gvis"       FILES ${gvis_header_files})
	
	source_group("Source Files\\gall"		FILES ${gall_source_files})
	source_group("Source Files\\gcoders"	FILES ${gcoders_source_files})
	source_group("Source Files\\gdata"		FILES ${gdata_source_files})
	source_group("Source Files\\gio"		FILES ${gio_source_files})
	source_group("Source Files\\gmodels"	FILES ${gmodels_source_files})
	source_group("Source Files\\gproc"		FILES ${gproc_source_files})
	source_group("Source Files\\gprod"		FILES ${gprod_source_files})
	source_group("Source Files\\gproj"		FILES ${gproj_source_files})
	source_group("Source Files\\gset"		FILES ${gset_source_files})
	source_group("Source Files\\gins"		FILES ${gins_source_files})	
	source_group("Source Files\\gvis"       FILES ${gvis_source_files})
else()
	message(STATUS "For Linux / unix, there is no need to group so detailed")
endif()

add_definitions(-D GREAT_VIS_LIBRARY)
add_library(${PROJECT_NAME} SHARED ${header_files} ${source_files})


if(CMAKE_SYSTEM_NAME MATCHES "Windows")
	include_directories(${zlib_ROOT})
	link_directories(${ROOT}/build/Lib/Debug)
	link_directories(${ROOT}/build/Lib/Release)
	link_directories(${ROOT}/build/Lib/RelWithDebInfo)
	link_directories(${ROOT}/build/Lib/MinSizeRel)
	target_link_libraries(${PROJECT_NAME} ${libBase})
	target_link_libraries(${PROJECT_NAME} ${libMat})
	target_link_libraries(${PROJECT_NAME} ${libAmb})
	if(USE_VIS)
		target_link_libraries(${PROJECT_NAME} 
			debug ${opencv_ROOT}/lib/opencv_world401d.lib 
			${GLFW_ROOT}/lib/glfw3.lib
			glu32.lib 
			opengl32.lib 
			#vcruntimed.lib 
			#ucrtd.lib
			#msvcrtd.lib
			optimized ${opencv_ROOT}/lib/opencv_world401.lib
			${GLFW_ROOT}/lib/glfw3.lib
			glu32.lib opengl32.lib
			#msvcrt.lib
			#ucrt.lib
			#vcruntime.lib
		)
	endif()
	if(USE_VIS AND USE_CERES)
		target_link_libraries(${PROJECT_NAME} 
			debug ${ceres_ROOT}/x64/Release/ceres.lib			
			${ceres_ROOT}/x64/Release/ceres_static.lib
			${ceres_ROOT}/x64/Release/libglog_static.lib
			optimized ${ceres_ROOT}/x64/Release/ceres.lib			
			${ceres_ROOT}/x64/Release/ceres_static.lib
			${ceres_ROOT}/x64/Release/libglog_static.lib
		)
	endif()
else()
	link_directories(${ROOT}/build/Lib)
	target_link_libraries(${PROJECT_NAME} ${libBase})
	target_link_libraries(${PROJECT_NAME} ${libMat})
	target_link_libraries(${PROJECT_NAME} ${libAmb})
	if(USE_VIS)
		target_link_libraries(${PROJECT_NAME} ${GLFW_LIBRARIES} ${OpenCV_LIBS} ${OPENGL_gl_LIBRARY})
	endif()
endif()

add_dependencies(${PROJECT_NAME} ${libBase} ${libMat} ${libAmb} ${libLSQ} )
