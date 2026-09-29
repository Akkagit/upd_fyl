#Minimum requirement of CMake version : 3.0.0
cmake_minimum_required(VERSION 3.0.0)

#Project name and version number
project(${libLidar} VERSION 0.1.0)

# set include path
include_directories(${libMatSrc})
include_directories(${libBaseSrc})
include_directories(${libRTKSrc})
include_directories(${libINSSrc})
include_directories(${libAmbSrc})
include_directories(${libLSQSrc})
include_directories(${libLidarSrc})

if(USE_LIDAR)
	include_directories(${opencv_ROOT}/include)
	include_directories(${GLFW_ROOT}/include)
	include_directories(${PCL_ROOT}/include)
endif()
if(USE_LIDAR AND USE_CERES)
	include_directories(${ceres_ROOT}/ceres-solver/include)
	include_directories(${ceres_ROOT}/win/include)
	include_directories(${ceres_ROOT}/glog/src/windows)
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
	file(GLOB glidar_header_files   glidar/*.h)

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
	file(GLOB glidar_source_files   glidar/*.cpp)

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
	source_group("Header Files\\glidar"     FILES ${glidar_header_files})
	
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
	source_group("Source Files\\glidar"     FILES ${glidar_source_files})
else()
	message(STATUS "For Linux / unix, there is no need to group so detailed")
endif()

add_definitions(-D GREAT_LIDAR_LIBRARY)
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
	if(USE_LIDAR)
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
			${pcl_ROOT}/lib/pcl_common_release.lib
			${pcl_ROOT}/lib/pcl_features_release.lib
			${pcl_ROOT}/lib/pcl_filters_release.lib
			${pcl_ROOT}/lib/pcl_io_ply_release.lib
			${pcl_ROOT}/lib/pcl_io_release.lib
			${pcl_ROOT}/lib/pcl_kdtree_release.lib
			${pcl_ROOT}/lib/pcl_keypoints_release.lib
			${pcl_ROOT}/lib/pcl_ml_release.lib
			${pcl_ROOT}/lib/pcl_octree_release.lib
			${pcl_ROOT}/lib/pcl_outofcore_release.lib
			${pcl_ROOT}/lib/pcl_people_release.lib
			${pcl_ROOT}/lib/pcl_recognition_release.lib
			${pcl_ROOT}/lib/pcl_registration_release.lib
			${pcl_ROOT}/lib/pcl_sample_consensus_release.lib
			${pcl_ROOT}/lib/pcl_search_release.lib
			${pcl_ROOT}/lib/pcl_segmentation_release.lib
			${pcl_ROOT}/lib/pcl_stereo_release.lib
			${pcl_ROOT}/lib/pcl_surface_release.lib
			${pcl_ROOT}/lib/pcl_tracking_release.lib
			${pcl_ROOT}/lib/pcl_visualization_release.lib
			${boost_ROOT}/lib/libboost_atomic-vc141-mt-1_64.lib
			${boost_ROOT}/lib/libboost_bzip2-vc141-mt-1_64.lib
			${boost_ROOT}/lib/libboost_chrono-vc141-mt-1_64.lib
			${boost_ROOT}/lib/libboost_container-vc141-mt-1_64.lib
			${boost_ROOT}/lib/libboost_context-vc141-mt-1_64.lib
			${boost_ROOT}/lib/libboost_coroutine-vc141-mt-1_64.lib
			${boost_ROOT}/lib/libboost_date_time-vc141-mt-1_64.lib
			${boost_ROOT}/lib/libboost_exception-vc141-mt-1_64.lib
			${boost_ROOT}/lib/libboost_fiber-vc141-mt-1_64.lib
			${boost_ROOT}/lib/libboost_filesystem-vc141-mt-1_64.lib
			${boost_ROOT}/lib/libboost_graph-vc141-mt-1_64.lib
			${boost_ROOT}/lib/libboost_graph_parallel-vc141-mt-1_64.lib
			${boost_ROOT}/lib/libboost_iostreams-vc141-mt-1_64.lib
			${boost_ROOT}/lib/libboost_locale-vc141-mt-1_64.lib
			${boost_ROOT}/lib/libboost_log-vc141-mt-1_64.lib
			${boost_ROOT}/lib/libboost_log_setup-vc141-mt-1_64.lib
			${boost_ROOT}/lib/libboost_math_c99-vc141-mt-1_64.lib
			${boost_ROOT}/lib/libboost_math_c99f-vc141-mt-1_64.lib
			${boost_ROOT}/lib/libboost_math_c99l-vc141-mt-1_64.lib
			${boost_ROOT}/lib/libboost_math_tr1-vc141-mt-1_64.lib
			${boost_ROOT}/lib/libboost_math_tr1f-vc141-mt-1_64.lib
			${boost_ROOT}/lib/libboost_math_tr1l-vc141-mt-1_64.lib
			${boost_ROOT}/lib/libboost_mpi-vc141-mt-1_64.lib
			${boost_ROOT}/lib/libboost_numpy-vc141-mt-1_64.lib
			${boost_ROOT}/lib/libboost_numpy3-vc141-mt-1_64.lib
			${boost_ROOT}/lib/libboost_prg_exec_monitor-vc141-mt-1_64.lib
			${boost_ROOT}/lib/libboost_program_options-vc141-mt-1_64.lib
			${boost_ROOT}/lib/libboost_python-vc141-mt-1_64.lib
			${boost_ROOT}/lib/libboost_python3-vc141-mt-1_64.lib
			${boost_ROOT}/lib/libboost_random-vc141-mt-1_64.lib
			${boost_ROOT}/lib/libboost_regex-vc141-mt-1_64.lib
			${boost_ROOT}/lib/libboost_serialization-vc141-mt-1_64.lib
			${boost_ROOT}/lib/libboost_signals-vc141-mt-1_64.lib
			${boost_ROOT}/lib/libboost_system-vc141-mt-1_64.lib
			${boost_ROOT}/lib/libboost_test_exec_monitor-vc141-mt-1_64.lib
			${boost_ROOT}/lib/libboost_thread-vc141-mt-1_64.lib
			${boost_ROOT}/lib/libboost_timer-vc141-mt-1_64.lib
			${boost_ROOT}/lib/libboost_type_erasure-vc141-mt-1_64.lib
			${boost_ROOT}/lib/libboost_unit_test_framework-vc141-mt-1_64.lib
			${boost_ROOT}/lib/libboost_wave-vc141-mt-1_64.lib
			${boost_ROOT}/lib/libboost_wserialization-vc141-mt-1_64.lib
			${boost_ROOT}/lib/libboost_zlib-vc141-mt-1_64.lib
			${flann_ROOT}/lib/flann.lib
			${flann_ROOT}/lib/flann_cpp.lib
			${flann_ROOT}/lib/flann_cpp_s.lib
			${flann_ROOT}/lib/flann_s.lib
			${openNI2_ROOT}/Lib/OpenNI2.lib
			${qhull_ROOT}/lib/qhull.lib
			${qhull_ROOT}/lib/qhullcpp.lib
			${qhull_ROOT}/lib/qhullstatic.lib
			${qhull_ROOT}/lib/qhullstatic_r.lib
			${qhull_ROOT}/lib/qhull_p.lib
			${qhull_ROOT}/lib/qhull_r.lib
			${vtk_ROOT}/lib/vtkalglib-8.0.lib
			${vtk_ROOT}/lib/vtkChartsCore-8.0.lib
			${vtk_ROOT}/lib/vtkCommonColor-8.0.lib
			${vtk_ROOT}/lib/vtkCommonComputationalGeometry-8.0.lib
			${vtk_ROOT}/lib/vtkCommonCore-8.0.lib
			${vtk_ROOT}/lib/vtkCommonDataModel-8.0.lib
			${vtk_ROOT}/lib/vtkCommonExecutionModel-8.0.lib
			${vtk_ROOT}/lib/vtkCommonMath-8.0.lib
			${vtk_ROOT}/lib/vtkCommonMisc-8.0.lib
			${vtk_ROOT}/lib/vtkCommonSystem-8.0.lib
			${vtk_ROOT}/lib/vtkCommonTransforms-8.0.lib
			${vtk_ROOT}/lib/vtkDICOMParser-8.0.lib
			${vtk_ROOT}/lib/vtkDomainsChemistry-8.0.lib
			${vtk_ROOT}/lib/vtkexoIIc-8.0.lib
			${vtk_ROOT}/lib/vtkexpat-8.0.lib
			${vtk_ROOT}/lib/vtkFiltersAMR-8.0.lib
			${vtk_ROOT}/lib/vtkFiltersCore-8.0.lib
			${vtk_ROOT}/lib/vtkFiltersExtraction-8.0.lib
			${vtk_ROOT}/lib/vtkFiltersFlowPaths-8.0.lib
			${vtk_ROOT}/lib/vtkFiltersGeneral-8.0.lib
			${vtk_ROOT}/lib/vtkFiltersGeneric-8.0.lib
			${vtk_ROOT}/lib/vtkFiltersGeometry-8.0.lib
			${vtk_ROOT}/lib/vtkFiltersHybrid-8.0.lib
			${vtk_ROOT}/lib/vtkFiltersHyperTree-8.0.lib
			${vtk_ROOT}/lib/vtkFiltersImaging-8.0.lib
			${vtk_ROOT}/lib/vtkFiltersModeling-8.0.lib
			${vtk_ROOT}/lib/vtkFiltersParallel-8.0.lib
			${vtk_ROOT}/lib/vtkFiltersParallelImaging-8.0.lib
			${vtk_ROOT}/lib/vtkFiltersPoints-8.0.lib
			${vtk_ROOT}/lib/vtkFiltersProgrammable-8.0.lib
			${vtk_ROOT}/lib/vtkFiltersSelection-8.0.lib
			${vtk_ROOT}/lib/vtkFiltersSMP-8.0.lib
			${vtk_ROOT}/lib/vtkFiltersSources-8.0.lib
			${vtk_ROOT}/lib/vtkFiltersStatistics-8.0.lib
			${vtk_ROOT}/lib/vtkFiltersTexture-8.0.lib
			${vtk_ROOT}/lib/vtkFiltersTopology-8.0.lib
			${vtk_ROOT}/lib/vtkFiltersVerdict-8.0.lib
			${vtk_ROOT}/lib/vtkfreetype-8.0.lib
			${vtk_ROOT}/lib/vtkGeovisCore-8.0.lib
			${vtk_ROOT}/lib/vtkgl2ps-8.0.lib
			${vtk_ROOT}/lib/vtkhdf5-8.0.lib
			${vtk_ROOT}/lib/vtkhdf5_hl-8.0.lib
			${vtk_ROOT}/lib/vtkImagingColor-8.0.lib
			${vtk_ROOT}/lib/vtkImagingCore-8.0.lib
			${vtk_ROOT}/lib/vtkImagingFourier-8.0.lib
			${vtk_ROOT}/lib/vtkImagingGeneral-8.0.lib
			${vtk_ROOT}/lib/vtkImagingHybrid-8.0.lib
			${vtk_ROOT}/lib/vtkImagingMath-8.0.lib
			${vtk_ROOT}/lib/vtkImagingMorphological-8.0.lib
			${vtk_ROOT}/lib/vtkImagingSources-8.0.lib
			${vtk_ROOT}/lib/vtkImagingStatistics-8.0.lib
			${vtk_ROOT}/lib/vtkImagingStencil-8.0.lib
			${vtk_ROOT}/lib/vtkInfovisCore-8.0.lib
			${vtk_ROOT}/lib/vtkInfovisLayout-8.0.lib
			${vtk_ROOT}/lib/vtkInteractionImage-8.0.lib
			${vtk_ROOT}/lib/vtkInteractionStyle-8.0.lib
			${vtk_ROOT}/lib/vtkInteractionWidgets-8.0.lib
			${vtk_ROOT}/lib/vtkIOAMR-8.0.lib
			${vtk_ROOT}/lib/vtkIOCore-8.0.lib
			${vtk_ROOT}/lib/vtkIOEnSight-8.0.lib
			${vtk_ROOT}/lib/vtkIOExodus-8.0.lib
			${vtk_ROOT}/lib/vtkIOExport-8.0.lib
			${vtk_ROOT}/lib/vtkIOExportOpenGL-8.0.lib
			${vtk_ROOT}/lib/vtkIOGeometry-8.0.lib
			${vtk_ROOT}/lib/vtkIOImage-8.0.lib
			${vtk_ROOT}/lib/vtkIOImport-8.0.lib
			${vtk_ROOT}/lib/vtkIOInfovis-8.0.lib
			${vtk_ROOT}/lib/vtkIOLegacy-8.0.lib
			${vtk_ROOT}/lib/vtkIOLSDyna-8.0.lib
			${vtk_ROOT}/lib/vtkIOMINC-8.0.lib
			${vtk_ROOT}/lib/vtkIOMovie-8.0.lib
			${vtk_ROOT}/lib/vtkIONetCDF-8.0.lib
			${vtk_ROOT}/lib/vtkIOParallel-8.0.lib
			${vtk_ROOT}/lib/vtkIOParallelXML-8.0.lib
			${vtk_ROOT}/lib/vtkIOPLY-8.0.lib
			${vtk_ROOT}/lib/vtkIOSQL-8.0.lib
			${vtk_ROOT}/lib/vtkIOTecplotTable-8.0.lib
			${vtk_ROOT}/lib/vtkIOVideo-8.0.lib
			${vtk_ROOT}/lib/vtkIOXML-8.0.lib
			${vtk_ROOT}/lib/vtkIOXMLParser-8.0.lib
			${vtk_ROOT}/lib/vtkjpeg-8.0.lib
			${vtk_ROOT}/lib/vtkjsoncpp-8.0.lib
			${vtk_ROOT}/lib/vtklibharu-8.0.lib
			${vtk_ROOT}/lib/vtklibxml2-8.0.lib
			${vtk_ROOT}/lib/vtklz4-8.0.lib
			${vtk_ROOT}/lib/vtkmetaio-8.0.lib
			${vtk_ROOT}/lib/vtkNetCDF-8.0.lib
			${vtk_ROOT}/lib/vtknetcdf_c++.lib
			${vtk_ROOT}/lib/vtkoggtheora-8.0.lib
			${vtk_ROOT}/lib/vtkParallelCore-8.0.lib
			${vtk_ROOT}/lib/vtkpng-8.0.lib
			${vtk_ROOT}/lib/vtkproj4-8.0.lib
			${vtk_ROOT}/lib/vtkRenderingAnnotation-8.0.lib
			${vtk_ROOT}/lib/vtkRenderingContext2D-8.0.lib
			${vtk_ROOT}/lib/vtkRenderingContextOpenGL-8.0.lib
			${vtk_ROOT}/lib/vtkRenderingCore-8.0.lib
			${vtk_ROOT}/lib/vtkRenderingFreeType-8.0.lib
			${vtk_ROOT}/lib/vtkRenderingGL2PS-8.0.lib
			${vtk_ROOT}/lib/vtkRenderingImage-8.0.lib
			${vtk_ROOT}/lib/vtkRenderingLabel-8.0.lib
			${vtk_ROOT}/lib/vtkRenderingLIC-8.0.lib
			${vtk_ROOT}/lib/vtkRenderingLOD-8.0.lib
			${vtk_ROOT}/lib/vtkRenderingOpenGL-8.0.lib
			${vtk_ROOT}/lib/vtkRenderingVolume-8.0.lib
			${vtk_ROOT}/lib/vtkRenderingVolumeOpenGL-8.0.lib
			${vtk_ROOT}/lib/vtksqlite-8.0.lib
			${vtk_ROOT}/lib/vtksys-8.0.lib
			${vtk_ROOT}/lib/vtktiff-8.0.lib
			${vtk_ROOT}/lib/vtkverdict-8.0.lib
			${vtk_ROOT}/lib/vtkViewsContext2D-8.0.lib
			${vtk_ROOT}/lib/vtkViewsCore-8.0.lib
			${vtk_ROOT}/lib/vtkViewsInfovis-8.0.lib
			${vtk_ROOT}/lib/vtkzlib-8.0.lib
			glu32.lib opengl32.lib
			#msvcrt.lib
			#ucrt.lib
			#vcruntime.lib
		)
	endif()
	if(USE_LIDAR AND USE_CERES)
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
endif()

add_dependencies(${PROJECT_NAME} ${libBase} ${libMat} ${libAmb} ${libLSQ} )